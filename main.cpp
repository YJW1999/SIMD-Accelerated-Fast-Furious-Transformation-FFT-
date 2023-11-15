#include <complex>
#include <iostream>
#include <valarray>
#include <vector>
#include <arm_neon.h>
#include <chrono>
#include <thread>

const double PI = 3.141592653589793238460;

class CNumber {
public:
    float real;
    float imag;
    CNumber() {}

    CNumber(float real) {
        this->real = real;
        this->imag = 0.0;
    }
    CNumber(float real, float imag) {
        this->real = real;
        this->imag = imag;
    }
    CNumber(const CNumber& other) {
        this->real = other.real;
        this->imag = other.imag;
    }
    CNumber(CNumber&& other) {
        this->real = other.real;
        this->imag = other.imag;
    }
    CNumber operator*(const CNumber& other) const {
        CNumber result(
            this->real * other.real - this->imag * other.imag,
            this->real * other.imag + this->imag * other.real
        );
        return result;
    }
    CNumber operator+(const CNumber& other) {
        return CNumber(this->real + other.real, this->imag + other.imag);
    }
    CNumber operator-(const CNumber& other) {
        return CNumber(this->real - other.real, this->imag - other.imag);
    }
    CNumber operator=(const CNumber& other) {
        if (this != &other) {
            this->real = other.real;
            this->imag = other.imag;
        }
        return *this;
    }
    friend std::ostream& operator<<(std::ostream& os, const CNumber& obj) {
        os << "Real: " << obj.real << " Imag: " << obj.imag;
        return os;
    }
};

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
typedef std::vector<CNumber> Cvector;

void fft(Cvector& x) {
    const unsigned int N = (unsigned int)x.size();
    if (N <= 1) return;
    // make sure it is a power of 2
    unsigned int log2N = (log2(N));
    if (N != (1 << log2N)) {
        throw std::invalid_argument("N must be a power of 2");
    }
    //
    Cvector temp(N);
    for (unsigned int i = 0; i < N; ++i) {
        unsigned int j = 0;
        for (unsigned int bit = 0; bit < log2N; ++bit) {
            j <<= 1;
            j |= (i >> bit) & 1;
        }
        temp[i] = x[j];
    }
    x = temp;
    // fft
    for (unsigned int s = 1; s <= log2N; s++) {
        unsigned int m = 1 << s;
        Complex wm_temp = std::polar(1.0, -2 * PI / m);
        CNumber wm(wm_temp.real(), wm_temp.imag());
        for (unsigned int k = 0; k < N; k += m) {
            CNumber w = 1.0;
            for (unsigned int j = 0; j < m / 2; j++) {
                CNumber t = w * x[k + j + m / 2];
                //std::cout<< t << std::endl;
                CNumber u = x[k + j];
                x[k + j] = u + t;
                x[k + j + m / 2] = u - t;
                w = wm * w;
            }
        }
    }
}

// A function for conducting multithread
void update_x(const unsigned int log2N, std::vector<float> &real_vector, std::vector<float> &imag_vector, Cvector &x, unsigned int start, unsigned int end){
    for (unsigned int i = start; i < end; ++i) {
        unsigned int j = 0;
        for (unsigned int bit = 0; bit < log2N; bit++) {
            j <<= 1;
            j |= (i >> bit) & 1;
        }
        //temp[i] = x[j];
        real_vector[i] = x[j].real;
        imag_vector[i] = x[j].imag;
    }
}

void fft_simd(Cvector& x) {
    const unsigned int N = (unsigned int)x.size();
    if (N <= 1) return;
    // make sure it is a power of 2
    unsigned int log2N = (log2(N));
    if (N != (1 << log2N)) {
        throw std::invalid_argument("N must be a power of 2");
    }
    //
    Cvector temp(N);
    std::vector<float> real_vector(N);
    std::vector<float> imag_vector(N);
    std::vector<std::thread> t_pool;
    for(int i=0; i<8; i++){
        t_pool.emplace_back(std::thread(update_x, log2N, std::ref(real_vector), std::ref(imag_vector), std::ref(x), N*i/8, N*(i+1)/8));
    }
    for(auto &thread : t_pool){
        thread.join();
    }
    x = temp;
    // fft
    for (unsigned int s = 1; s <= log2N; s++) {
        unsigned int m = 1 << s;
        CNumber wm(cos(-2 * PI / m), sin(-2 * PI / m));
        for (unsigned int k = 0; k < N; k += m) {
            CNumber w = 1.0;
            CNumber w1 = w*wm;
            CNumber w2 = w1*wm;
            CNumber w3 = w2*wm;
            float32x4_t w_real = {w.real, w1.real, w2.real, w3.real};
            float32x4_t w_imag = {w.imag, w1.imag, w2.imag, w3.imag};
            unsigned int j = 0;
            unsigned int num_simd_elements = 4;
            //std::cout << (m/2) << std::endl;
            for(; (m / 2)-j>=num_simd_elements; j+=num_simd_elements) {
                //CNumber t = w * x[k + j + m / 2];
                float* t_real_ptr = &real_vector[k+j+m/2];
                float32x4_t xkjm2_real = vld1q_f32(t_real_ptr);
                float* t_imag_ptr = &imag_vector[k+j+m/2];
                float32x4_t xkjm2_imag = vld1q_f32(t_imag_ptr);
                float32x4_t t_real = vmulq_f32(w_real, xkjm2_real);
                float32x4_t t_imag = vmulq_f32(w_real, xkjm2_imag);
                float32x4_t t_real_p2 = vmulq_f32(w_imag, xkjm2_imag);
                float32x4_t t_imag_p2 = vmulq_f32(w_imag, xkjm2_real);
                t_real = vsubq_f32(t_real, t_real_p2);
                t_imag = vaddq_f32(t_imag, t_imag_p2);
                //CNumber u = x[k + j];
                float* u_real_ptr = &real_vector[k+j];
                float32x4_t u_real = vld1q_f32(u_real_ptr);
                float* u_imag_ptr = &imag_vector[k+j];
                float32x4_t u_imag = vld1q_f32(u_imag_ptr);
                //x[k + j] = u + t;
                float32x4_t result_real = vaddq_f32(t_real, u_real);
                float32x4_t result_imag = vaddq_f32(t_imag, u_imag);
                float* xkj_real_ptr = &real_vector[k+j];
                vst1q_f32(xkj_real_ptr, result_real);
                float* xkj_imag_ptr = &imag_vector[k+j];
                vst1q_f32(xkj_imag_ptr, result_imag);
                //x[k + j + m / 2] = u - t;
                result_real = vsubq_f32(u_real, t_real);
                result_imag = vsubq_f32(u_imag, t_imag);
                float* xkjm2_real_ptr = &real_vector[k+j+m/2];
                vst1q_f32(xkjm2_real_ptr, result_real);
                float* xkjm2_imag_ptr = &imag_vector[k+j+m/2];
                vst1q_f32(xkjm2_imag_ptr, result_imag);
                //w = wm * w;
                w = w3*wm;
                w1 = w*wm;
                w2 = w1*wm;
                w3 = w2*wm;
                w_real = {w.real, w1.real, w2.real, w3.real};
                w_imag = {w.imag, w1.imag, w2.imag, w3.imag};
            }
            for(; j<m/2; ++j){
                //CNumber t = w * x[k + j + m / 2];
                float t_real = w.real * real_vector[k + j + m / 2] - w.imag * imag_vector[k + j + m / 2];
                float t_imag = w.real * imag_vector[k + j + m / 2] + w.imag * real_vector[k + j + m / 2];
                //std::cout<< t_real << " " << t_imag << std::endl;
                //CNumber u = x[k + j];
                float u_real = real_vector[k+j];
                float u_imag = imag_vector[k+j];
                //x[k + j] = u + t;
                real_vector[k+j] = u_real + t_real;
                imag_vector[k+j] = u_imag + t_imag;
                //x[k + j + m / 2] = u - t;
                real_vector[k+j+m/2] = u_real - t_real;
                imag_vector[k + j + m / 2] = u_imag - t_imag;
                //w = wm * w;
                w = wm * w;
            }
        }
    }

    /*for(int i=0; i<N; ++i){
        x[i] = CNumber(real_vector[i], imag_vector[i]);
    }*/
}

int main()
{
    //const Complex test[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,10.0, 11.0,12.0,13.0,14.0,15.0,16.0 };
    //CArray data(test, 16);
    Cvector data1;
    Cvector data2;
    //4194304
    for(int i=0; i< 33554432; i++){
        float data = (float)(rand()%1000);
        data1.emplace_back(data);
        data2.emplace_back(data);
    }
    
    // fft test
    std::cout << "fft" << std::endl;
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    fft(data1);
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> period = end_time - start_time;
    std::cout << period.count() << std::endl;
    
    // simd version of fft test
    std::cout << "fft_simd" << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
    fft_simd(data2);
    end_time = std::chrono::high_resolution_clock::now();
    period = end_time - start_time;
    std::cout << period.count() << std::endl;
    
    for(int i=0; i< 200; i++){
        //std::cout << data1[i] << std::endl;
        //std::cout << data2[i] << std::endl;
    }

    return 0;
}
