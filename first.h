#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <random>
#include <omp.h>

template <typename T>
double *sample(T dist, int sample_size){
    std::random_device rd;
    boost::random::mt19937 gen(rd());
    double *sample = new double[sample_size];
    for(int i=0; i<sample_size; ++i){
        sample[i] = dist(gen);
    }
    return sample;
}


//
void print_sample(double *sample, int sample_size);


//
double compute_etest(std::function<double(double)> g, double *X, double *Y, int sample_size);


//
double compute_asymptotic_power(double alpha, double b, double a);

//
template <typename T>
double compute_empirical_power(int n, int N, double crit_val, T d1, T d2, std::function<double(double)> g, boost::random::mt19937 gen){
    double cnt_reject = 0;


    #pragma omp parallel for
    for(int i=0; i<N; ++i){
        double *X = sample(d1, n);
        double *Y = sample(d2, n);
        double etest_val = compute_etest(g, X, Y, n);
        cnt_reject += ( (n * etest_val) >= crit_val );
        delete[] X;
        delete[] Y;
    }

    
    return cnt_reject / N;
}


//
template<typename T>
double quantile(const std::vector<T>& data, double probability) {
    if (data.empty()) return 0.0;
    if (probability < 0.0 || probability > 1.0) {
        throw std::invalid_argument("Probability must be in [0, 1]");
    }
    
    std::vector<T> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());
    
    if (probability == 0.0) return sorted_data.front();
    if (probability == 1.0) return sorted_data.back();
    
    double index = probability * (sorted_data.size() - 1);
    size_t lower_index = static_cast<size_t>(std::floor(index));
    size_t upper_index = static_cast<size_t>(std::ceil(index));
    
    if (lower_index == upper_index) {
        return sorted_data[lower_index];
    }
    
    // Линейная интерполяция
    double weight = index - lower_index;
    return (1 - weight) * sorted_data[lower_index] + weight * sorted_data[upper_index];
}


//
template<typename T>
void random_split_direct(std::vector<T> Z, size_t n, double *X, double *Y) {
    std::random_device rd;
    boost::random::mt19937 gen(rd());
    std::shuffle(Z.begin(), Z.end(), gen);

   
    for(int i=0; i<n; ++i){
        X[i] = Z[i];
        Y[i] = Z[i+n];
    }

}


//
template <typename T>
double compute_crit_val(int n, int M, double alpha, T d1, std::function<double(double)> g){
    
    double *Z_ = sample(d1, 2*n);
    std::vector<double> Z;
    for(int i=0; i<2*n; ++i){
        Z.push_back(Z_[i]);
    }
    
    std::vector<double> test_vals;
    double *X = new double[n];
    double *Y = new double[n];

    #pragma omp parallel for
    for(int i=0; i<M; ++i){
        random_split_direct(Z, n, X, Y);
        double etest_val = compute_etest(g, X, Y, n);
        test_vals.push_back(n * etest_val);
    }
    
    delete[] X;
    delete[] Y;

    return quantile(test_vals, 1 - alpha);
}

