#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <random>
#include <omp.h>

template <typename T>
void sample(T dist, int sample_size, double *X){
    std::random_device rd;
    boost::random::mt19937 gen(rd());
    for(int i=0; i<sample_size; ++i){
        X[i] = dist(gen);
    }
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


    #pragma omp parallel reduction(+:cnt_reject)
    {   
        double *X = new double[n];
        double *Y = new double[n];
        double etest_val;

        #pragma omp for 
        for(int i=0; i<N; ++i){
            sample(d1, n, X);
            sample(d2, n, Y);
            cnt_reject += ( (n * compute_etest(g, X, Y, n)) >= crit_val );
            
        }
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
    
    double *Z_ = new double[2*n];
    sample(d1, 2*n, Z_);
    std::vector<double> Z;
    for(int i=0; i<2*n; ++i){
        Z.push_back(Z_[i]);
    }
    delete[] Z_;
    
    std::vector<double> test_vals(M);
    
    #pragma omp parallel
    {

        double *X = new double[n];
        double *Y = new double[n];

        #pragma omp for
        for(int i=0; i<M; ++i){
            random_split_direct(Z, n, X, Y);
            double etest_val = compute_etest(g, X, Y, n);
            test_vals[i] = n * etest_val;
        }

        delete[] X;
        delete[] Y;
    }

    return quantile(test_vals, 1 - alpha);
}


// 0 for normal
// 1 for cauchy
std::vector<double> read_integrals(int index){
    switch (index)
    {
    case 0:
        break;
    case 1:
        break;
    
    default:
        std::cout << "WRONG INDEX FOR READ_INTEGRALS\n";
        exit(1);
    }
}



void experiment_step(std::function<double(double)> g,
                     std::function<double(double)> d2_g,
                    double h1, double h2, double alpha,
                    int N, int M, std::vector<int> sample_sizes,
                    std::vector<double> integrals,
                    std::vector<double> crit_vals){

    double J1 = integrals[0];
    double J2 = integrals[1];
    double J3 = integrals[2];
    double J1_star = integrals[3] * h1 * h1;
    double J2_star = integrals[4] * h2 * h2;

    double b1 = sqrt(abs(J1_star));
    double b2 = sqrt(abs(J2_star));
    double b = sqrt( b1 * b1 + b2 * b2);
    
    double a = pow( J2 + J1 * J1 - 2 * J3, 0.25);

}


template <typename T>
void run_experiment(std::function<double(double)> g,
                    std::function<double(double)> d2_g,
                    double h1, std::vector<double> h2_vals,
                    double alpha, int N, int M,
                    std::vector<int> sample_sizes,
                    T dist_template, int exp_index)
{
    auto d1 = ...;
    auto f = ...;

    // J1 J2 J3 J1_ J2_
    std::vector<double> integrals = read_integrals(exp_index);
    
    std::vector<double> crit_vals;
    for(int n : sample_sizes){
        double cv = compute_crit_val(n, M, alpha, d1, g);
        crit_vals.push_back(cv);
    }

    
    for(double h2 : h2_vals){
        auto [e_pow, a_pow] = ...;
    }


}