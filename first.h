#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <random>
#include <omp.h>
#include <fstream>


#include <chrono>


class Timer {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    
public:
    Timer() : start(std::chrono::high_resolution_clock::now()) {}
    
    ~Timer() {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Время выполнения: " << duration.count()/ 1000000.0 << " сек" << std::endl;
    }
};


template <typename T>
void sample(T dist, int sample_size, double *X){
    std::random_device rd;
    boost::random::mt19937 gen(rd());
    for(int i=0; i<sample_size; ++i){
        X[i] = dist(gen);
    }
}


void print_vector(std::vector<double> V);


//
void print_sample(double *sample, int sample_size);

double compute_edf(std::vector<double>& ordered_sample, double t);


double compute_cm(double *X, double *Y, int sample_size);

double compute_etest(std::function<double(double)> g, double *X, double *Y, int sample_size);

double compute_wmw(double *X, double *Y, int sample_size);

double compute_ad(double *X, double *Y, int sample_size);

double compute_ks(double *X, double *Y, int sample_size);
//

double compute_asymptotic_power(double alpha, double b, double a);

//
template <typename T>
double compute_empirical_power(int n, int N, double crit_val, T d1, T d2,
                               std::function<double(double*, double*, int)> compute_test)
{
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
            cnt_reject += ( (n * compute_test(X, Y, n)) >= crit_val );
            
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
double compute_crit_val(int n, int M, double alpha, T d1,
                        std::function<double(double*, double*, int)> compute_test)
{
    
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
            double test_val = compute_test(X, Y, n);
            test_vals[i] = n * test_val;
        }

        delete[] X;
        delete[] Y;
    }

    return quantile(test_vals, 1 - alpha);
}


enum class DistributionType {
    NORMAL = 0,
    CAUCHY = 1
};

std::vector<double> read_numbers_from_file(const std::string& filename);


std::vector<double> read_integrals(DistributionType d_type);


template <typename T>
std::pair<std::vector<double>, double> experiment_step(
                    double h1, double h2, double alpha,
                    int N, int M, std::vector<int> sample_sizes,
                    std::vector<double> integrals,
                    std::vector<double> crit_vals,
                    T dist_template,
                    std::function<double(double*, double*, int)> compute_test)
{

    auto d1 = dist_template(1,0,0);

    double J1 = integrals[0];
    double J2 = integrals[1];
    double J3 = integrals[2];
    double J1_star = integrals[3] * h1 * h1;
    double J2_star = integrals[4] * h2 * h2;

    double b1 = sqrt(abs(J1_star));
    double b2 = sqrt(abs(J2_star));
    
    double a = pow( J2 + J1 * J1 - 2 * J3, 0.25);
    double b = sqrt( b1 * b1 + b2 * b2);
    
    std::vector<double> emp_powers;

    for (int i=0; i<sample_sizes.size(); ++i){
        Timer t1;
        int n = sample_sizes[i];
        std::cout << "n = " << n << "  ||  ";
        auto d2_n = dist_template(n, h1, h2);
        double e_pow = compute_empirical_power(n, N, crit_vals[i], d1, d2_n, compute_test);
        emp_powers.push_back(e_pow);
    }

    double a_pow = compute_asymptotic_power(alpha, b, a);

    return {emp_powers, a_pow};

}



boost::random::normal_distribution<double> get_normal(int n, double h1, double h2);


boost::random::cauchy_distribution<double> get_cauchy(int n, double h1, double h2);

template <typename T>
void run_experiment(double h1, std::vector<double> h2_vals,
                    double alpha, int N, int M,
                    std::vector<int> sample_sizes,
                    T dist_template, std::vector<double> integrals,
                    std::function<double(double*, double*, int)> compute_test,
                    bool compute_AP)
{   

    std::cout << "N = " << N << "\n"; 
    std::cout << "M = " << M << "\n"; 
    std::cout << "alpha = " << alpha << "\n"; 
    std::cout << "h1 = " << h1 << "\n\n"; 
    auto d1 = dist_template(1,0,0);

    std::vector<double> crit_vals;
    std::cout << "computing crit_vals...\n";
    for(int n : sample_sizes){
        Timer t1;
        std::cout << "n = " << n << "  ||  ";
        double cv = compute_crit_val(n, M, alpha, d1, compute_test);
        crit_vals.push_back(cv);
    }
    std::cout << "\n";
    
    for(double h2 : h2_vals){
        {
        Timer t1;
        std::cout << "h2 = " << h2 << "\n";
        auto [e_pow, a_pow] = experiment_step(h1, h2,
            alpha, N, M, sample_sizes, integrals, crit_vals, dist_template, compute_test);
        std::cout << "emp_powers = ";
        print_vector(e_pow);
        if(compute_AP){
            std::cout << "asp_power = " << a_pow << "\n";
        }
        }
        std::cout << "\n\n";  
    }
}




template <typename T>
void run_experiment(std::vector<double> h1_vals, double h2, 
                    double alpha, int N, int M,
                    std::vector<int> sample_sizes,
                    T dist_template, std::vector<double> integrals,
                    std::function<double(double*, double*, int)> compute_test,
                    bool compute_AP)
{   

    std::cout << "N = " << N << "\n"; 
    std::cout << "M = " << M << "\n"; 
    std::cout << "alpha = " << alpha << "\n"; 
    std::cout << "h2 = " << h2 << "\n\n"; 
    auto d1 = dist_template(1,0,0);

    std::vector<double> crit_vals;
    std::cout << "computing crit_vals...\n";
    for(int n : sample_sizes){
        Timer t1;
        std::cout << "n = " << n << "  ||  ";
        double cv = compute_crit_val(n, M, alpha, d1, compute_test);
        crit_vals.push_back(cv);
    }
    std::cout << "\n";
    
    for(double h1 : h1_vals){
        {
        Timer t1;
        std::cout << "h1 = " << h1 << "\n";
        auto [e_pow, a_pow] = experiment_step(h1, h2,
            alpha, N, M, sample_sizes, integrals, crit_vals, dist_template, compute_test);
        std::cout << "emp_powers = ";
        print_vector(e_pow);
        if (compute_AP){
            std::cout << "asp_power = " << a_pow << "\n";
        }
        }
        std::cout << "\n\n";  
    }
}