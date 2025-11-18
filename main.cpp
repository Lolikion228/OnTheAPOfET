#include "boost/random/cauchy_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <iostream>
#include "first.h"
#include <random>
#include <chrono>
#include "omp.h"

// 117 
// 451
// 68
boost::random::mt19937 gen(68);

double g(double x){
    return log(1 + x*x);
}

// double d2_g(double x){
//     return 2 * (1 - x*x) / pow((1 + x * x),2);
// }

double g2(double x){
    return x>=0 ? x : -x;
}

// double d2_g2(double x){
//     return 0;
// }

void testing(){
    // int sample_size = 10;
    // double *X = new double[sample_size]{14, -10, 1, 2, 3, 4, 4, 2, 2, 3};
    // double *Y = new double[sample_size]{1, 20, 9, 9, 9, 9, 1, 3, 8, 8};
    // double etest_val = compute_etest(g, X, Y, sample_size);
    // std::cout << etest_val << "\n";
    

    // std::vector<double> X{2, 2, 3, 11, 12, 17, 6, 20, 14};
    // std::cout << quantile(X, 0.39) << "\n";

    // std::cout << compute_asymptotic_power(0.14, 0.7, 0.1) << "\n";


    // std::vector<double> Z{10,20,30,40,50,60};
    // int sample_size = Z.size()/2;
    // double *X = new double[sample_size];
    // double *Y = new double[sample_size];
    // print_vector(Z);
    // random_split_direct(Z,sample_size , X, Y);
    // print_sample(X, sample_size);
    // print_sample(Y, sample_size);
    // delete[] X;
    // delete[] Y;

    int n = 400;
    int M = 5000;
    double alpha = 0.12;
    boost::random::cauchy_distribution<> d1(0,1);
    // boost::random::cauchy_distribution<> d2(-3/sqrt(n), 1);
    double cv1 = compute_crit_val(n, M, alpha, d1,
        [](double* X, double *Y, int n ){return compute_etest(g2, X, Y, n);}, gen);
    double cv2 = compute_crit_val_v2(n, M, alpha, d1,
        [](double* X, double *Y, int n ){return compute_etest(g2, X, Y, n);}, gen);
    std::cout << "cv1 = " << cv1 << "\n";
    std::cout << "cv2 = " << cv2 << "\n";

    int N = 5000;
    double *X = new double[2 * n];
    double cnt1 = 0;
    double cnt2 = 0; 
    for(int i=0; i<N; ++i){
        sample(d1, 2*n, X, gen);
        if(n * compute_etest(g2, X, X + n, n) >= cv1){
            ++cnt1;
        }
        if(n * compute_etest(g2, X, X + n, n) >= cv2){
            ++cnt2;
        }
    }
    std::cout << "(1) type 1 err: " << cnt1 / N << "\n";
    std::cout << "(2) type 1 err: " << cnt2 / N << "\n";
    

    // std::cout << "(1) EP = " << compute_empirical_power(n, N, cv1, d1, d2,
    // [](double* X, double *Y, int n ){return compute_etest(g2, X, Y, n);}, gen ) << "\n";
    // std::cout << "(2) EP = " << compute_empirical_power(n, N, cv2, d1, d2,
    // [](double* X, double *Y, int n ){return compute_etest(g2, X, Y, n);}, gen ) << "\n";

}


template <typename T>
void ex_tmp(std::vector<double> h1_vals, double h2, double alpha, int N, int M,
            std::vector<int>sample_sizes, std::vector<double> integrals, T get_dist, const char name[])
{

    // std::cout<< "ET_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals,
    //      [](double* X, double *Y, int n ){return compute_etest(g, X, Y, n);}, false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    std::cout<< "HT_" << name << "\n";
    run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals,
    [](double* X, double *Y, int n ){return compute_etest(g2, X, Y, n);}, false, gen);
    std::cout << "================================\n";
    std::cout << "================================\n\n";

    // std::cout<< "WMW_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_wmw, false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "AD_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_ad, false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "KS_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_ks,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "CM_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_cm,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";
}


template <typename T>
void ex_tmp(double h1, std::vector<double> h2_vals,  double alpha, int N, int M,
            std::vector<int>sample_sizes, std::vector<double> integrals, T get_dist, const char name[])
{

    // std::cout<< "ET_" << name << "\n";
    // run_experiment(  h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals,
    //      [](double* X, double *Y, int n ){return compute_etest(g, X, Y, n);}, false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    std::cout<< "HT_" << name << "\n";
    run_experiment(h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals,
    [](double* X, double *Y, int n ){return compute_etest(g2, X, Y, n);}, false, gen);
    std::cout << "================================\n";
    std::cout << "================================\n\n";

    // std::cout<< "WMW_" << name << "\n";
    // run_experiment(h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals, compute_wmw, false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "AD_" << name << "\n";
    // run_experiment(h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals, compute_ad, false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "KS_" << name << "\n";
    // run_experiment(h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals, compute_ks,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "CM_" << name << "\n";
    // run_experiment(h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals, compute_cm,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";
}


// normal with h1=0
void ex1(){
    std::vector<double>  h2_vals{1,  3, 5};
    double h1 = 0.0;
    double alpha = 0.05;
    std::vector<int> sample_sizes{100, 400, 900};
    int N = 5000;
    int M = 5000;
    std::vector<double> integrals = read_integrals(DistributionType::NORMAL);
    ex_tmp(h1, h2_vals, alpha, N, M, sample_sizes, integrals, get_normal, "NORMAL");
}


// normal with h2=0
void ex2(){
    std::vector<double>  h1_vals{1,  3,  5};
    double h2 = 0.0;
    double alpha = 0.05;
    std::vector<int> sample_sizes{100, 400, 900};
    int N = 5000;
    int M = 5000;
    std::vector<double> integrals = read_integrals(DistributionType::NORMAL);
    ex_tmp(h1_vals, h2, alpha, N, M, sample_sizes, integrals, get_normal, "NORMAL");
}


// cauchy with h2=0
void ex3(){
    std::vector<double>  h1_vals{1, 3, 5, 7, 9};
    double h2 = 0.0;
    double alpha = 0.05;
    std::vector<int> sample_sizes{100, 400, 900};
    int N = 5000;
    int M = 5000;
    std::vector<double> integrals = read_integrals(DistributionType::CAUCHY);
    ex_tmp(h1_vals, h2, alpha, N, M, sample_sizes, integrals, get_cauchy, "CAUCHY");
}


// cauchy with h1=0
void ex4(){
    std::vector<double>  h2_vals{1, 3, 5, 7, 9};
    double h1 = 0.0;
    double alpha = 0.05;
    std::vector<int> sample_sizes{100, 400, 900};
    int N = 5000;
    int M = 5000;
    std::vector<double> integrals = read_integrals(DistributionType::CAUCHY);
    ex_tmp(h1, h2_vals, alpha, N, M, sample_sizes, integrals, get_cauchy, "CAUCHY");
}


// ДОМНОЖЕНИЕ НА n ОСТАВИТЬ ТОЛЬКО ДЛЯ ET???
// заменить метод перестановок
int main() {
    testing();
    // ex1();
    // ex2();
    // ex3();
    // ex4();
}