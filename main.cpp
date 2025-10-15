#include "boost/random/cauchy_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <iostream>
#include "first.h"
#include <random>
#include <chrono>
#include "omp.h"



double g(double x){
    return log(1 + x*x);
}

double d2_g(double x){
    return 2 * (1 - x*x) / pow((1 + x * x),2);
}


int fat(){
    int x = 2;
    for(int i=0; i<10000; ++i){
        for(int j=0; j<10000; ++j){
           x = (x+1) % 2;
        }
    }
    return x;
}

int main() {
    std::random_device rd;
    boost::random::mt19937 gen(rd());
    auto start = std::chrono::high_resolution_clock::now();
    

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

    // int n = 20;
    // int M = 10;
    // double alpha = 0.01;
    // boost::random::normal_distribution<> d1(0,1);
    // double cv = compute_crit_val(n, M, alpha, d1, g);
    // std::cout << cv << "\n";


    
    // const int N = 100;
    // std::vector<int> data(N, 1);
    // #pragma omp parallel for
    // for(int i=0; i<N; ++i){
    //     data[i] = fat();
    // }
    // std :: cout<<data[0];


    // int n = 20;
    // int N = 20;
    // double cv = 51;
    // boost::random::normal_distribution<> d1(0, 1);
    // boost::random::normal_distribution<> d2(0.2, 5);
    // double ep = compute_empirical_power(n, N, cv, d1, d2, g);
    // std::cout << ep << "\n";

    double h1 = 0.0;
    std::vector<double> h2_vals{1, 2, 23};
    double alpha = 0.05;
    int N = 100;
    int M = 100;
    std::vector<int> sample_sizes{50, 100};

    run_experiment(g, d2_g, h1, h2_vals, alpha, N, M, sample_sizes);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Время выполнения: " << duration.count() / 1000000.0 << " секунд" << std::endl;

}