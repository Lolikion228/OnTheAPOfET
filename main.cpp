
#include "boost/random/cauchy_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <iostream>
#include "first.h"
#include <random>

double g(double x){
    return log(1 + x*x);
}


void print_vector(std::vector<double> V){
    for(int i=0; i<V.size(); ++i){
        std::cout << V[i] << " ";
    }
    std::cout << "\n";
}

int main() {
    std::random_device rd;
    boost::random::mt19937 gen(rd());

    
    // int sample_size = 8;
    // double *X = new double[sample_size]{1, 2, 3, 4, 4, 2, 2, 3};
    // double *Y = new double[sample_size]{9, 9, 9, 9, 1, 3, 8, 8};
    // double etest_val = compute_etest(g, X, Y, sample_size);
    // std::cout << etest_val << "\n";
    

    // std::vector<double> X{1, 2, 3, 4, 4, 2, 2, 3, 11, 12, 17, 6, 20, 14};
    // std::cout << quantile(X, 0.87) << "\n";

    // std::cout << compute_asymptotic_power(0.02, 0.1, 0.3) << "\n";


    // std::vector<double> Z{10,20,30,40,50,60};
    // int sample_size = Z.size()/2;
    // double *X = new double[sample_size];
    // double *Y = new double[sample_size];
    // print_vector(Z);
    // random_split_direct(Z,sample_size , X, Y, gen);
    // print_sample(X, sample_size);
    // print_sample(Y, sample_size);

    // int n = 300;
    // int M = 500;
    // double alpha = 0.03;
    // boost::random::normal_distribution<> d1(0,1);
    // double cv = compute_crit_val(n, M, alpha, d1, g);
    // std::cout << cv << "\n";

    // int n = 100;
    // int N = 5000;
    // double cv = 55;
    // boost::random::normal_distribution<> d1(0, 1);
    // boost::random::normal_distribution<> d2(0.2 ,5);
    // double ep = compute_empirical_power(n, N, cv, d1, d2, g, gen);
    // std::cout << ep << "\n";

}