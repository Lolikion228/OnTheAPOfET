
#include "boost/random/cauchy_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <iostream>
#include "first.h"

double g(double x){
    return log(1 + x*x);
}


int main() {
    const int random_seed = 112;
    boost::random::mt19937 gen(random_seed);

    
    // int sample_size = 8;
    // double *X = new double[sample_size]{1, 2, 3, 4, 4, 2, 2, 3};
    // double *Y = new double[sample_size]{9, 9, 9, 9, 1, 3, 8, 8};
    // double etest_val = compute_etest(g, X, Y, sample_size);
    // std::cout << etest_val << "\n";
    

    // std::vector<double> X{1, 2, 3, 4, 4, 2, 2, 3, 11, 12, 17, 6, 20, 14};
    // std::cout << quantile(X, 0.87) << "\n";

    std::cout << compute_asymptotic_power(0.02, 0.1, 0.3) << "\n";
}