#include "boost/random/cauchy_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/math/distributions.hpp"
#include <boost/random.hpp>
#include <iostream>
#include "omp.h"



double compute_etest(std::function<double(double)> g, double *X, double *Y, int sample_size){
    double phi_a = 0, phi_b = 0, phi_ab = 0;

    #pragma omp parallel for
    for(int i=0; i<sample_size; ++i){
        for(int j=0; j<=i; ++j){
            phi_ab += g(X[i] - Y[j]);
        }
        for(int j=i+1; j<sample_size; ++j){
            phi_a += g(X[i] - X[j]);
            phi_b += g(Y[i] - Y[j]);
            phi_ab += g(X[i] - Y[j]);
        }

    }

    return (phi_ab - phi_a - phi_b) / (sample_size * sample_size);
}


double compute_asymptotic_power(double alpha, double b, double a){
    boost::math::normal_distribution<> dist(0, 1);
    double z = boost::math::quantile(dist, 1 - alpha / 2);
    return 1 - boost::math::cdf(dist, z - b / a) + boost::math::cdf(dist, -z - b / a);
}


void print_sample(double *sample, int sample_size){
    for(int i=0; i<sample_size; ++i){
        std::cout << sample[i] << " ";
    }
    std::cout << "\n";
}
