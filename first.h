#include "boost/random/cauchy_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <iostream>

template <typename T>
double *sample(boost::random::mt19937 gen, T dist, int sample_size){
    double *sample = new double[sample_size];
    for(int i=0; i<sample_size; ++i){
        sample[i] = dist(gen);
    }
    return sample;
}

void print_sample(double *sample, int sample_size);