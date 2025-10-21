#include "boost/random/cauchy_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/math/distributions.hpp"
#include <boost/random.hpp>
#include <iostream>
#include "omp.h"
#include "first.h"
#include <algorithm>

boost::random::normal_distribution<double> get_normal(int n, double h1, double h2){
    return boost::random::normal_distribution<double>(-h1 / (h2 + sqrt(n)), 1 / (1 + h2/sqrt(n) ) );
}

boost::random::cauchy_distribution<double> get_cauchy(int n, double h1, double h2){
    return boost::random::cauchy_distribution<double>(-h1 / (h2 + sqrt(n)), 1 / (1 + h2/sqrt(n) ) );
}




void print_vector(std::vector<double> V){
    for(int i=0; i<V.size(); ++i){
        std::cout << V[i] << " ";
    }
    std::cout << "\n";
}

std::vector<double> read_numbers_from_file(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<double> numbers;
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    double number;
    while (file >> number) {
        numbers.push_back(number);
    }
    
    file.close();
    return numbers;
}

std::vector<double> read_integrals(DistributionType d_type){
    switch (d_type)
    {
    case DistributionType::NORMAL:
        return read_numbers_from_file("/home/lolikion/Документы/study/нир5сем/code/ex/precomputed_integrals/normal.txt");
        break;
    case DistributionType::CAUCHY:
        return read_numbers_from_file("/home/lolikion/Документы/study/нир5сем/code/ex/precomputed_integrals/cauchy.txt");
        break;
    
    default:
        std::cout << "WRONG INDEX FOR READ_INTEGRALS\n";
        exit(1);
    }
}


double compute_ad(std::function<double(double)> g, double *X, double *Y, int sample_size, std::function<double(double)> F1){
    std::vector<double> Z;
    Z.push_back(0);
    for(int i=0; i<sample_size; ++i){
        Z.push_back(X[i]);
        Z.push_back(Y[i]);
    }
    std::sort(Z.begin()+1, Z.end());
    int N = 2 * sample_size;

    double res = 0;

    # pragma omp parallel for reduction(+: res)
    for(int k=1; k<=N; ++k){
        res += (2 * k - 1) * ( log(F1(Z[k])) + log(1 - F1(Z[N - k + 1])) );
    }

    return -N - res / N;
}



double compute_wmw(std::function<double(double)> g, double *X, double *Y, int sample_size, std::function<double(double)> F1){
    double res = 0;

    #pragma omp parallel for reduction(+:res)
    for(int i=0; i<sample_size; ++i){
        for(int j=0; j<sample_size; ++j){
            res += (Y[j] <= X[i]);
        }
    }

    return res;
}


double compute_etest(std::function<double(double)> g, double *X, double *Y, int sample_size, std::function<double(double)> F1){
    double phi_a = 0, phi_b = 0, phi_ab = 0;

    #pragma omp parallel for reduction(+:phi_a, phi_b, phi_ab)
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



