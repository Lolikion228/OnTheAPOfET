#include "boost/random/cauchy_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/math/distributions.hpp"
#include <boost/random.hpp>
#include <iostream>
#include "omp.h"
#include "first.h"


double compute_etest(std::function<double(double)> g, double *X, double *Y, int sample_size){
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



void run_experiment(std::function<double(double)> g,
                    std::function<double(double)> d2_g,
                    double h1, std::vector<double> h2_vals,
                    double alpha, int N, int M,
                    std::vector<int> sample_sizes)
{   


    // HARDCODED!!!
    auto dist_template = get_normal; 

    // HARDCODED!!!
    // J1 J2 J3 J1_ J2_
    std::vector<double> integrals = read_integrals(DistributionType::NORMAL);


    auto d1 = dist_template(1,0,0);

    std::vector<double> crit_vals;
    for(int n : sample_sizes){
        double cv = compute_crit_val(n, M, alpha, d1, g);
        crit_vals.push_back(cv);
    }

    
    for(double h2 : h2_vals){
        auto [e_pow, a_pow] = experiment_step(g, d2_g, h1, h2,
            alpha, N, M, sample_sizes, integrals, crit_vals, dist_template);
        std::cout << "h2 = " << h2 << "\n";
        std::cout << "emp_powers = ";
        print_vector(e_pow);
        std::cout << "asp_power = " << a_pow << "\n";
        std::cout << "\n\n";
    }


}