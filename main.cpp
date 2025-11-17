#include "boost/random/cauchy_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
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

double g2(double x){
    return fabs(x);
}

double d2_g2(double x){
    return 0;
}



// cauchy with h2=0
void ex3(){
    std::vector<double>  h1_vals{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double h2 = 0.0;
    double alpha = 0.05;

    std::vector<int> sample_sizes{100, 400, 900, 1600, 2000, 3000};
    // std::vector<int> sample_sizes{100, 300, 400, 500, 600, 750, 900, 1100, 1300, 1450, 1600, 1800, 2000, 2500, 3000};
    int N = 5000;
    int M = 5000;
    

    std::vector<double> integrals = read_integrals(DistributionType::CAUCHY);

    std::cout<< "ET_CAUCHY:\n";
    run_experiment(g, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_etest, true);
    std::cout << "================================\n";
    std::cout << "================================\n\n";

    // std::cout<< "HT_CAUCHY:\n";
    // run_experiment(g2, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_etest,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "WMW_CAUCHY:\n";
    // run_experiment(g, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_wmw, false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "AD_CAUCHY:\n";
    // run_experiment(g, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_ad,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "KS_CAUCHY:\n";
    // run_experiment(g, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_ks,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "CM_CAUCHY:\n";
    // run_experiment(g, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_cm,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";


}

// cauchy with h1=0
void ex4(){
    std::vector<double>  h2_vals{1, 2, 3, 4, 5, 6, 7, 8, 9};
    double h1 = 0.0;
    double alpha = 0.05;

    std::vector<int> sample_sizes{100};
    // std::vector<int> sample_sizes{100, 300, 400, 500, 600, 750, 900, 1100, 1300, 1450, 1600, 1800, 2000, 2500, 3000};
    int N = 5000;
    int M = 5000;
    

    std::vector<double> integrals = read_integrals(DistributionType::CAUCHY);

    // std::cout<< "ET_CAUCHY:\n";
    // run_experiment(g, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_etest, true);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "HT_CAUCHY:\n";
    // run_experiment(g2, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_etest,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "WMW_CAUCHY:\n";
    // run_experiment(g, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_wmw,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "AD_CAUCHY:\n";
    // run_experiment(g, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_ad,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "KS_CAUCHY:\n";
    // run_experiment(g, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_ks, false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "CM_CAUCHY:\n";
    // run_experiment(g, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_cauchy, integrals, compute_cm,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";
}

// normal with h1=0
void ex1(){
    std::vector<double>  h2_vals{1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5};
    double h1 = 0.0;
    double alpha = 0.05;
// 
    std::vector<int> sample_sizes{100};
    // std::vector<int> sample_sizes{100, 300, 400, 500, 600, 750, 900, 1100, 1300, 1450, 1600, 1800, 2000, 2500, 3000};
    int N = 5000;
    int M = 5000;
    


    std::vector<double> integrals = read_integrals(DistributionType::NORMAL);

    // std::cout<< "ET_NORMAL:\n";
    // run_experiment(g, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_normal, integrals, compute_etest, true);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    std::cout<< "HT_NORMAL:\n";
    run_experiment(g2, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_normal, integrals, compute_etest,  true);
    std::cout << "================================\n";
    std::cout << "================================\n\n";

    // std::cout<< "WMW_NORMAL:\n";
    // run_experiment(g, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_normal, integrals, compute_wmw,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "AD_NORMAL:\n";
    // run_experiment(g, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_normal, integrals, compute_ad, false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "KS_NORMAL:\n";
    // run_experiment(g, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_normal, integrals, compute_ks,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "CM_NORMAL:\n";
    // run_experiment(g, d2_g, h1, h2_vals, alpha, N, M, sample_sizes, get_normal, integrals, compute_cm,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";
}




// normal with h2=0
void ex2(){
    std::vector<double>  h1_vals{1,3, 5};
    double h2 = 0.0;
    double alpha = 0.05;

    std::vector<int> sample_sizes{100, 200};
    // std::vector<int> sample_sizes{100, 300, 400, 500, 600, 750, 900, 1100, 1300, 1450, 1600, 1800, 2000, 2500, 3000};
    int N = 5000;
    int M = 5000;
    

    std::vector<double> integrals = read_integrals(DistributionType::NORMAL);

    std::cout<< "ET_NORMAL:\n";
    run_experiment(g, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_normal, integrals, compute_etest,  true);
    std::cout << "================================\n";
    std::cout << "================================\n\n";

    // std::cout<< "HT_NORMAL:\n";
    // run_experiment(g2, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_normal, integrals, compute_etest,  true);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";



    // std::cout<< "WMW_NORMAL:\n";
    // run_experiment(g, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_normal, integrals, compute_wmw, false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "AD_NORMAL:\n";
    // run_experiment(g, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_normal, integrals, compute_ad, false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "KS_NORMAL:\n";
    // run_experiment(g, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_normal, integrals, compute_ks,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "CM_NORMAL:\n";
    // run_experiment(g, d2_g, h1_vals, h2, alpha, N, M, sample_sizes, get_normal, integrals, compute_cm,  false);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";
}

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

}


int main() {
    // ex1();
    ex2();
    // ex3();
    // ex4();
}