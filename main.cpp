
#include "boost/random/cauchy_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <iostream>
#include "first.h"



int main() {
    const int random_seed = 112;
    boost::random::mt19937 gen(random_seed);

    boost::random::normal_distribution<> dist1(0, 1);
    double *sample1 = sample(gen, dist1, 10);
 

    print_sample(sample1, 10);

    
}