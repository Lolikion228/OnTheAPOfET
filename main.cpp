#include "boost/math/distributions/cauchy.hpp"
#include "boost/random/cauchy_distribution.hpp"
#include <boost/random.hpp>
#include <iostream>



int main() {
    boost::random::mt19937 gen(std::time(0));

    boost::random::cauchy_distribution<> mydist(0, 1);

    std::cout << mydist(gen) << "\n";
}