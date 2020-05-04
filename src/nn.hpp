#ifndef nn_hpp
#define nn_hpp

#include <RcppArmadillo.h>
#include <ANN/ANN.h>

using namespace arma;

umat nn2_cpp(mat data, const int k);

#endif
