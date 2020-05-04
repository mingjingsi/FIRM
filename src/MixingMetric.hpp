#ifndef MixingMetric_hpp
#define MixingMetric_hpp

#include <RcppArmadillo.h>
#include "nn.hpp"

using namespace arma;

vec MixingMetric(mat embedding, uvec dataset_list, const int k, const int max_k);

#endif
