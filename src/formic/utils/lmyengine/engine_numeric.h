///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file cqmc/numeric/numeric.h
///
/// \brief   header file for miscellaneous functions related to numbers
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef CQMC_NUMERIC_HEADER
#define CQMC_NUMERIC_HEADER

#include<string>
#include<cassert>
#include<cmath>
#include<vector>
#include<set>

namespace cqmc {

  template<class S> void unbiased_ratio_of_means(const int n, const S * const p, const S * const f, const S * const g, const bool correct, S & r, S & v);
  template<class S> void mpi_unbiased_ratio_of_means(const int n, const S * const p, const S * const f, const S * const g, const bool correct, S & r, S & v);
  int my_round(const double d);

}

#endif
