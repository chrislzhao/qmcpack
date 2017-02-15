///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file cqmc/numeric/numeric.cpp
///
/// \brief   implementation file for miscellaneous functions related to numbers
///
///////////////////////////////////////////////////////////////////////////////////////////////////

//#include<mpi.h>
#include<formic/utils/mpi_interface.h>

#include<formic/utils/lmyengine/engine_numeric.h>


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   computes the unbiased estimate of a ratio of means:  <f>_p / <g>_p  in which the
///          numerator and denominator values are sampled from the same probability distribution p
///
/// \param[in]     n        the number of samples
/// \param[in]     p        the probability weight for each sample
/// \param[in]     f        the numerator samples
/// \param[in]     g        the denominator samples
/// \param[out]    r        on exit, the estimate of the ratio <f>_p / <g>_p
/// \param[out]    v        on exit, the estimate of the variance in the ratio
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void cqmc::unbiased_ratio_of_means(const int n, const S * const p, const S * const f, const S * const g, const bool correct, S & r, S & v) {

  // compute the normalization, the numerator and denominator means, the means of the squares, and the mean of the products
  S nm = 0.0; // normalization constant
  S mf = 0.0; // mean of numerator
  S mg = 0.0; // mean of denominator
  S sf = 0.0; // mean of the square of the numerator terms
  S sfnc = 0.0; // mean of the square of the numerator used in exact sampling 
  S sg = 0.0; // mean of the square of the denominator terms
  S mp = 0.0; // mean of the product of numerator times denominator
  for (int i = 0; i < n; i++) {
    nm += p[i];
    S x = p[i] * f[i];
    mf += x;
    if ( correct ) {
      sf += x * f[i];
      mp += x * g[i];
    }
    else 
      sfnc += p[i] * (f[i] / g[i]) * (f[i] / g[i]) * g[i];
    x = p[i] * g[i];
    mg += x;
    if ( correct )
      sg += x * g[i];
  }
  mf /= nm;
  mg /= nm;
  sf /= nm;
  sfnc /= nm;
  sg /= nm;
  mp /= nm;

  // compute the numerator and denominator variances and the covariance
  const S vf = ( sf - mf * mf ) * double(n) / double(n-1);
  const S vg = ( sg - mg * mg ) * double(n) / double(n-1);
  const S cv = ( mp - mf * mg ) * double(n) / double(n-1);

  double z = ( correct ? unity(S()) : zero(S()));

  // compute the unbiased estimate of the ratio of means
  r = ( mf / mg ) / ( unity(S()) + z * ( vg / mg / mg - cv / mf / mg ) / double(n) );

  // compute the unbiased estimate of the variance of the ratio of means
  if ( correct ) 
    v = ( mf * mf / mg / mg / double(n) ) * ( vf / mf / mf + vg / mg / mg - 2.0 * cv / mf / mg );
  else 
    v = sfnc - mf * mf; 

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   computes the unbiased estimate of a ratio of means:  <f>_p / <g>_p  in which the
///          numerator and denominator values are sampled from the same probability distribution p
///          and samples are combined across all processors
///
/// \param[in]     n        the number of samples on this process
/// \param[in]     p        the probability weight for each sample
/// \param[in]     f        the numerator samples
/// \param[in]     g        the denominator samples
/// \param[out]    r        on exit, the estimate of the ratio <f>_p / <g>_p
/// \param[out]    v        on exit, the estimate of the variance in the ratio
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void cqmc::mpi_unbiased_ratio_of_means(const int n, const S * const p, const S * const f, const S * const g, const bool correct, S & r, S & v) {

  // compute the normalization, the numerator and denominator means, the means of the squares, and the mean of the products
  S y[8];
  y[0] = zero(S()); // normalization constant
  y[1] = zero(S()); // mean of numerator
  y[2] = zero(S()); // mean of denominator
  y[3] = zero(S()); // mean of the square of the numerator terms
  y[4] = zero(S()); // mean of the square of the denominator terms
  y[5] = zero(S()); // mean of the product of numerator times denominator
  y[6] = double(n) * unity(S()); // number of samples
  y[7] = zero(S()); // mean of square of numerator terms used in exact sampling(when correct is false) 
  for (int i = 0; i < n; i++) {
    y[0] += p[i];
    S x = p[i] * f[i];
    y[1] += x;
    if ( correct ) {
      y[3] += x * f[i];
      y[5] += x * g[i];
    }
    else 
      y[7] += p[i] * f[i] * f[i];
    x = p[i] * g[i];
    y[2] += x;
    if ( correct )
      y[4] += x * g[i];
  }
  S z[8];
  formic::mpi::allreduce(&y[0], &z[0], 8, MPI::SUM);
  const S mf = z[1] / z[0]; // mean of numerator
  const S mg = z[2] / z[0]; // mean of denominator
  const S sf = z[3] / z[0]; // mean of the square of the numerator terms
  const S sg = z[4] / z[0]; // mean of the square of the denominator terms
  const S mp = z[5] / z[0]; // mean of the product of numerator times denominator
  const S sfnc = z[7] / z[0]; // mean of square of the numerator used in exact sampling  
  const S ns = z[6];        // number of samples

  // compute the numerator and denominator variances and the covariance
  const S vf = ( sf - mf * mf ) * ns / ( ns - 1.0 );
  const S vg = ( sg - mg * mg ) * ns / ( ns - 1.0 );
  const S cv = ( mp - mf * mg ) * ns / ( ns - 1.0 );

  S c = (correct ? unity(S()) : zero(S()));

  // compute the unbiased estimate of the ratio of means
  r = ( mf / mg ) / ( unity(S()) + c * ( vg / mg / mg - cv / mf / mg ) / ns );

  // compute the unbiased estimate of the variance of the ratio of means
  if ( correct ) 
    v = ( mf * mf / mg / mg ) * ( vf / mf / mf + vg / mg / mg - 2.0 * cv / mf / mg );
  else 
    v = sfnc - mf * mf;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   rounds a double to the nearest integer
///
/// \param[in]    d          the number to round
///
/// \return the integer nearst to d
///////////////////////////////////////////////////////////////////////////////////////////////////
int cqmc::my_round(const double d) { return ( d >= 0 ? int(d + 0.5) : int(d - 0.5) ); }
