// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <queue>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <vector>
#include <queue>
#include <list>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace Rcpp;

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

//' Jensen–Shannon distance metric (i.e. the square root of the Jensen–Shannon divergence) between the columns of a dense matrix m
//'
//' @param m Input matrix
//' @param ncores integer Number of cores (default=1)
//' @return Vectorized version of the lower triangle as an R distance object, stats::dist()
//' @examples
//' ex = matrix(1:9, nrow = 3, ncol = 3)
//' jsDist(ex)
//'
// [[Rcpp::export]]
arma::mat jsDist(const arma::mat& m, int ncores=1) {
  //arma::vec d(m.n_cols*(m.n_cols-1)/2);
  arma::mat d(m.n_cols,m.n_cols,arma::fill::zeros);
//#pragma omp parallel for num_threads(ncores) shared(d)
  for(int i=0;i<(m.n_cols-1);i++) {
    arma::vec li=log(m.col(i));
    for(int j=i+1;j<m.n_cols;j++) {
      arma::vec lj=log(m.col(j));
      arma::vec ji=m.col(j)+m.col(i);
      ji=m.col(j)%lj + m.col(i)%li - ji%(log(ji)-log(2.0)); 
      double v=arma::accu(ji.elem(arma::find_finite(ji))); 
      //d[m.n_cols*i - i*(i-1)/2 + j-i]=sqrt(v/2.0);
      d(j,i)=d(i,j)=v;
    }
  }

  return(d);
}
