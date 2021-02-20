// [[Rcpp::depends(RcppArmadillo)]]
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


// calculates factor-stratified sums for each column
// rowSel is an integer factor; 
// note that the 0-th column will return sums for any NA values; 0 or negative values will be omitted
// note2: trailing levels that are not mentioned in the vector will be omitted, resulting in a smaller matrix

//' Calculates factor-stratified sums for each column
//' 
//' @param sY sparse matrix (dgCmatrix)
//' @param rowSel integer factor. Note that the 0-th column will return sums for any NA values; 0 or negative values will be omitted
//' @return Matrix
// [[Rcpp::export]]
arma::mat colSumByFac(SEXP sY,  SEXP rowSel) {
// need to do this as SEXP, modify the slots on the fly
  S4 mat(sY);  
  const arma::uvec i(( unsigned int *)INTEGER(mat.slot("i")),LENGTH(mat.slot("i")),false,true); 
  const arma::ivec dims(INTEGER(mat.slot("Dim")),LENGTH(mat.slot("Dim")),false,true); 
  const arma::ivec p(INTEGER(mat.slot("p")),LENGTH(mat.slot("p")),false,true); 
  arma::vec Y(REAL(mat.slot("x")),LENGTH(mat.slot("x")),false,true); 
  
  const arma::ivec rs=arma::ivec(INTEGER(rowSel),LENGTH(rowSel),false,true);

  int ncols=p.size()-1;
  int nlevels=0;
  for(int j=0;j<rs.size();j++) { 
    if(rs[j]!=NA_INTEGER) {
      if(rs[j]>nlevels) { nlevels=rs[j]; }
    }
  }
  if(nlevels==0) { stop("colSumByFac(): supplied factor doesn't have any levels!"); }
  arma::mat sumM(nlevels+1,ncols,arma::fill::zeros);

  // for each gene
//#pragma omp parallel for shared(sumM) 
  for(int g=0;g<ncols;g++) {
    int p0=p[g]; int p1=p[g+1]; 
    if(p1-p0 <1) { continue; }
    for(int j=p0;j<p1;j++) {
      int row=i[j];
      int f=rs[row];
      if(f==NA_INTEGER) {
	sumM(0,g)+=Y[j];
      } else if(f>0) {
      	sumM(f,g)+=Y[j];
      }
    }
  }
  return sumM;
}

