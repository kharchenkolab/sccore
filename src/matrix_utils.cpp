#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;

//' Jensen–Shannon distance metric (i.e. the square root of the Jensen–Shannon divergence) between the columns of a dense matrix m
//'
//' @param m Input matrix
//' @param ncores Number of threads to be set via omp_set_num_threads() for RcppArmadillo 
//' @return Vectorized version of the lower triangle as an R distance object, stats::dist()
//' @examples
//' ex = matrix(1:9, nrow = 3, ncol = 3)
//' # JS distance calculated between columns of input matrix
//' jsDist(ex)
//'
//'
//' # To demonstrate how the above JS Distance to the JS Divergence, 
//' # we use the third-party function 'philentropy::JSD()', which 
//' # computes the JS divergence between rows of the input matrix. 
//' # The following will give the same results as 'jsDist(ex)':
//' sqrt(philentropy::JSD(t(ex), est.prob = "empirical"))
//'
//'
//' # Conversely, we can use the column-normalized matrix, 
//' # and ignore the argument 'est.prob = "empirical"' from 'philentropy::JSD()', 
//' # which calculates the relative frequencies of each vector are computed internally). 
//' # This again will give the same results as 'jsDist(ex)':
//' ex_cnorm = t( t(ex)/colSums(ex) )
//' sqrt(philentropy::JSD(t(ex_cnorm)))
//'
//'
//' # Again obviously 'jsDist(ex)**2' will be the JS divergence, 
//' # equaling 'philentropy::JSD(t(ex_cnorm))' and 'philentropy::JSD(t(ex), est.prob = "empirical")'
//' jsDist(ex)**2 
//' philentropy::JSD(t(ex_cnorm))
//' philentropy::JSD(t(ex), est.prob = "empirical")
//'
//'
// [[Rcpp::export]]
arma::mat jsDist(const arma::mat& m, int ncores=1) {


  // internals of void armadillo_set_number_of_omp_threads(int n), 
  // cf https://github.com/RcppCore/RcppArmadillo/blob/43708298d73f0e23328c6e5f20575c980147e3c5/src/RcppArmadillo.cpp#L107-L115
  // 
  #ifdef _OPENMP
    omp_set_num_threads(ncores); 
  #else
    (void)(ncores);                  // prevent unused variable warning
  #endif

  arma::mat x = m;                  
  arma::rowvec s = arma::sum(x, 0); 
  if (arma::any(s <= 0)){
    Rcpp::stop("All columns must have positive sum.");
  }                            
  x.each_row() /= s;

  arma::mat d(x.n_cols, x.n_cols, arma::fill::zeros);
  const double log2 = std::log(2.0);

  for(int i=0; i<(x.n_cols-1); i++) {
    arma::vec li = log(x.col(i));
    li.elem(arma::find_nonfinite(li)).zeros();

    for(int j=i+1; j<x.n_cols; j++) {
      arma::vec lj = log(x.col(j));
      lj.elem(arma::find_nonfinite(lj)).zeros();

      arma::vec ji = x.col(j) + x.col(i);
      arma::vec lji = log(ji);
      lji.elem(arma::find_nonfinite(lji)).zeros();

      ji = x.col(j)%lj + x.col(i)%li - ji%(lji - log2);
      double v = arma::accu(ji);

      d(j,i) = d(i,j) = std::sqrt(std::max(v, 0.0) / (2.0 * log2));
    }
  }

  return(d); 
}

//' Calculates factor-stratified sums for each column
//'
//' @param sY sparse matrix (dgCmatrix)
//' @param rowSel integer factor. Note that the 0-th column will return sums for any NA values; 0 or negative values will be omitted
//' @return Matrix
// [[Rcpp::export]]
NumericMatrix colSumByFactor(SEXP sY,  IntegerVector rowSel) {
  // note: trailing levels that are not mentioned in the vector will be omitted, resulting in a smaller matrix
  // need to do this as SEXP, modify the slots on the fly
  S4 mat(sY);
  const arma::uvec i(( unsigned int *)INTEGER(mat.slot("i")),LENGTH(mat.slot("i")),false,true);
  const arma::ivec dims(INTEGER(mat.slot("Dim")),LENGTH(mat.slot("Dim")),false,true);
  const arma::ivec p(INTEGER(mat.slot("p")),LENGTH(mat.slot("p")),false,true);
  arma::vec Y(REAL(mat.slot("x")),LENGTH(mat.slot("x")),false,true);

  List dimnames(mat.slot("Dimnames"));
  CharacterVector geneNames(dimnames[1]);


  CharacterVector factorLevels = rowSel.attr("levels");
  int nlevels=factorLevels.size();
  CharacterVector expandedFactorLevels(nlevels+1);
  expandedFactorLevels[0]="NA";
  for(int i=0;i<nlevels;i++) { expandedFactorLevels[i+1]=factorLevels[i]; }
  const arma::ivec rs=arma::ivec(INTEGER(rowSel),LENGTH(rowSel),false,true);

  int ncols=p.size()-1;
  if(nlevels==0) { stop("colSumByFactor(): supplied factor doesn't have any levels!"); }
  NumericMatrix sumM(nlevels+1, ncols);

  // for each gene
  for(int g=0;g<ncols;g++) {
    int p0=p[g]; int p1=p[g+1];
    if(p1-p0 <1) { continue; }
    for(int j=p0;j<p1;j++) {
      int row=i[j];
      int f=rs[row];
      if(f==NA_INTEGER) {
        sumM(0,g) += Y[j];
      } else if(f>0) {
        sumM(f,g) += Y[j];
      }
    }
  }

  colnames(sumM) = geneNames;
  rownames(sumM) = expandedFactorLevels;
  return sumM;
}
