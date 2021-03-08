#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;

//' Jensen–Shannon distance metric (i.e. the square root of the Jensen–Shannon divergence) between the columns of a dense matrix m
//'
//' @param m Input matrix
//' @return Vectorized version of the lower triangle as an R distance object, stats::dist()
//' @examples
//' ex = matrix(1:9, nrow = 3, ncol = 3)
//' jsDist(ex)
//'
// [[Rcpp::export]]
arma::mat jsDist(const arma::mat& m) {
  //arma::vec d(m.n_cols*(m.n_cols-1)/2);
  arma::mat d(m.n_cols, m.n_cols, arma::fill::zeros);
  for(int i=0; i<(m.n_cols-1); i++) {
    arma::vec li = log(m.col(i));
    for(int j=i+1; j<m.n_cols; j++) {
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
