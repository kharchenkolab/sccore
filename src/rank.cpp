#include <RcppArmadillo.h>
#include <algorithm>

using namespace std;
using namespace Rcpp;

// Average-rank helpers for the Wilcoxon/Mann-Whitney differential-expression core (matrixDE()).
// Vendored from pagoda2 (rank.cpp) so the same rank statistic can be shared by pagoda2 and conos.
// thanks to nrussel of stackoverflow!
class RankComparator {
private:
    const Rcpp::NumericVector& ref;
    bool is_na(double x) const { return Rcpp::traits::is_na<REALSXP>(x); }
public:
    RankComparator(const Rcpp::NumericVector& ref_) : ref(ref_) {}
    bool operator()(const int ilhs, const int irhs) const {
        double lhs = ref[ilhs], rhs = ref[irhs];
        if (is_na(lhs)) return false;
        if (is_na(rhs)) return true;
        return lhs < rhs;
    }
};

//' Average ranks of a numeric vector (ties resolved by mid-rank)
//'
//' @param x numeric vector
//' @return numeric vector of average ranks (1-based), NA values ranked last
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector avg_rank(Rcpp::NumericVector x) {
    R_xlen_t sz = x.size();
    Rcpp::IntegerVector w = Rcpp::seq(0, sz - 1);
    std::sort(w.begin(), w.end(), RankComparator(x));

    Rcpp::NumericVector r = Rcpp::no_init_vector(sz);
    for (R_xlen_t n, i = 0; i < sz; i += n) {
        n = 1;
        while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
        for (R_xlen_t k = 0; k < n; k++) {
            r[w[i + k]] = i + (n + 1) / 2.;
        }
    }
    return r;
}

class BlockRankComparator {
private:
  const Rcpp::NumericVector::const_iterator it;
  bool is_na(double x) const { return Rcpp::traits::is_na<REALSXP>(x); }
public:
  BlockRankComparator(const Rcpp::NumericVector::const_iterator it_) : it(it_) {}
  bool operator()(const int ilhs, const int irhs) const {
    double lhs = *(it+ilhs), rhs = *(it+irhs);
    if (is_na(lhs)) return false;
    if (is_na(rhs)) return true;
    return lhs < rhs;
  }
};

//' Per-column average ranks of a sparse (dgCMatrix) matrix, with implicit zeros ranked correctly
//'
//' Ranks the entries within each column (gene), shifting the non-zero ranks up by the number of zeros so
//' that the (unstored) zeros occupy the lowest ranks. Used by matrixDE() for the Mann-Whitney statistic.
//'
//' @param sY a sparse column-compressed (dgCMatrix) matrix
//' @return a dgCMatrix of the same sparsity pattern whose stored values are the within-column ranks
//' @keywords internal
// [[Rcpp::export]]
S4 sparse_matrix_column_ranks(const SEXP sY) {
  const S4 mat(sY);
  const NumericVector x = mat.slot("x");
  const IntegerVector p = mat.slot("p");
  const IntegerVector dims = mat.slot("Dim");

  int ncols = p.size() - 1;
  int nrows = dims[0];

  NumericVector xbr(x.size());

  // iterate over columns
  for (int g = 0; g < ncols; g++) {
    int p0 = p[g]; int p1 = p[g+1];
    int sz = p1 - p0;
    if (sz < 1) { continue; }
    Rcpp::IntegerVector w = Rcpp::seq(0, sz - 1);
    std::sort(w.begin(), w.end(), BlockRankComparator(x.begin() + p0));

    // determine the number of zeroes
    int nzeros = nrows - sz;
    for (R_xlen_t n, i = 0; i < sz; i += n) {
      n = 1;
      while (i + n < sz && x[p0 + w[i]] == x[p0 + w[i + n]]) ++n;
      for (R_xlen_t k = 0; k < n; k++) {
        xbr[p0 + w[i + k]] = nzeros + i + (n + 1) / 2.;
      }
    }
  }

  S4 r("dgCMatrix");
  r.slot("i") = mat.slot("i");
  r.slot("p") = mat.slot("p");
  r.slot("Dim") = mat.slot("Dim");
  r.slot("Dimnames") = mat.slot("Dimnames");
  r.slot("x") = xbr;
  r.slot("factors") = mat.slot("factors");
  return r;
}
