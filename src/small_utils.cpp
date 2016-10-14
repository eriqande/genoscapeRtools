#include <Rcpp.h>
using namespace Rcpp;

//' return a matrix of cumumlative sums along the rows or colums
//' @param x a numeric matrix
//' @param dim 1 for rows, 2 for columns
//' @export
// [[Rcpp::export]]
NumericMatrix mat_cumul_cpp(NumericMatrix x, int dim) {
  int i,j;
  int R = x.nrow();
  int C = x.ncol();
  double sum;
  NumericMatrix ret(R, C);


  if(dim == 1) {
    for(i=0;i<R;i++) {
      sum = 0.0;
      for(j=0;j<C;j++) {
        sum += x(i,j);
        ret(i,j) =  sum;
      }
    }
  } else if(dim == 2) {
    for(j=0;j<C;j++) {
      sum = 0.0;
      for(i=0;i<R;i++) {
        sum += x(i,j);
        ret(i,j) = sum;
      }
    }
  } else {
    stop("Dim must be 1 or 2");
  }

  return(ret);
}


