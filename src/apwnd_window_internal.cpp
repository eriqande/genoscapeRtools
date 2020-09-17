#include <Rcpp.h>
using namespace Rcpp;

//' Compute average pairwise nucleotide diversity statistics for positions from a to b
//'
//' This is the internal Rcpp function that will be called once for each window.
//' @param x an 0,1,2,-1 integer matrix
//' @param a the lowest (base-1) index of the SNPs in the window
//' @param b the highest (base-1) index of the SNPs in the window
//' @return This returns a list of vectors.  First vector is the average pairwise nucleotide difference
//' for each pair of individuals.  Second vector is the fraction of sites missing in the pair.
//' @export
// [[Rcpp::export]]
List apwnd_window_internal(IntegerMatrix x, int a, int b) {

  int i, j, l;
  double sum;  // for summing up the average nucleotide diffs at each SNP
  double msum;  // for summing up the number of sites missing data at at least on member of the pair
  int N = x.nrow();
  NumericVector nd(N * (N - 1) / 2);  // to return all the average pairwise nucleotide diversity values
  NumericVector fm(N * (N - 1) / 2);  // to return fraction missing for each pair
  int y1, y2;  // to hold the genotypes
  double tmp;  // to hold a value before making it positive if need be
  int pair = -1; // to keep track of the index of the pair

  // cycle over all pairs of individuals
  for(i = 0; i < (N - 1); i++) {
    for(j = i + 1; j < N; j++) {
      pair++; // increment the pair index
      sum = 0.0; // initialize to zero to sum up over all SNPs
      msum = 0.0;

      for(l = (a - 1); l <= (b - 1); l++) {  // cycle over positions
        y1 = x(i, l);
        y2 = x(j, l);

        if(y1 == -1 || y2 == -1) {
          msum += 1.0;
        }
        else {
          if(y1 == 1 && y2 == 1) {
            sum += 0.5;
          }
          else {
            tmp = 0.5 * (y1 - y2);
            if(tmp < 0.0) tmp = -1.0 * tmp;  // make it positive
            sum += tmp;
          }
        }
      } // close loop over l

      // transfer values to the output vectors
      nd[pair] = sum;
      fm[pair] = msum;

    }  // close loop over j
  } // close loop over i

  // return that as a named list
  List ret;
  ret["nd"] = nd;
  ret["msum"] = fm;
  return(ret);

}
