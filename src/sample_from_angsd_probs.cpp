#include <Rcpp.h>
using namespace Rcpp;

//' Sample a genotype (0,1,2) from probabilities for all the snps and individuals
//' in x.
//' @param x a numeric vector that holds the posterior probs from an
//' ANGSD geno.gz file (from doGenos 32)
//' @param samples character vector of sample names
//' @param snps character vector of SNP names
//' @param thresh any snp in an individual must have one posterior greater than T in order to
//' have a non-missing value simulated for it.  For example, if T is 0.6 and the posteriors
//' are 0.33, 0.33, 0.33, this locus will just get a missing (-1).
//' @export
// [[Rcpp::export]]
IntegerMatrix sample_from_angsd_probs(NumericVector x, CharacterVector samples, CharacterVector snps, double thresh) {

  int i,j,k, ri=0, start=0, notMissing;
  double r, cumul;
  int N = samples.length();
  int L = snps.length();
  IntegerMatrix ret(N, L);  // matrix to return the values in
  NumericVector randos(N * L);  // place to store simulated uniform(0,1).

  randos = runif(N * L, 0.0, 1.0); // get all you uniform rvs

  for(j=0;j<L;j++) {  // cycle over the snps
    for(i=0;i<N;i++) {  // cycle over the individuals
      r = randos[ri++];  // get the random number that we are going to work with here
      cumul = 0.0;
      notMissing = 0;
      for(k=0;k<3;k++) {  // cycle over the three possible SNP genotypes first and confirm that one is over T
        if(x[start + k] > thresh) {
          notMissing = 1;
          break;
        }
      }
      if(notMissing) {
        ret(i,j) = -2;  // set this so that we can see if we ever fail to give value here...
        for(k=0;k<3;k++) {
          cumul += x[start + k];
          if(cumul>=r) {
            ret(i,j) = k;
            break;
          }
        }
      } else {
        ret(i,j) = -1;
      }
      start += 3;  // go to the next set of 3 probs for the next one
    }
  }

  rownames(ret) = samples;
  colnames(ret) = snps;
  return(ret);
}
