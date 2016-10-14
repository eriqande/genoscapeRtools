# small utilities in R. mostly not exported

mat_cumul_R <- function(x, dim) {
  ret <- mat_cumul_cpp(x, dim)
  dimnames(ret) <- dimnames(x)
  ret
}
