// dissimilarity.cpp
// -----------------------------------------------------------------------------
// Pairwise posterior dissimilarity (1 - posterior similarity matrix) from a
// matrix of MCMC cluster-label draws. Moved out of an inline cppFunction so it
// is compiled once at package install time and exported into the namespace.
// -----------------------------------------------------------------------------

// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

//' Posterior dissimilarity matrix from MCMC label draws
//'
//' For every pair of genes, returns one minus the fraction of MCMC draws in
//' which the two genes share a cluster label (i.e. one minus their posterior
//' co-assignment probability). The diagonal is zero.
//'
//' @param zmat An integer matrix of label draws with one row per retained draw
//'   and one column per gene.
//' @return A symmetric numeric matrix of dissimilarities (genes by genes).
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix calcDissimilarityMatrixCpp(IntegerMatrix zmat)
{
    int n = zmat.ncol();
    int s = zmat.nrow();
    NumericMatrix dissimMatrix(n, n);
    std::fill(dissimMatrix.begin(), dissimMatrix.end(), 0.0);

    for (int i = 0; i < n; i++)
    {
        dissimMatrix(i, i) = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            int same = 0;
            for (int k = 0; k < s; k++)
            {
                if (zmat(k, i) == zmat(k, j))
                    same++;
            }
            double dissim = 1.0 - (double)same / s;
            dissimMatrix(i, j) = dissim;
            dissimMatrix(j, i) = dissim;
        }
    }
    return dissimMatrix;
}