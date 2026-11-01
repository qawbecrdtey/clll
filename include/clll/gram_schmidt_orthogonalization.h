/**
 * @file
 * @brief Integer Gram-Schmidt orthogonalization.
 */

#ifndef CLLL_GRAM_SCHMIDT_ORTHOGONALIZATION_H
#define CLLL_GRAM_SCHMIDT_ORTHOGONALIZATION_H

#include "vector.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes the scaled Gram-Schmidt coefficients of a basis.
 *
 * This is Algorithm 2 (GSO computation for the @f$\tilde\mu@f$-values only)
 * from Thiemann et al., "Formalizing the LLL Basis Reduction Algorithm and the
 * LLL Factorization Algorithm in Isabelle/HOL".
 *
 * For vectors @f$f_0, \ldots, f_{m-1} \in \mathbb{Z}^n@f$ it computes
 * @f[
 *     \tilde\mu_{i,j} = d_{j+1}\,\mu_{i,j} \qquad (j \le i),
 * @f]
 * where @f$\mu_{i,j}@f$ are the Gram-Schmidt coefficients, @f$d_0 = 1@f$, and
 * @f$d_k@f$ for @f$k \ge 1@f$ is the determinant of the Gram matrix of
 * @f$f_0, \ldots, f_{k-1}@f$. These values are integers, and
 * @f$\tilde\mu_{i,i} = d_{i+1}@f$.
 *
 * Row i of the result goes to mu[i]: entries 0 through i are written and the
 * rest are left unchanged.
 *
 * @param[out] mu  The rows mu[0], ..., mu[m-1].
 * @param[in]  fs  The vectors fs[0], ..., fs[m-1].
 * @param[in]  m   Number of vectors.
 *
 * @pre fs[0], ..., fs[m-1] are linearly independent and all have the same
 *      length n, so @f$m \le n@f$.
 * @pre mu[0], ..., mu[m-1] are initialized with length n.
 * @pre No mu[i] shares entries with any fs[k]. In particular, mu and fs must
 *      not be the same array.
 *
 * @warning The preconditions are not checked. If they do not hold, the
 *          behavior is undefined.
 */
void clll_gram_schmidt_orthogonalization(clll_vector *mu, clll_vector const *fs, unsigned long m);

#ifdef __cplusplus
}
#endif

#endif
