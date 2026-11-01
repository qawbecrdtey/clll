/**
 * @file
 * @brief Packed storage for @f$\tilde\mu@f$ and the Gramian determinants, and
 *        the row step of Algorithm 2.
 *
 * The lower triangle of @f$\tilde\mu@f$ and the determinants
 * @f$d_0, \ldots, d_m@f$ share one array of @f$m(m+1)/2 + 1@f$ entries. Row i
 * holds @f$d_i@f$ followed by @f$\tilde\mu_{i,0}, \ldots, \tilde\mu_{i,i-1}@f$,
 * and @f$d_m@f$ comes last:
 * @f[
 *     d_0 \;\big|\; d_1\ \tilde\mu_{1,0} \;\big|\;
 *     d_2\ \tilde\mu_{2,0}\ \tilde\mu_{2,1} \;\big|\; \cdots \;\big|\; d_m
 * @f]
 * So @f$d_k@f$ is entry @f$k(k+1)/2@f$ and @f$\tilde\mu_{i,j}@f$ is entry
 * @f$i(i+1)/2 + 1 + j@f$. The diagonal entry @f$\tilde\mu_{i,i}@f$ lands on
 * the same index as @f$d_{i+1}@f$, so the two are stored once.
 *
 * Indices are computed in size_t. Since @f$m \le n@f$, the array has at most
 * @f$mn + 1@f$ entries, which is no more than the basis itself.
 */

#ifndef CLLL_GRAM_SCHMIDT_ORTHOGONALIZATION_ROW_H
#define CLLL_GRAM_SCHMIDT_ORTHOGONALIZATION_ROW_H

#include <clll/vector.h>

#include "assert_and_assume.h"

#include <gmp.h>

#include <stddef.h>

/**
 * @brief Index of @f$d_k@f$ in the packed array.
 *
 * @param[in] k  The index of the determinant, @f$0 \le k \le m@f$.
 * @return @f$k(k+1)/2@f$.
 */
static inline CLLL_CONST_FUNCTION size_t d_index(unsigned long const k) {
    return (size_t)k * (k + 1) / 2;
}

/**
 * @brief Index of @f$\tilde\mu_{i,j}@f$ in the packed array.
 *
 * @param[in] i  The row.
 * @param[in] j  The column, @f$0 \le j \le i@f$.
 * @return @f$i(i+1)/2 + 1 + j@f$.
 */
static inline CLLL_CONST_FUNCTION size_t mu_index(unsigned long const i, unsigned long const j) {
    return d_index(i) + 1 + j;
}

/**
 * @brief Computes row i of @f$\tilde\mu@f$ with Algorithm 2.
 *
 * Writes the entries mu_index(i, 0), ..., mu_index(i, i) of @p mu. The last
 * of them is @f$d_{i+1}@f$.
 *
 * @param[in,out] mu     The packed array.
 * @param[in]     fs     The vectors. Only fs[0], ..., fs[i] are read.
 * @param[in]     i      The row to compute.
 * @param[in,out] sigma  Scratch space.
 * @param[in,out] t      Scratch space.
 *
 * @pre Rows 0, ..., i-1 of @p mu, and with them @f$d_1, \ldots, d_i@f$, are
 *      up to date.
 * @pre @p sigma and @p t are initialized and are not entries of @p mu or
 *      @p fs. Their values are ignored and are left unspecified.
 */
void clll_gram_schmidt_orthogonalization_row(mpz_t *mu, clll_vector const *fs, unsigned long i,
                                             mpz_t sigma, mpz_t t);

#endif
