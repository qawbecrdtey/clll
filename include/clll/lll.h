/**
 * @file
 * @brief LLL lattice basis reduction.
 */

#ifndef CLLL_LLL_H
#define CLLL_LLL_H

#include "vector.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Reduces a lattice basis in place with the LLL algorithm.
 *
 * This is Algorithm 5 (the LLL algorithm, verified integer version) from
 * Thiemann et al., "Formalizing the LLL Basis Reduction Algorithm and the LLL
 * Factorization Algorithm in Isabelle/HOL". It uses exact integer arithmetic
 * only.
 *
 * On return, basis[0], ..., basis[m-1] is a basis of the same lattice, and it
 * is reduced with respect to @f$\alpha@f$. That is, its Gram-Schmidt vectors
 * @f$g_i@f$ and coefficients @f$\mu_{i,j}@f$ satisfy
 * @f[
 *     |\mu_{i,j}| \le \frac{1}{2} \quad (j < i)
 *     \qquad \text{and} \qquad
 *     \|g_{i-1}\|^2 \le \alpha \|g_i\|^2 \quad (0 < i < m).
 * @f]
 * A lattice usually has many such bases. Which one is returned is not
 * specified beyond these conditions.
 *
 * Vectors are reordered by swapping whole structs, as clll_vector_swap() does.
 * Afterwards basis[i] is the i-th reduced vector, but a pointer such as
 * basis[i].v saved before the call may now point to a different vector.
 *
 * @param[in,out] basis  The vectors basis[0], ..., basis[m-1].
 * @param[in]     m      Number of vectors.
 * @param[in]     alpha  The approximation factor @f$\alpha@f$. The paper's
 *                       examples use 3/2.
 *
 * @pre basis[0], ..., basis[m-1] are linearly independent and all have the
 *      same length n, so @f$m \le n@f$.
 * @pre @f$\alpha \ge 4/3@f$. The running time is polynomial only for
 *      @f$\alpha > 4/3@f$.
 * @pre @f$\alpha@f$ is in canonical form, as GMP requires of every mpq_t
 *      operand: its denominator is positive and shares no factor with its
 *      numerator. Results of GMP's rational arithmetic are canonical. After
 *      mpq_set_ui() or mpq_set_si(), call mpq_canonicalize() unless the two
 *      arguments are already coprime.
 *
 * @warning The preconditions are not checked. If they do not hold, the
 *          behavior is undefined: the call may crash, never return, or corrupt
 *          memory. A build with assertions enabled catches some violations,
 *          but not all.
 */
void clll_lll(clll_vector *basis, unsigned long m, mpq_t const alpha);

#ifdef __cplusplus
}
#endif

#endif
