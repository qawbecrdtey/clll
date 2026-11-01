/**
 * @file
 * @brief Integer vectors with GMP entries.
 */

#ifndef CLLL_VECTOR_H
#define CLLL_VECTOR_H

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief A vector in @f$\mathbb{Z}^{len}@f$.
 *
 * The functions follow GMP's conventions. Initialize a vector before use and
 * clear it afterwards. The result comes first, and it may be the same vector
 * as an operand. All vectors passed to one call must have the same length;
 * otherwise the behavior is undefined.
 */
typedef struct {
    unsigned long len; /**< Number of entries. */
    mpz_t *v;          /**< The entries v[0], ..., v[len-1]. Read and write them directly. */
} clll_vector;

/**
 * @brief Initializes @p x to the zero vector of length @p len.
 *
 * @param[out] x    The vector to initialize.
 * @param[in]  len  Its length.
 */
void clll_vector_init(clll_vector *x, unsigned long len);

/**
 * @brief Frees the memory used by @p x.
 *
 * @param[in,out] x  An initialized vector. Initialize it again before reusing it.
 */
void clll_vector_clear(clll_vector *x);

/**
 * @brief Copies @p op into @p rop.
 *
 * @param[out] rop  The destination.
 * @param[in]  op   The source.
 */
void clll_vector_set(clll_vector *rop, clll_vector const *op);

/**
 * @brief Swaps @p x and @p y in constant time.
 *
 * Only the structs are exchanged, so each entry array moves along with its
 * vector.
 *
 * @param[in,out] x  The first vector.
 * @param[in,out] y  The second vector.
 */
void clll_vector_swap(clll_vector *x, clll_vector *y);

/**
 * @brief Sets @p rop to the inner product @f$u \cdot v@f$.
 *
 * @param[out] rop  The result. It must not be an entry of @p u or @p v.
 * @param[in]  u    The first vector.
 * @param[in]  v    The second vector.
 *
 * @warning The aliasing rule for @p rop is not checked. Breaking it is
 *          undefined behavior.
 */
void clll_vector_inner_product(mpz_t rop, clll_vector const *u, clll_vector const *v);

/**
 * @brief Sets @p rop to @f$u + v@f$.
 *
 * @param[out] rop  The result.
 * @param[in]  u    The first summand.
 * @param[in]  v    The second summand.
 */
void clll_vector_add(clll_vector *rop, clll_vector const *u, clll_vector const *v);

/**
 * @brief Sets @p rop to @f$u - v@f$.
 *
 * @param[out] rop  The result.
 * @param[in]  u    The minuend.
 * @param[in]  v    The subtrahend.
 */
void clll_vector_sub(clll_vector *rop, clll_vector const *u, clll_vector const *v);

/**
 * @brief Sets @p rop to the entrywise product of @p u and @p v.
 *
 * @param[out] rop  The result.
 * @param[in]  u    The first factor.
 * @param[in]  v    The second factor.
 */
void clll_vector_mul(clll_vector *rop, clll_vector const *u, clll_vector const *v);

#ifdef __cplusplus
}
#endif

#endif
