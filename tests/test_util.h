/**
 * @file
 * @brief Allocation, a fixed random generator and arrays of vectors for the
 *        tests.
 */

#ifndef CLLL_TESTS_TEST_UTIL_H
#define CLLL_TESTS_TEST_UTIL_H

#include <clll.h>

#include <gmp.h>

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

/**
 * @brief Allocates @p size bytes of zeroed memory, or aborts.
 *
 * The memory is zeroed because otherwise GCC warns that the returned arrays
 * may be read uninitialized.
 */
static inline void *test_alloc(size_t const size) {
    void *const p = calloc(size > 0 ? size : 1, 1);
    if (p == NULL) {
        abort();
    }
    return p;
}

/** @brief State of a xorshift generator with a fixed seed, so every run uses the same inputs. */
static uint32_t rng_state = 0x9e3779b9u;

/** @brief Restarts the generator from @p seed. */
static inline void rng_seed(uint32_t const seed) {
    rng_state = seed;
}

/** @brief Returns the next 32 random bits. */
static inline uint32_t rng_next(void) {
    uint32_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_state = x;
    return x;
}

/**
 * @brief Sets @p x to a random integer with @f$|x| < 2^{bits}@f$.
 *
 * @f$|x|@f$ and the sign are uniform, so 0 is twice as likely as any other value.
 */
static inline void rng_entry(mpz_t x, unsigned long const bits) {
    mpz_set_ui(x, 0);
    for (unsigned long k = 0; k < bits; k += 32) {
        mpz_mul_2exp(x, x, 32);
        mpz_add_ui(x, x, rng_next());
    }
    mpz_fdiv_r_2exp(x, x, bits);
    if (rng_next() & 1) {
        mpz_neg(x, x);
    }
}

/** @brief Returns @p m zero vectors of length @p n. */
static inline clll_vector *basis_new(unsigned long const m, unsigned long const n) {
    clll_vector *const b = test_alloc(m * sizeof(clll_vector));
    for (unsigned long i = 0; i < m; i++) {
        clll_vector_init(&b[i], n);
    }
    return b;
}

/** @brief Clears and frees @p m vectors from basis_new(). */
static inline void basis_free(clll_vector *const b, unsigned long const m) {
    for (unsigned long i = 0; i < m; i++) {
        clll_vector_clear(&b[i]);
    }
    free(b);
}

/** @brief Returns a copy of the @p m vectors of length @p n in @p b. */
static inline clll_vector *basis_copy(clll_vector const *const b, unsigned long const m, unsigned long const n) {
    clll_vector *const c = basis_new(m, n);
    for (unsigned long i = 0; i < m; i++) {
        clll_vector_set(&c[i], &b[i]);
    }
    return c;
}

#endif
