/**
 * @file
 * @brief Checks that a clll_lll() result is a reduced basis of the input
 *        lattice.
 *
 * Any correct reduced basis passes. Lattice membership is decided by exact
 * rational elimination, without using the library.
 */

#ifndef CLLL_TESTS_LLL_CONDITIONS_H
#define CLLL_TESTS_LLL_CONDITIONS_H

#include "test_util.h"

#include <clll.h>

#include <gmp.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief The values of @f$\alpha@f$ every test uses, as {numerator, denominator}.
 *
 * 4/3 is the smallest allowed, 7/5 is close to fplll's default of 1250/887, and
 * 3/2 is the paper's.
 */
static unsigned long const test_alphas[][2] = {{4, 3}, {7, 5}, {3, 2}, {2, 1}, {100, 3}};

/**
 * @brief Returns the rows of clll_gram_schmidt_orthogonalization() for the
 *        @p m vectors of length @p n in @p b.
 */
static inline clll_vector *gso_rows(clll_vector const *const b, unsigned long const m, unsigned long const n) {
    clll_vector *const mu = basis_new(m, n);
    clll_gram_schmidt_orthogonalization(mu, b, m);
    return mu;
}

/** @brief Sets @p d to @f$d_k@f$, read from the rows returned by gso_rows(). */
static inline void gram_det(mpz_t d, clll_vector const *const mu, unsigned long const k) {
    if (k == 0) {
        mpz_set_ui(d, 1);
    } else {
        mpz_set(d, mu[k - 1].v[k - 1]);
    }
}

/**
 * @brief Runs Gaussian elimination over @f$\mathbb{Q}@f$ on @f$[G \mid R]@f$.
 *
 * @f$G_{k,l} = f_k \cdot f_l@f$ is the Gram matrix of @p f and
 * @f$R_{k,s} = v_s \cdot f_k@f$. The @f$m \times (m + t)@f$ result is stored
 * row by row in @p a, which the caller clears with gram_clear().
 *
 * @return Whether @p f is linearly independent.
 */
static inline bool gram_eliminate(mpq_t *const a, clll_vector const *const f, unsigned long const m,
                           clll_vector const *const v, unsigned long const t) {
    unsigned long const w = m + t;
    mpz_t ip;
    mpz_init(ip);
    for (unsigned long k = 0; k < m; k++) {
        for (unsigned long c = 0; c < w; c++) {
            mpq_init(a[k * w + c]);
            clll_vector_inner_product(ip, &f[k], c < m ? &f[c] : &v[c - m]);
            mpq_set_z(a[k * w + c], ip);
        }
    }
    mpz_clear(ip);

    mpq_t q, p;
    mpq_init(q);
    mpq_init(p);
    bool independent = true;
    for (unsigned long r = 0; r < m; r++) {
        if (mpq_sgn(a[r * w + r]) == 0) {
            independent = false;
            break;
        }
        for (unsigned long k = r + 1; k < m; k++) {
            mpq_div(q, a[k * w + r], a[r * w + r]);
            for (unsigned long c = r; c < w; c++) {
                mpq_mul(p, q, a[r * w + c]);
                mpq_sub(a[k * w + c], a[k * w + c], p);
            }
        }
    }
    mpq_clear(q);
    mpq_clear(p);
    return independent;
}

/** @brief Clears and frees the matrix filled by gram_eliminate(). */
static inline void gram_clear(mpq_t *const a, unsigned long const m, unsigned long const t) {
    for (size_t k = 0; k < (size_t)m * (m + t); k++) {
        mpq_clear(a[k]);
    }
    free(a);
}

/** @brief Returns whether the @p m vectors in @p f are linearly independent. */
static inline bool independent(clll_vector const *const f, unsigned long const m) {
    mpq_t *const a = test_alloc((size_t)m * m * sizeof(mpq_t));
    bool const ok = gram_eliminate(a, f, m, NULL, 0);
    gram_clear(a, m, 0);
    return ok;
}

/**
 * @brief Returns whether each of v[0], ..., v[t-1] is an integer combination
 *        of the independent f[0], ..., f[m-1].
 */
static inline bool in_lattice(clll_vector const *const v, unsigned long const t, clll_vector const *const f,
                       unsigned long const m, unsigned long const n) {
    unsigned long const w = m + t;
    mpq_t *const a = test_alloc((size_t)m * w * sizeof(mpq_t));
    bool ok = gram_eliminate(a, f, m, v, t);

    mpq_t *const x = test_alloc(m * sizeof(mpq_t));
    for (unsigned long k = 0; k < m; k++) {
        mpq_init(x[k]);
    }
    mpq_t p;
    mpq_init(p);
    mpz_t sum;
    mpz_init(sum);
    for (unsigned long s = 0; s < t && ok; s++) {
        for (unsigned long k = m; k-- > 0;) {
            mpq_set(x[k], a[k * w + m + s]);
            for (unsigned long c = k + 1; c < m; c++) {
                mpq_mul(p, a[k * w + c], x[c]);
                mpq_sub(x[k], x[k], p);
            }
            mpq_div(x[k], x[k], a[k * w + k]);
            ok = ok && mpz_cmp_ui(mpq_denref(x[k]), 1) == 0;
        }
        for (unsigned long l = 0; l < n && ok; l++) {
            mpz_set_ui(sum, 0);
            for (unsigned long k = 0; k < m; k++) {
                mpz_addmul(sum, mpq_numref(x[k]), f[k].v[l]);
            }
            ok = mpz_cmp(sum, v[s].v[l]) == 0;
        }
    }
    mpz_clear(sum);
    mpq_clear(p);
    for (unsigned long k = 0; k < m; k++) {
        mpq_clear(x[k]);
    }
    free(x);
    gram_clear(a, m, t);
    return ok;
}

/** @brief Returns whether @f$|\mu_{i,j}| \le 1/2@f$ for all @f$j < i@f$. */
static inline bool size_reduced(clll_vector const *const mu, unsigned long const m) {
    mpz_t t;
    mpz_init(t);
    bool ok = true;
    for (unsigned long i = 0; i < m && ok; i++) {
        for (unsigned long j = 0; j < i && ok; j++) {
            mpz_mul_2exp(t, mu[i].v[j], 1);
            ok = mpz_cmpabs(t, mu[j].v[j]) <= 0;
        }
    }
    mpz_clear(t);
    return ok;
}

/** @brief Returns whether @f$\|g_{i-1}\|^2 \le \alpha \|g_i\|^2@f$ for all @f$0 < i < m@f$. */
static inline bool lovasz_holds(clll_vector const *const mu, unsigned long const m, mpq_t const alpha) {
    mpz_t lhs, rhs, d;
    mpz_init(lhs);
    mpz_init(rhs);
    mpz_init(d);
    bool ok = true;
    for (unsigned long i = 1; i < m && ok; i++) {
        gram_det(d, mu, i);
        mpz_mul(lhs, d, d);
        mpz_mul(lhs, lhs, mpq_denref(alpha));
        gram_det(d, mu, i - 1);
        mpz_mul(rhs, d, mu[i].v[i]);
        mpz_mul(rhs, rhs, mpq_numref(alpha));
        ok = mpz_cmp(lhs, rhs) <= 0;
    }
    mpz_clear(lhs);
    mpz_clear(rhs);
    mpz_clear(d);
    return ok;
}

/**
 * @brief Checks that @p out is a reduced basis of the lattice spanned by @p in,
 *        as lll.h defines it.
 *
 * Prints a line for every condition that fails.
 *
 * @return Whether every condition holds.
 */
static inline bool check_reduced(char const *const label, clll_vector const *const in, clll_vector const *const out,
                          unsigned long const m, unsigned long const n, mpq_t const alpha) {
    if (m == 0) {
        return true;
    }
    for (unsigned long i = 0; i < m; i++) {
        if (out[i].len != n) {
            printf("MISMATCH: %s: vector %lu has length %lu, not %lu\n", label, i, out[i].len, n);
            return false;
        }
    }
    if (!independent(out, m)) {
        printf("MISMATCH: %s: the output is linearly dependent\n", label);
        return false;
    }

    bool ok = true;
    clll_vector *const mu_out = gso_rows(out, m, n);
    clll_vector *const mu_in = gso_rows(in, m, n);
    if (!size_reduced(mu_out, m)) {
        printf("MISMATCH: %s: not size-reduced, some |mu_{i,j}| > 1/2\n", label);
        ok = false;
    }
    if (!lovasz_holds(mu_out, m, alpha)) {
        printf("MISMATCH: %s: some ||g_{i-1}||^2 > alpha ||g_i||^2\n", label);
        ok = false;
    }
    if (mpz_cmp(mu_out[m - 1].v[m - 1], mu_in[m - 1].v[m - 1]) != 0) {
        printf("MISMATCH: %s: d_m changed, so the lattice did\n", label);
        ok = false;
    }
    if (!in_lattice(out, m, in, m, n)) {
        printf("MISMATCH: %s: an output vector is not in the input lattice\n", label);
        ok = false;
    }
    if (!in_lattice(in, m, out, m, n)) {
        printf("MISMATCH: %s: an input vector is not in the output lattice\n", label);
        ok = false;
    }
    basis_free(mu_in, m);
    basis_free(mu_out, m);
    return ok;
}

#endif
