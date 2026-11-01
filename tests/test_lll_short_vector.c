/**
 * @file
 * @brief Checks the short-vector bounds of a reduced basis on Example 1 of the
 *        paper and on subset-sum lattices.
 *
 * For a basis @f$b_0, \ldots, b_{m-1}@f$ reduced with respect to
 * @f$\alpha@f$, with @f$T = m(m-1)/2@f$:
 *
 * 1. @f$\|b_0\|^2 \le \alpha^{m-1} \|v\|^2@f$ for every nonzero lattice
 *    vector @f$v@f$,
 * 2. @f$\|b_0\|^{2m} \le \alpha^T d_m@f$, and
 * 3. @f$\prod_i \|b_i\|^2 \le \alpha^T d_m@f$.
 *
 * Whether the reduced basis contains a shortest vector, or exposes the planted
 * subset, is printed but not checked.
 */

#include <clll.h>

#include "lll_conditions.h"

#include <gmp.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

static void norm2(mpz_t r, clll_vector const *const x) {
    clll_vector_inner_product(r, x, x);
}

static bool first_vector_bound(clll_vector const *const b, unsigned long const m, mpq_t const alpha,
                               clll_vector const *const v) {
    mpz_t lhs, rhs, t;
    mpz_init(lhs);
    mpz_init(rhs);
    mpz_init(t);
    norm2(t, &b[0]);
    mpz_pow_ui(lhs, mpq_denref(alpha), m - 1);
    mpz_mul(lhs, lhs, t);
    norm2(t, v);
    mpz_pow_ui(rhs, mpq_numref(alpha), m - 1);
    mpz_mul(rhs, rhs, t);
    bool const ok = mpz_cmp(lhs, rhs) <= 0;
    mpz_clear(lhs);
    mpz_clear(rhs);
    mpz_clear(t);
    return ok;
}

static bool determinant_bounds(char const *const label, clll_vector const *const b, unsigned long const m,
                               unsigned long const n, mpq_t const alpha) {
    unsigned long const t = m * (m - 1) / 2;
    clll_vector *const mu = gso_rows(b, m, n);
    mpz_t rhs, lhs, den_t, x;
    mpz_init(rhs);
    mpz_init(lhs);
    mpz_init(den_t);
    mpz_init(x);
    mpz_pow_ui(rhs, mpq_numref(alpha), t);
    mpz_mul(rhs, rhs, mu[m - 1].v[m - 1]);
    mpz_pow_ui(den_t, mpq_denref(alpha), t);

    bool ok = true;
    norm2(x, &b[0]);
    mpz_pow_ui(lhs, x, m);
    mpz_mul(lhs, lhs, den_t);
    if (mpz_cmp(lhs, rhs) > 0) {
        printf("MISMATCH: %s: ||b_0||^{2m} > alpha^T d_m\n", label);
        ok = false;
    }

    mpz_set(lhs, den_t);
    for (unsigned long i = 0; i < m; i++) {
        norm2(x, &b[i]);
        mpz_mul(lhs, lhs, x);
    }
    if (mpz_cmp(lhs, rhs) > 0) {
        printf("MISMATCH: %s: prod_i ||b_i||^2 > alpha^T d_m\n", label);
        ok = false;
    }

    mpz_clear(rhs);
    mpz_clear(lhs);
    mpz_clear(den_t);
    mpz_clear(x);
    basis_free(mu, m);
    return ok;
}

static bool reduce_and_check(char const *const label, clll_vector const *const in, unsigned long const m,
                             unsigned long const n, mpq_t const alpha, clll_vector const *const v,
                             clll_vector **const out) {
    if (!in_lattice(v, 1, in, m, n)) {
        printf("MISMATCH: %s: the planted vector is not in the lattice (test setup)\n", label);
        return false;
    }
    *out = basis_copy(in, m, n);
    clll_lll(*out, m, alpha);
    if (!check_reduced(label, in, *out, m, n, alpha)) {
        return false;
    }

    bool ok = true;
    for (unsigned long i = 0; i < m; i++) {
        if (!first_vector_bound(*out, m, alpha, &in[i])) {
            printf("MISMATCH: %s: ||b_0||^2 > alpha^{m-1} ||f_%lu||^2\n", label, i);
            ok = false;
        }
    }
    if (!first_vector_bound(*out, m, alpha, v)) {
        printf("MISMATCH: %s: ||b_0||^2 > alpha^{m-1} ||v||^2 for the planted v\n", label);
        ok = false;
    }
    return determinant_bounds(label, *out, m, n, alpha) && ok;
}

static bool example1(unsigned long const num, unsigned long const den) {
    clll_vector *const f = basis_new(3, 3);
    mpz_set_ui(f[0].v[0], 1);
    mpz_set_ui(f[0].v[1], 1894885908);
    mpz_set_ui(f[1].v[1], 1);
    mpz_set_ui(f[1].v[2], 1894885908);
    mpz_set_str(f[2].v[2], "2147483648", 10);
    clll_vector *const v = basis_new(1, 3);
    mpz_set_si(v[0].v[0], -3);
    mpz_set_si(v[0].v[1], 17);
    mpz_set_si(v[0].v[2], 4);

    mpq_t alpha;
    mpq_init(alpha);
    mpq_set_ui(alpha, num, den);
    char label[64];
    snprintf(label, sizeof label, "Example 1, alpha = %lu/%lu", num, den);
    clll_vector *b = NULL;
    bool const ok = reduce_and_check(label, f, 3, 3, alpha, v, &b);
    if (ok) {
        mpz_t x;
        mpz_init(x);
        norm2(x, &b[0]);
        gmp_printf("OK: %s: conditions and bounds (1)-(3) hold; b_0 = (%Zd, %Zd, %Zd), ||b_0||^2 = %Zd%s\n", label,
                   b[0].v[0], b[0].v[1], b[0].v[2], x,
                   mpz_cmp_ui(x, 314) == 0 ? " (a shortest vector)" : "");
        mpz_clear(x);
    }
    if (b != NULL) {
        basis_free(b, 3);
    }
    mpq_clear(alpha);
    basis_free(v, 1);
    basis_free(f, 3);
    return ok;
}

static bool exposes_subset(clll_vector const *const b_r, clll_vector const *const f, unsigned long const k) {
    if (mpz_sgn(b_r->v[k]) != 0) {
        return false;
    }
    mpz_t sum;
    mpz_init(sum);
    int sign = 0;
    bool shape = true;
    for (unsigned long i = 0; i < k && shape; i++) {
        int const e = mpz_sgn(b_r->v[i]);
        if (e != 0) {
            shape = mpz_cmpabs_ui(b_r->v[i], 1) == 0 && (sign == 0 || sign == e);
            sign = e;
            mpz_add(sum, sum, f[i].v[k]);
        }
    }
    bool const found = shape && sign != 0 && mpz_cmp(sum, f[k].v[k]) == 0;
    mpz_clear(sum);
    return found;
}

static bool subset_sum(unsigned long const k, unsigned long const draw, bool const show) {
    unsigned long const m = k + 1;
    clll_vector *const f = basis_new(m, m);
    clll_vector *const v = basis_new(1, m);

    mpz_t s;
    mpz_init(s);
    bool any = false;
    for (unsigned long i = 0; i < k; i++) {
        mpz_set_ui(f[i].v[i], 1);
        rng_entry(f[i].v[k], 40);
        mpz_abs(f[i].v[k], f[i].v[k]);
        mpz_setbit(f[i].v[k], 39);
        if ((rng_next() & 1) || (i + 1 == k && !any)) {
            mpz_set_ui(v[0].v[i], 1);
            mpz_add(s, s, f[i].v[k]);
            any = true;
        }
    }
    if (show) {
        printf("Subset-sum instance with %lu weights:", k);
        for (unsigned long i = 0; i < k; i++) {
            gmp_printf(" %Zd", f[i].v[k]);
        }
        gmp_printf("; target %Zd\n", s);
    }
    mpz_set(f[k].v[k], s);
    for (unsigned long i = 0; i < m; i++) {
        mpz_mul_2exp(f[i].v[k], f[i].v[k], 20);
    }

    bool ok = true;
    mpq_t alpha;
    mpq_init(alpha);
    for (size_t a = 0; a < sizeof test_alphas / sizeof test_alphas[0]; a++) {
        mpq_set_ui(alpha, test_alphas[a][0], test_alphas[a][1]);
        char label[80];
        snprintf(label, sizeof label, "subset sum, %lu weights, draw %lu, alpha = %lu/%lu", k, draw, test_alphas[a][0],
                 test_alphas[a][1]);
        clll_vector *b = NULL;
        bool const this_ok = reduce_and_check(label, f, m, m, alpha, v, &b);
        if (this_ok) {
            unsigned long r = 0;
            while (r < m && !exposes_subset(&b[r], f, k)) {
                r++;
            }
            if (r < m) {
                printf("OK: %s: conditions and bounds (1)-(3) hold; b_%lu exposes a subset with the target sum: {",
                       label, r);
                bool first = true;
                for (unsigned long i = 0; i < k; i++) {
                    if (mpz_sgn(b[r].v[i]) != 0) {
                        printf(first ? "%lu" : ", %lu", i);
                        first = false;
                    }
                }
                printf("}\n");
            } else {
                printf("OK: %s: conditions and bounds (1)-(3) hold; no reduced vector exposes a subset\n", label);
            }
        }
        if (b != NULL) {
            basis_free(b, m);
        }
        ok = this_ok && ok;
    }
    mpq_clear(alpha);
    mpz_clear(s);
    basis_free(v, 1);
    basis_free(f, m);
    return ok;
}

int main(void) {
    bool ok = true;

    for (size_t a = 0; a < sizeof test_alphas / sizeof test_alphas[0]; a++) {
        ok = example1(test_alphas[a][0], test_alphas[a][1]) && ok;
    }

    for (unsigned long k = 6; k <= 14; k += 2) {
        for (unsigned long draw = 0; draw < 2; draw++) {
            ok = subset_sum(k, draw, k == 12 && draw == 0) && ok;
        }
    }

    return ok ? 0 : 1;
}
