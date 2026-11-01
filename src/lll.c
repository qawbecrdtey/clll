#include <clll/lll.h>

#include <clll/vector.h>

#include "assert_and_assume.h"
#include "gram_schmidt_orthogonalization_row.h"

#include <gmp.h>

#include <limits.h>
#include <stdbool.h>
#include <stddef.h>

static void submul_entries(mpz_t *const rop, mpz_t const c, mpz_t *const v, unsigned long const len) {
    mp_limb_t const c0 = mpz_getlimbn(c, 0);
    if (mpz_size(c) == 1 && c0 <= ULONG_MAX) {
        unsigned long const u = (unsigned long)c0;
        if (mpz_sgn(c) > 0) {
            for (unsigned long l = 0; l < len; l++) {
                mpz_submul_ui(rop[l], v[l], u);
            }
        }
        else {
            for (unsigned long l = 0; l < len; l++) {
                mpz_addmul_ui(rop[l], v[l], u);
            }
        }
    }
    else {
        for (unsigned long l = 0; l < len; l++) {
            mpz_submul(rop[l], c, v[l]);
        }
    }
}

void clll_lll(clll_vector *const basis, unsigned long const m, mpq_t const alpha_in) {
    mpq_t alpha;
    mpq_init(alpha);
    mpq_set(alpha, alpha_in);
    mpq_canonicalize(alpha);
    ASSERT_AND_ASSUME(mpq_cmp_ui(alpha, 4, 3) >= 0);

    if (m == 0) {
        mpq_clear(alpha);
        return;
    }

    void *(*allocate)(size_t);
    void (*deallocate)(void *, size_t);
    mp_get_memory_functions(&allocate, NULL, &deallocate);

    // Line 1
    ASSERT_AND_ASSUME(m <= basis[0].len);
    size_t const entries = d_index(m) + 1;
    mpz_t *const mu = allocate(entries * sizeof(mpz_t));
    for (size_t k = 0; k < entries; k++) {
        mpz_init(mu[k]);
    }
    mpz_set_ui(mu[d_index(0)], 1);
    mpz_srcptr const num = mpq_numref(alpha);
    mpz_srcptr const den = mpq_denref(alpha);
    bool const den_one = mpz_size(den) == 1 && mpz_getlimbn(den, 0) == 1;

    mpz_t *const sq = allocate((m + 1) * sizeof(mpz_t));
    mpz_t *const pr = allocate((m + 1) * sizeof(mpz_t));
    bool *const sq_ok = allocate((m + 1) * sizeof(bool));
    bool *const pr_ok = allocate((m + 1) * sizeof(bool));
    for (unsigned long k = 0; k <= m; k++) {
        mpz_init(sq[k]);
        mpz_init(pr[k]);
        sq_ok[k] = false;
        pr_ok[k] = false;
    }

    mpz_t *const half = allocate((m + 1) * sizeof(mpz_t));
    bool *const half_ok = allocate((m + 1) * sizeof(bool));
    for (unsigned long k = 0; k <= m; k++) {
        mpz_init(half[k]);
        half_ok[k] = false;
    }

    mpz_t c, rhs, a, b;
    mpz_inits(c, rhs, a, b, NULL);

    unsigned long i = 0;
    unsigned long rows = 0;
    bool upw = true;
    // Line 2
    while (i < m) {
        if (i == rows) {
            clll_gram_schmidt_orthogonalization_row(mu, basis, i, a, b);
            rows += 1;
        }

        if (upw) {
            // Lines 3-4: size reduction.
            for (unsigned long j = i; j-- > 0;) {
                mpz_ptr const x = mu[mu_index(i, j)];
                mpz_srcptr const dj = mu[d_index(j + 1)];
                if (!half_ok[j + 1]) {
                    mpz_fdiv_q_2exp(half[j + 1], dj, 1);
                    half_ok[j + 1] = true;
                }
                // round(mu_{i,j}) is 0 iff -h <= mu~_{i,j} < h + (d_{j+1} mod 2), with h = half[j + 1].
                int const cmp = mpz_cmpabs(x, half[j + 1]);
                if (cmp < 0 || (cmp == 0 && (mpz_sgn(x) < 0 || mpz_odd_p(dj)))) {
                    continue;
                }
                mpz_add(rhs, x, half[j + 1]);
                mpz_fdiv_qr(c, x, rhs, dj);
                mpz_sub(x, x, half[j + 1]);
                submul_entries(basis[i].v, c, basis[j].v, basis[i].len);
                submul_entries(mu + mu_index(i, 0), c, mu + mu_index(j, 0), j);
            }
        }

        // Line 5: Lovasz test.
        bool violated = false;
        if (i > 0) {
            if (!sq_ok[i]) {
                mpz_mul(sq[i], mu[d_index(i)], mu[d_index(i)]);
                if (!den_one) {
                    mpz_mul(sq[i], sq[i], den);
                }
                sq_ok[i] = true;
            }
            if (!pr_ok[i]) {
                mpz_mul(pr[i], mu[d_index(i - 1)], mu[d_index(i + 1)]);
                pr_ok[i] = true;
            }
            mpz_mul(rhs, pr[i], num);
            violated = mpz_cmp(sq[i], rhs) > 0;
        }

        if (violated) {
            // Line 6: swap.
            for (unsigned long j = 0; j + 1 < i; j++) {
                mpz_swap(mu[mu_index(i - 1, j)], mu[mu_index(i, j)]);
            }
            for (unsigned long k = i + 1; k < rows; k++) {
                mpz_mul(a, mu[mu_index(i, i - 1)], mu[mu_index(k, i - 1)]);
                mpz_addmul(a, mu[mu_index(k, i)], mu[d_index(i - 1)]);
                mpz_mul(b, mu[d_index(i + 1)], mu[mu_index(k, i - 1)]);
                mpz_submul(b, mu[mu_index(i, i - 1)], mu[mu_index(k, i)]);
                ASSERT_AND_ASSUME(mpz_divisible_p(a, mu[d_index(i)]));
                mpz_divexact(mu[mu_index(k, i - 1)], a, mu[d_index(i)]);
                ASSERT_AND_ASSUME(mpz_divisible_p(b, mu[d_index(i)]));
                mpz_divexact(mu[mu_index(k, i)], b, mu[d_index(i)]);
            }
            mpz_mul(a, mu[mu_index(i, i - 1)], mu[mu_index(i, i - 1)]);
            mpz_add(a, a, pr[i]);
            ASSERT_AND_ASSUME(mpz_divisible_p(a, mu[d_index(i)]));
            mpz_divexact(mu[d_index(i)], a, mu[d_index(i)]);
            sq_ok[i] = false;
            pr_ok[i - 1] = false;
            pr_ok[i + 1] = false;
            half_ok[i] = false;

            clll_vector_swap(&basis[i - 1], &basis[i]);
            i -= 1;
            upw = false;
        }
        else {
            // Line 7
            i += 1;
            upw = true;
        }
    }

    // Line 8
    mpz_clears(c, rhs, a, b, NULL);
    for (unsigned long k = 0; k <= m; k++) {
        mpz_clear(half[k]);
    }
    deallocate(half, (m + 1) * sizeof(mpz_t));
    deallocate(half_ok, (m + 1) * sizeof(bool));
    for (unsigned long k = 0; k <= m; k++) {
        mpz_clear(sq[k]);
        mpz_clear(pr[k]);
    }
    deallocate(sq, (m + 1) * sizeof(mpz_t));
    deallocate(pr, (m + 1) * sizeof(mpz_t));
    deallocate(sq_ok, (m + 1) * sizeof(bool));
    deallocate(pr_ok, (m + 1) * sizeof(bool));
    for (size_t k = 0; k < entries; k++) {
        mpz_clear(mu[k]);
    }
    deallocate(mu, entries * sizeof(mpz_t));
    mpq_clear(alpha);
}
