#include <clll/gram_schmidt_orthogonalization.h>

#include <clll/vector.h>

#include "assert_and_assume.h"
#include "gram_schmidt_orthogonalization_row.h"

#include <gmp.h>

#include <stddef.h>

void clll_gram_schmidt_orthogonalization_row(
        mpz_t *const mu,
        clll_vector const *const fs,
        unsigned long const i,
        mpz_t sigma,
        mpz_t t) {
    ASSERT_AND_ASSUME(i < fs[0].len);
    ASSERT_AND_ASSUME(fs[i].len == fs[0].len);

    clll_vector_inner_product(mu[mu_index(i, 0)], &fs[i], &fs[0]);
    for (unsigned long j = 1; j <= i; j++) {
        mpz_mul(sigma, mu[mu_index(i, 0)], mu[mu_index(j, 0)]);
        for (unsigned long l = 1; l < j; l++) {
            mpz_mul(t, mu[mu_index(l, l)], sigma);
            mpz_addmul(t, mu[mu_index(i, l)], mu[mu_index(j, l)]);
            ASSERT_AND_ASSUME(mpz_divisible_p(t, mu[mu_index(l - 1, l - 1)]));
            mpz_divexact(sigma, t, mu[mu_index(l - 1, l - 1)]);
        }
        clll_vector_inner_product(t, &fs[i], &fs[j]);
        mpz_mul(mu[mu_index(i, j)], mu[mu_index(j - 1, j - 1)], t);
        mpz_sub(mu[mu_index(i, j)], mu[mu_index(i, j)], sigma);
    }
}

void clll_gram_schmidt_orthogonalization(
        clll_vector *const mu,
        clll_vector const *const fs,
        unsigned long const m) {
    ASSERT_AND_ASSUME(m == 0 || m <= fs[0].len);

    if (m == 0) {
        return;
    }

    void *(*allocate)(size_t);
    void (*deallocate)(void *, size_t);
    mp_get_memory_functions(&allocate, NULL, &deallocate);

    size_t const entries = d_index(m) + 1;
    mpz_t *const packed = allocate(entries * sizeof(mpz_t));
    for (size_t k = 0; k < entries; k++) {
        mpz_init(packed[k]);
    }
    for (unsigned long i = 0; i < m; i++) {
        ASSERT_AND_ASSUME(mu[i].len == fs[0].len);
        for (unsigned long j = 0; j <= i; j++) {
            mpz_swap(packed[mu_index(i, j)], mu[i].v[j]);
        }
    }

    mpz_t sigma, t;
    mpz_inits(sigma, t, NULL);
    for (unsigned long i = 0; i < m; i++) {
        clll_gram_schmidt_orthogonalization_row(packed, fs, i, sigma, t);
    }
    mpz_clears(sigma, t, NULL);

    for (unsigned long i = 0; i < m; i++) {
        for (unsigned long j = 0; j <= i; j++) {
            mpz_swap(packed[mu_index(i, j)], mu[i].v[j]);
        }
    }
    for (size_t k = 0; k < entries; k++) {
        mpz_clear(packed[k]);
    }
    deallocate(packed, entries * sizeof(mpz_t));
}
