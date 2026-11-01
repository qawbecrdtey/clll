/**
 * @file
 * @brief Compares clll_gram_schmidt_orthogonalization() with Examples 2 and 3
 *        of the paper and with rational Gram-Schmidt on random bases.
 */

#include <clll.h>

#include "test_util.h"

#include <gmp.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static char const *const example1_f[3][3] = {
    {"1", "1894885908", "0"},
    {"0", "1", "1894885908"},
    {"0", "0", "2147483648"},
};
static char const *const example1_g[3][3] = {
    {"1", "1894885908", "0"},
    {"-1894885908/3590592604336984465", "1/3590592604336984465", "1894885908"},
    {"7710738904443408018070044672/12892355250319448667906759645314351761",
     "-4069236502255632384/12892355250319448667906759645314351761",
     "2147483648/12892355250319448667906759645314351761"},
};

static char const *const example3_f[3][3] = {
    {"-3", "17", "4"},
    {"-8480", "-811", "-2908"},
    {"1290", "3351", "-13268"},
};
static char const *const example3_g[3][3] = {
    {"-3", "17", "4"},
    {"-2662657/314", "-255011/314", "-456598/157"},
    {"99196564668416/25441719249", "91577292685312/25441719249", "-314806070411264/25441719249"},
};

static bool check_paper_example(char const *const name, char const *const f_str[3][3],
                                char const *const g_str[3][3]) {
    clll_vector *const f = basis_new(3, 3);
    for (unsigned long i = 0; i < 3; i++) {
        for (unsigned long l = 0; l < 3; l++) {
            mpz_set_str(f[i].v[l], f_str[i][l], 10);
        }
    }
    clll_vector *const mu = basis_new(3, 3);
    clll_gram_schmidt_orthogonalization(mu, f, 3);

    printf("mu~ of %s (row i: mu~_{i,0} ... mu~_{i,i} = d_{i+1}):\n", name);
    for (unsigned long i = 0; i < 3; i++) {
        printf(" ");
        for (unsigned long j = 0; j <= i; j++) {
            gmp_printf(" %Zd", mu[i].v[j]);
        }
        printf("\n");
    }

    mpq_t g[3][3], c, p, expected;
    mpq_init(c);
    mpq_init(p);
    mpq_init(expected);
    bool ok = true;
    for (unsigned long i = 0; i < 3; i++) {
        for (unsigned long l = 0; l < 3; l++) {
            mpq_init(g[i][l]);
            mpq_set_z(g[i][l], f[i].v[l]);
        }
        for (unsigned long j = 0; j < i; j++) {
            mpq_set_num(c, mu[i].v[j]);
            mpq_set_den(c, mu[j].v[j]);
            mpq_canonicalize(c);
            for (unsigned long l = 0; l < 3; l++) {
                mpq_mul(p, c, g[j][l]);
                mpq_sub(g[i][l], g[i][l], p);
            }
        }
        for (unsigned long l = 0; l < 3; l++) {
            mpq_set_str(expected, g_str[i][l], 10);
            mpq_canonicalize(expected);
            if (!mpq_equal(g[i][l], expected)) {
                printf("MISMATCH: %s: entry %lu of g_%lu is not %s\n", name, l, i, g_str[i][l]);
                ok = false;
            }
        }
    }

    mpz_t d3;
    mpz_init(d3);
    mpz_set_ui(d3, 1);
    mpz_mul_2exp(d3, d3, 62);
    if (mpz_cmp(mu[2].v[2], d3) != 0) {
        printf("MISMATCH: %s: d_3 is not 2^62\n", name);
        ok = false;
    }
    if (ok) {
        printf("OK: %s: the Gram-Schmidt vectors match the paper exactly, and d_3 = 2^62\n", name);
    }

    mpz_clear(d3);
    for (unsigned long i = 0; i < 3; i++) {
        for (unsigned long l = 0; l < 3; l++) {
            mpq_clear(g[i][l]);
        }
    }
    mpq_clear(c);
    mpq_clear(p);
    mpq_clear(expected);
    basis_free(mu, 3);
    basis_free(f, 3);
    return ok;
}

static unsigned long sentinel(unsigned long const i, unsigned long const j) {
    return 1000003 + 100 * i + j;
}

typedef struct {
    unsigned long bases;
    unsigned long skipped;
    unsigned long entries;
} tally;

static bool check_random(tally *const t, unsigned long const m, unsigned long const n, unsigned long const bits) {
    clll_vector *const f = basis_new(m, n);
    for (unsigned long i = 0; i < m; i++) {
        for (unsigned long l = 0; l < n; l++) {
            rng_entry(f[i].v[l], bits);
        }
    }

    mpq_t *const g = test_alloc((size_t)m * n * sizeof(mpq_t));
    mpq_t *const mu_q = test_alloc((size_t)m * m * sizeof(mpq_t));
    mpq_t *const d = test_alloc((m + 1) * sizeof(mpq_t));
    mpq_t ip, p, b;
    mpq_init(ip);
    mpq_init(p);
    mpq_init(b);
    for (size_t k = 0; k < (size_t)m * n; k++) {
        mpq_init(g[k]);
    }
    for (size_t k = 0; k < (size_t)m * m; k++) {
        mpq_init(mu_q[k]);
    }
    for (unsigned long k = 0; k <= m; k++) {
        mpq_init(d[k]);
    }

    mpq_set_ui(d[0], 1, 1);
    bool dependent = false;
    for (unsigned long i = 0; i < m && !dependent; i++) {
        for (unsigned long l = 0; l < n; l++) {
            mpq_set_z(g[i * n + l], f[i].v[l]);
        }
        for (unsigned long j = 0; j < i; j++) {
            mpq_set_ui(ip, 0, 1);
            for (unsigned long l = 0; l < n; l++) {
                mpq_set_z(p, f[i].v[l]);
                mpq_mul(p, p, g[j * n + l]);
                mpq_add(ip, ip, p);
            }
            mpq_div(b, d[j + 1], d[j]);
            mpq_div(mu_q[i * m + j], ip, b);
            for (unsigned long l = 0; l < n; l++) {
                mpq_mul(p, mu_q[i * m + j], g[j * n + l]);
                mpq_sub(g[i * n + l], g[i * n + l], p);
            }
        }
        mpq_set_ui(b, 0, 1);
        for (unsigned long l = 0; l < n; l++) {
            mpq_mul(p, g[i * n + l], g[i * n + l]);
            mpq_add(b, b, p);
        }
        dependent = mpq_sgn(b) == 0;
        mpq_mul(d[i + 1], d[i], b);
    }

    bool ok = true;
    if (dependent) {
        t->skipped += 1;
    } else {
        clll_vector *const mu = basis_new(m, n);
        for (unsigned long i = 0; i < m; i++) {
            for (unsigned long j = 0; j < n; j++) {
                mpz_set_ui(mu[i].v[j], j <= i ? 77 : sentinel(i, j));
            }
        }
        clll_gram_schmidt_orthogonalization(mu, f, m);

        for (unsigned long i = 0; i < m; i++) {
            for (unsigned long j = 0; j < i; j++) {
                mpq_mul(p, d[j + 1], mu_q[i * m + j]);
                if (mpq_cmp_z(p, mu[i].v[j]) != 0) {
                    printf("MISMATCH: %lu x %lu, %lu-bit: mu~_{%lu,%lu} is not d_%lu mu_{%lu,%lu}\n", m, n, bits, i,
                           j, j + 1, i, j);
                    ok = false;
                }
            }
            if (mpq_cmp_z(d[i + 1], mu[i].v[i]) != 0) {
                printf("MISMATCH: %lu x %lu, %lu-bit: mu~_{%lu,%lu} is not d_%lu\n", m, n, bits, i, i, i + 1);
                ok = false;
            }
            for (unsigned long j = i + 1; j < n; j++) {
                if (mpz_cmp_ui(mu[i].v[j], sentinel(i, j)) != 0) {
                    printf("MISMATCH: %lu x %lu, %lu-bit: entry %lu of row %lu, above the diagonal, was written\n",
                           m, n, bits, j, i);
                    ok = false;
                }
            }
            t->entries += n;
        }
        basis_free(mu, m);
        t->bases += 1;
    }

    for (size_t k = 0; k < (size_t)m * n; k++) {
        mpq_clear(g[k]);
    }
    for (size_t k = 0; k < (size_t)m * m; k++) {
        mpq_clear(mu_q[k]);
    }
    for (unsigned long k = 0; k <= m; k++) {
        mpq_clear(d[k]);
    }
    free(g);
    free(mu_q);
    free(d);
    mpq_clear(ip);
    mpq_clear(p);
    mpq_clear(b);
    basis_free(f, m);
    return ok;
}

int main(void) {
    bool ok = true;

    ok = check_paper_example("Example 1 (paper Example 2)", example1_f, example1_g) && ok;
    ok = check_paper_example("Example 3 (paper Example 3)", example3_f, example3_g) && ok;

    rng_seed(0x2545f491u);
    tally t = {0, 0, 0};
    bool random_ok = true;
    for (unsigned long m = 1; m <= 6; m++) {
        for (unsigned long n = m; n <= m + 2; n++) {
            for (unsigned long draw = 0; draw < 4; draw++) {
                random_ok = check_random(&t, m, n, 3) && random_ok;
                random_ok = check_random(&t, m, n, 12) && random_ok;
                random_ok = check_random(&t, m, n, 60) && random_ok;
            }
        }
    }
    if (random_ok) {
        printf("OK: %lu random bases (m = 1..6, n = m..m+2): all %lu entries agree with rational Gram-Schmidt"
               " (dependent draws skipped: %lu)\n",
               t.bases, t.entries, t.skipped);
    }
    ok = random_ok && ok;

    {
        clll_vector *const f = basis_new(1, 3);
        clll_vector *const mu = basis_new(1, 3);
        for (unsigned long j = 0; j < 3; j++) {
            mpz_set_ui(f[0].v[j], j + 1);
            mpz_set_ui(mu[0].v[j], sentinel(0, j));
        }
        clll_gram_schmidt_orthogonalization(mu, f, 0);
        bool untouched = true;
        for (unsigned long j = 0; j < 3; j++) {
            untouched = untouched && mpz_cmp_ui(mu[0].v[j], sentinel(0, j)) == 0;
        }
        printf("%s: m = 0 leaves mu untouched\n", untouched ? "OK" : "MISMATCH");
        ok = untouched && ok;
        basis_free(mu, 1);
        basis_free(f, 1);
    }

    return ok ? 0 : 1;
}
