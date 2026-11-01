/**
 * @file
 * @brief Checks that clll_lll() returns a reduced basis of the same lattice
 *        for random, knapsack, identity and large-multiplier inputs.
 */

#include <clll.h>

#include "lll_conditions.h"

#include <gmp.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

typedef struct {
    unsigned long reductions;
    unsigned long skipped;
    bool ok;
} tally;

static void reduce_and_check(tally *const t, char const *const name, clll_vector const *const in,
                             unsigned long const m, unsigned long const n) {
    char label[128];
    mpq_t alpha;
    mpq_init(alpha);
    for (size_t a = 0; a < sizeof test_alphas / sizeof test_alphas[0]; a++) {
        mpq_set_ui(alpha, test_alphas[a][0], test_alphas[a][1]);
        clll_vector *const out = basis_copy(in, m, n);
        clll_lll(out, m, alpha);
        snprintf(label, sizeof label, "%s, alpha = %lu/%lu", name, test_alphas[a][0], test_alphas[a][1]);
        t->ok = check_reduced(label, in, out, m, n, alpha) && t->ok;
        basis_free(out, m);
        t->reductions += 1;
    }
    mpq_clear(alpha);
}

static void random_basis(tally *const t, unsigned long const m, unsigned long const n, unsigned long const bits,
                         unsigned long const draw) {
    clll_vector *const in = basis_new(m, n);
    for (unsigned long i = 0; i < m; i++) {
        for (unsigned long l = 0; l < n; l++) {
            rng_entry(in[i].v[l], bits);
        }
    }
    if (independent(in, m)) {
        char name[64];
        snprintf(name, sizeof name, "random %lu x %lu, %lu-bit, draw %lu", m, n, bits, draw);
        reduce_and_check(t, name, in, m, n);
    } else {
        t->skipped += 1;
    }
    basis_free(in, m);
}

static void knapsack_basis(tally *const t, unsigned long const m, unsigned long const bits, unsigned long const draw) {
    clll_vector *const in = basis_new(m, m + 1);
    for (unsigned long i = 0; i < m; i++) {
        mpz_set_ui(in[i].v[i], 1);
        rng_entry(in[i].v[m], bits);
        mpz_abs(in[i].v[m], in[i].v[m]);
    }
    char name[64];
    snprintf(name, sizeof name, "knapsack %lu x %lu, %lu-bit, draw %lu", m, m + 1, bits, draw);
    reduce_and_check(t, name, in, m, m + 1);
    basis_free(in, m);
}

static void identity_basis(tally *const t, unsigned long const n) {
    clll_vector *const in = basis_new(n, n);
    for (unsigned long i = 0; i < n; i++) {
        mpz_set_ui(in[i].v[i], 1);
    }
    char name[64];
    snprintf(name, sizeof name, "identity %lu x %lu", n, n);
    reduce_and_check(t, name, in, n, n);
    basis_free(in, n);
}

static void large_multiplier_basis(tally *const t, mpz_t const c) {
    clll_vector *const in = basis_new(3, 3);
    mpz_set_ui(in[0].v[0], 1);
    mpz_set(in[1].v[0], c);
    mpz_set_ui(in[1].v[1], 1);
    mpz_set(in[2].v[0], c);
    mpz_set(in[2].v[1], c);
    mpz_set_ui(in[2].v[2], 1);
    char name[128];
    gmp_snprintf(name, sizeof name, "multiplier %Zd", c);
    reduce_and_check(t, name, in, 3, 3);
    basis_free(in, 3);
}

static void mixed_basis(tally *const t, unsigned long const m, unsigned long const rounds, unsigned long const draw) {
    clll_vector *const in = basis_new(m, m);
    for (unsigned long i = 0; i < m; i++) {
        for (unsigned long l = 0; l < m; l++) {
            rng_entry(in[i].v[l], 8);
        }
    }
    if (independent(in, m)) {
        mpz_t x;
        mpz_init(x);
        for (unsigned long r = 0; r < rounds; r++) {
            unsigned long const i = rng_next() % m;
            unsigned long const j = (i + 1 + rng_next() % (m - 1)) % m;
            rng_entry(x, 64 + rng_next() % 37);
            for (unsigned long l = 0; l < m; l++) {
                mpz_addmul(in[i].v[l], x, in[j].v[l]);
            }
        }
        mpz_clear(x);
        char name[64];
        snprintf(name, sizeof name, "mixed %lu x %lu, %lu rounds, draw %lu", m, m, rounds, draw);
        reduce_and_check(t, name, in, m, m);
    } else {
        t->skipped += 1;
    }
    basis_free(in, m);
}

static bool report(char const *const group, tally const *const t) {
    if (t->ok) {
        printf("OK: %s: %lu reductions satisfy every condition", group, t->reductions);
    } else {
        printf("MISMATCH: %s: some of %lu reductions fail a condition", group, t->reductions);
    }
    if (t->skipped > 0) {
        printf(" (dependent draws skipped: %lu)", t->skipped);
    }
    printf("\n");
    return t->ok;
}

int main(void) {
    bool ok = true;

    tally square = {0, 0, true};
    for (unsigned long m = 2; m <= 8; m++) {
        for (unsigned long draw = 0; draw < 2; draw++) {
            random_basis(&square, m, m, 4, draw);
            random_basis(&square, m, m, 20, draw);
            random_basis(&square, m, m, 60, draw);
        }
    }
    ok = report("square bases, m = n = 2..8, 4/20/60-bit entries", &square) && ok;

    static unsigned long const shapes[][2] = {{2, 6}, {3, 5}, {4, 7}};
    tally wide = {0, 0, true};
    for (size_t s = 0; s < sizeof shapes / sizeof shapes[0]; s++) {
        for (unsigned long draw = 0; draw < 3; draw++) {
            random_basis(&wide, shapes[s][0], shapes[s][1], 4, draw);
            random_basis(&wide, shapes[s][0], shapes[s][1], 20, draw);
            random_basis(&wide, shapes[s][0], shapes[s][1], 60, draw);
        }
    }
    ok = report("non-square bases, m x n = 2 x 6, 3 x 5, 4 x 7", &wide) && ok;

    tally knapsack = {0, 0, true};
    for (unsigned long m = 4; m <= 12; m += 4) {
        for (unsigned long draw = 0; draw < 2; draw++) {
            knapsack_basis(&knapsack, m, 30, draw);
        }
    }
    ok = report("knapsack bases (e_i | a_i), m = 4, 8, 12, 30-bit a_i", &knapsack) && ok;

    tally large = {0, 0, true};
    {
        mpz_t c;
        mpz_init(c);
        unsigned long const exponents[] = {32, 64};
        for (size_t e = 0; e < sizeof exponents / sizeof exponents[0]; e++) {
            for (long sign = 1; sign >= -1; sign -= 2) {
                mpz_ui_pow_ui(c, 2, exponents[e]);
                mpz_sub_ui(c, c, 1);
                mpz_mul_si(c, c, sign);
                large_multiplier_basis(&large, c);
                mpz_ui_pow_ui(c, 2, exponents[e]);
                mpz_mul_si(c, c, sign);
                large_multiplier_basis(&large, c);
            }
        }
        for (long sign = 1; sign >= -1; sign -= 2) {
            mpz_ui_pow_ui(c, 2, 100);
            mpz_add_ui(c, c, 3);
            mpz_mul_si(c, c, sign);
            large_multiplier_basis(&large, c);
        }
        mpz_clear(c);
    }
    for (unsigned long m = 3; m <= 5; m += 2) {
        for (unsigned long draw = 0; draw < 3; draw++) {
            mixed_basis(&large, m, 2 * m, draw);
        }
    }
    ok = report("large multipliers: 3 x 3 with c up to 2^100, mixed 3 x 3 and 5 x 5", &large) && ok;

    tally small = {0, 0, true};
    identity_basis(&small, 1);
    identity_basis(&small, 3);
    identity_basis(&small, 6);
    for (unsigned long n = 1; n <= 5; n += 2) {
        random_basis(&small, 1, n, 20, 0);
    }
    {
        clll_vector *const none = basis_new(0, 3);
        mpq_t alpha;
        mpq_init(alpha);
        mpq_set_ui(alpha, 3, 2);
        clll_lll(none, 0, alpha);
        mpq_clear(alpha);
        basis_free(none, 0);
        small.reductions += 1;
    }
    ok = report("identity bases, m = 1 and m = 0", &small) && ok;

    return ok ? 0 : 1;
}
