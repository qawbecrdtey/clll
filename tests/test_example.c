/**
 * @file
 * @brief Reduces the basis of Example 1 from the paper with @f$\alpha = 3/2@f$
 *        and checks that the result is a reduced basis of the same lattice.
 *
 * A lattice usually has many reduced bases, so the result is not compared with
 * Example 3 of the paper. Any basis that meets the conditions in lll.h passes.
 */

#include <clll.h>

#include "lll_conditions.h"

#include <gmp.h>

#include <stdbool.h>
#include <stdio.h>

int main(void) {
    clll_vector *const in = basis_new(3, 3);
    mpz_set_ui(in[0].v[0], 1);
    mpz_set_ui(in[0].v[1], 1894885908);
    mpz_set_ui(in[1].v[1], 1);
    mpz_set_ui(in[1].v[2], 1894885908);
    mpz_set_str(in[2].v[2], "2147483648", 10);

    mpq_t alpha;
    mpq_init(alpha);
    mpq_set_ui(alpha, 3, 2);
    clll_vector *const out = basis_copy(in, 3, 3);
    clll_lll(out, 3, alpha);

    printf("Reduced basis (alpha = 3/2):\n");
    for (unsigned long r = 0; r < 3; r++) {
        gmp_printf("  (%Zd, %Zd, %Zd)\n", out[r].v[0], out[r].v[1], out[r].v[2]);
    }

    char const *const label = "Example 1, alpha = 3/2";
    bool const ok = check_reduced(label, in, out, 3, 3, alpha);
    if (ok) {
        printf("OK: %s: the result is a reduced basis of the same lattice\n", label);
    }

    mpq_clear(alpha);
    basis_free(out, 3);
    basis_free(in, 3);
    return ok ? 0 : 1;
}
