#include <clll/vector.h>

#include "assert_and_assume.h"

#include <gmp.h>

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

void clll_vector_init(clll_vector *const x, unsigned long const len) {
    x->len = len;
    x->v = NULL;
    if (len == 0) {
        return;
    }
    if (len > SIZE_MAX / sizeof(mpz_t)) {
        abort();
    }

    void *(*allocate)(size_t);
    mp_get_memory_functions(&allocate, NULL, NULL);
    x->v = allocate(len * sizeof(mpz_t));
    for (unsigned long i = 0; i < len; i++) {
        mpz_init(x->v[i]);
    }
}

void clll_vector_clear(clll_vector *const x) {
    if (x->len == 0) {
        return;
    }
    for (unsigned long i = 0; i < x->len; i++) {
        mpz_clear(x->v[i]);
    }

    void (*deallocate)(void *, size_t);
    mp_get_memory_functions(NULL, NULL, &deallocate);
    deallocate(x->v, x->len * sizeof(mpz_t));
}

void clll_vector_set(clll_vector *const rop, clll_vector const *const op) {
    ASSERT_AND_ASSUME(rop->len == op->len);
    for (unsigned long i = 0; i < rop->len; i++) {
        mpz_set(rop->v[i], op->v[i]);
    }
}

void clll_vector_swap(clll_vector *const x, clll_vector *const y) {
    clll_vector const t = *x;
    *x = *y;
    *y = t;
}

void clll_vector_inner_product(mpz_t rop, clll_vector const *const u, clll_vector const *const v) {
    ASSERT_AND_ASSUME(u->len == v->len);
    // Leave rop alone if it is already 0, so its buffer is kept.
    if (mpz_sgn(rop) != 0) {
        mpz_set_ui(rop, 0);
    }
    for (unsigned long i = 0; i < u->len; i++) {
        mpz_addmul(rop, u->v[i], v->v[i]);
    }
}

void clll_vector_add(clll_vector *const rop, clll_vector const *const u, clll_vector const *const v) {
    ASSERT_AND_ASSUME(rop->len == u->len && u->len == v->len);
    for (unsigned long i = 0; i < rop->len; i++) {
        mpz_add(rop->v[i], u->v[i], v->v[i]);
    }
}

void clll_vector_sub(clll_vector *const rop, clll_vector const *const u, clll_vector const *const v) {
    ASSERT_AND_ASSUME(rop->len == u->len && u->len == v->len);
    for (unsigned long i = 0; i < rop->len; i++) {
        mpz_sub(rop->v[i], u->v[i], v->v[i]);
    }
}

void clll_vector_mul(clll_vector *const rop, clll_vector const *const u, clll_vector const *const v) {
    ASSERT_AND_ASSUME(rop->len == u->len && u->len == v->len);
    for (unsigned long i = 0; i < rop->len; i++) {
        mpz_mul(rop->v[i], u->v[i], v->v[i]);
    }
}
