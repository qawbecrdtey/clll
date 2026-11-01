/**
 * @file
 * @brief Assertions that become optimizer hints when NDEBUG is defined.
 */

#ifndef CLLL_ASSERT_AND_ASSUME_H
#define CLLL_ASSERT_AND_ASSUME_H

/**
 * @def CLLL_CONST_FUNCTION
 * @brief Marks a function as depending only on its arguments, with no side
 *        effects.
 *
 * Clang ignores an assumption that contains a call it cannot prove to be free
 * of side effects (-Wassume), even a call to a visible static inline function.
 * Mark a function with this so that calls to it can appear in
 * ASSERT_AND_ASSUME().
 */
#if defined(__GNUC__) || defined(__clang__)
#define CLLL_CONST_FUNCTION __attribute__((const))
#else
#define CLLL_CONST_FUNCTION
#endif

/**
 * @def ASSERT_AND_ASSUME(expr)
 * @brief Asserts @p expr, or lets the optimizer assume it under NDEBUG.
 *
 * Without NDEBUG this is assert(expr). With NDEBUG, @p expr is passed to the
 * optimizer as an assumption and is not evaluated, so it costs nothing at run
 * time.
 *
 * @note Clang is tested first because it also defines __GNUC__, and clang-cl
 *       defines _MSC_VER. It uses __builtin_assume, since Clang before 19 parses
 *       __attribute__((assume(...))) as an OpenMP attribute. GCC before 13 has
 *       no assumption that leaves its argument unevaluated, so there @p expr is
 *       evaluated, and the optimizer removes it only when it can see that it
 *       has no side effects.
 */
#ifdef NDEBUG

#if defined(__clang__)
#define ASSERT_AND_ASSUME(expr) \
    do { \
        __builtin_assume(expr); \
    } while (0)
#elif defined(_MSC_VER)
#define ASSERT_AND_ASSUME(expr) \
    do { \
        __assume(expr); \
    } while (0)
#elif defined(__GNUC__) && __GNUC__ >= 13
#define ASSERT_AND_ASSUME(expr) \
    do { \
        __attribute__((assume(expr))); \
    } while (0)
#elif defined(__GNUC__)
#define ASSERT_AND_ASSUME(expr) \
    do { \
        if (!(expr)) { \
            __builtin_unreachable(); \
        } \
    } while (0)
#else
#define ASSERT_AND_ASSUME(expr) \
    do { \
    } while (0)
#endif

#else

#include <assert.h>

#define ASSERT_AND_ASSUME(expr) \
    do { \
        assert(expr); \
    } while (0)
#endif

#endif
