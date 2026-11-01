ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(top_srcdir)/include
AM_CFLAGS = $(CLLL_STD_CFLAGS) $(CLLL_WARN_CFLAGS) $(CLLL_PLT_CFLAGS)

AM_DISTCHECK_CONFIGURE_FLAGS = CFLAGS='-g -O2 -Werror'

# ---------------------------------------------------------------------------
# libclll
# ---------------------------------------------------------------------------

include_HEADERS = include/clll.h

clllincludedir = $(includedir)/clll

clllinclude_HEADERS = \
	include/clll/gram_schmidt_orthogonalization.h \
	include/clll/lll.h \
	include/clll/vector.h

lib_LTLIBRARIES = libclll.la

libclll_la_SOURCES = \
	src/assert_and_assume.h \
	src/gram_schmidt_orthogonalization.c \
	src/gram_schmidt_orthogonalization_row.h \
	src/lll.c \
	src/vector.c

libclll_la_LIBADD = $(GMP_LIBS)
libclll_la_LDFLAGS = -version-info $(CLLL_VERSION_INFO) -no-undefined

# ---------------------------------------------------------------------------
# pkg-config and CMake
# ---------------------------------------------------------------------------

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = clll.pc

cmakedir = $(libdir)/cmake/clll
cmake_DATA = cmake/clllConfig.cmake cmake/clllConfigVersion.cmake

# Filled in by make, not configure: CMake would read ${prefix} as its own variable.
cmake/clllConfig.cmake: $(srcdir)/cmake/clllConfig.cmake.in Makefile
	$(AM_V_GEN)$(MKDIR_P) cmake && $(SED) \
	  -e 's|@includedir[@]|$(includedir)|g' \
	  -e 's|@libdir[@]|$(libdir)|g' \
	  $(srcdir)/cmake/clllConfig.cmake.in > $@

CLEANFILES = cmake/clllConfig.cmake
EXTRA_DIST = Doxyfile README.md cmake/clllConfig.cmake.in

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

check_PROGRAMS = \
	tests/test_example \
	tests/test_gram_schmidt \
	tests/test_lll_reduced \
	tests/test_lll_short_vector

TESTS = $(check_PROGRAMS)

tests_test_example_SOURCES = tests/test_example.c tests/lll_conditions.h tests/test_util.h
tests_test_example_LDADD = libclll.la $(GMP_LIBS)

tests_test_gram_schmidt_SOURCES = tests/test_gram_schmidt.c tests/test_util.h
tests_test_gram_schmidt_LDADD = libclll.la $(GMP_LIBS)

tests_test_lll_reduced_SOURCES = tests/test_lll_reduced.c tests/lll_conditions.h tests/test_util.h
tests_test_lll_reduced_LDADD = libclll.la $(GMP_LIBS)

tests_test_lll_short_vector_SOURCES = tests/test_lll_short_vector.c tests/lll_conditions.h tests/test_util.h
tests_test_lll_short_vector_LDADD = libclll.la $(GMP_LIBS)
