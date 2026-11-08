# clll

Exact integer LLL lattice basis reduction in C99, built on
[GMP](https://gmplib.org/).

`clll` implements Algorithm 2 (Gram-Schmidt orthogonalization over the integers)
and Algorithm 5 (the LLL algorithm, verified integer version) from R. Thiemann
et al., *Formalizing the LLL Basis Reduction Algorithm and the LLL
Factorization Algorithm in Isabelle/HOL*.
All arithmetic is exact: there is no floating point, and every intermediate
value is an integer.

## Requirements

- A C99 compiler (GCC, Clang, ...) and `make`
- GMP and its development headers:

  | Distribution    | Package        |
  |-----------------|----------------|
  | Debian, Ubuntu  | `libgmp-dev`   |
  | Fedora, RHEL    | `gmp-devel`    |
  | Arch            | `gmp`          |
  | Gentoo          | `dev-libs/gmp` |

Building from a git checkout also needs autoconf >= 2.71, automake >= 1.16
and libtool.
A release tarball does not.

## Building from a release tarball

`VERSION` is the release's version number.

```sh
tar xzf clll-VERSION.tar.gz
cd clll-VERSION
./configure
make
make check
sudo make install
```

`make install` puts everything under `/usr/local` by default.
If programs then cannot find `libclll.so`, run `sudo ldconfig` once.
`sudo make uninstall` removes everything again.

## Building from a git checkout

The repository holds only the hand-written files, so `configure` has to be
generated first:

```sh
autoreconf -fi
```

After that the steps are the same as for a tarball.
Building in a separate directory keeps the source tree clean:

```sh
mkdir build
cd build
../configure
make
make check
```

## Configure options

`./configure --help` lists every option.

When the compiler accepts it, `configure` adds `-fno-plt`, so that the library
calls GMP without going through a PLT stub.
This measured 0.3-1% faster.
`CFLAGS='-g -O2 -fplt'` turns it off.
`-fplt` on its own would also drop the default `-O2`.

## Using the library

```c
#include <clll.h>
#include <stdio.h>

int main(void) {
    /* A basis of a 2-dimensional lattice. */
    clll_vector basis[2];
    clll_vector_init(&basis[0], 2);
    clll_vector_init(&basis[1], 2);
    mpz_set_si(basis[0].v[0], 201);
    mpz_set_si(basis[0].v[1], 37);
    mpz_set_si(basis[1].v[0], 1648);
    mpz_set_si(basis[1].v[1], 297);

    /* Reduce it in place with alpha = 3/2. */
    mpq_t alpha;
    mpq_init(alpha);
    mpq_set_ui(alpha, 3, 2);
    clll_lll(basis, 2, alpha);
    mpq_clear(alpha);

    for (int i = 0; i < 2; i++) {
        gmp_printf("(%Zd, %Zd)\n", basis[i].v[0], basis[i].v[1]);
        clll_vector_clear(&basis[i]);
    }
    return 0;
}
```

This prints a basis of the same lattice that is reduced with respect to
α = 3/2, for example:

```
(1, 32)
(40, 1)
```

A lattice usually has many reduced bases.
`clll_lll` guarantees only the conditions documented in `include/clll/lll.h`,
not these particular vectors.

Vectors follow GMP's conventions:
initialize each one with `clll_vector_init` before use, clear it with
`clll_vector_clear` afterwards, and read or write entries directly as
`basis[i].v[j]`.
Like every GMP rational, `alpha` must be in canonical form.
`mpq_set_ui(alpha, 3, 2)` already is, but a fraction with common factors,
such as 6/4, needs `mpq_canonicalize` before the call.
C++ code can include `<clll.h>` directly.

### Documentation

The headers in `include/clll/` document every function.
To build the HTML version, run `doxygen` in the source directory and open
`doc/html/index.html`.

### With pkg-config

```sh
cc app.c $(pkg-config --cflags --libs clll) -o app
```

If `clll` is installed somewhere pkg-config does not look, point it there with
`PKG_CONFIG_PATH=DIR/lib/pkgconfig`.
Run the program with `LD_LIBRARY_PATH=DIR/lib` too, unless `DIR` is a standard
location.

### With CMake

```cmake
find_package(clll 1.0 REQUIRED)
target_link_libraries(app PRIVATE clll::clll)
```

For a non-standard prefix, configure your project with
`-DCMAKE_PREFIX_PATH=DIR`.

## Making a release

```sh
make distcheck
```

This builds (with `-Werror`), tests, installs and uninstalls a fresh copy of
the project, then writes `clll-VERSION.tar.gz`.
That tarball contains the generated `configure`, so the people who download it
need no autotools.
`VERSION` is the version number set in `configure.ac`.
