# Build instructions

Right now what you need is more or less:

    $ mkdir build
    $ cd build
    $ cmake ../src
    $ make

It will hopefully build if everything else is ok (see below). The
finished binaries will be found in `build/bin/`.

The build directory is self-contained. To clean up a build, you can
delete the entire build directory or its contents. You can make
several independent builds with different settings (see below for
options) simply by doing them in separate directories.

## Requirements

You need at least the following libraries installed:

- FFTW, built with Fortran bindings
- NetCDF, built with Fortran bindings
- XRAYLIB

You need at least GFortran 4.3, but for some reason GFortran 4.5
does not work. This is some kind of GFortran issue as far as I know.

All the libraries should be built with the same compilers which you 
use to compile the programs. Using a different GFortran caused some 
problems at least with NetCDF.

You also need CMake if it's not already installed on your system.

### Build options

There are some options you can give CMake.

Setting `-DBUILD_SHARED_LIBS=ON` or `-DWITH_OPENMP=ON` will cause
the build to fail on my system. I haven't been able to find out why
yet.

If you have an AMD processor and the optimized AMD Core Math Library
(ACML), you can use the ACML random number generators by setting
`-DWITH_ACML=ON`. This should work.

To enable FITS output, you can use `-DWITH_CFITSIO=ON`, which
requires the `cfitsio` library.


## Installing on OS X

On OS X, you should install the latest version of GFortran (as I
write this it's version 4.8). Get it from:
http://gcc.gnu.org/wiki/GFortranBinaries

I recommend you install Homebrew, a package manager for installing
all sorts of UNIX tools on a Mac.

See instructions for installing Homebrew at:
http://mxcl.github.com/homebrew/

With Homebrew, installing the required libraries is easy:

    $ brew install cmake cfitsio
    $ brew install fftw --with-fortran
    $ brew install netcdf --enable-fortran

I will make a Homebrew formula for xraylib as well. When it exists,
you can also do:

    $ brew install xraylib

But before that, you need to install xraylib yourself. It should not
be harder than downloading the latest source tarball from 
https://github.com/tschoonj/xraylib/downloads
and then installing it the old-fashioned way (see its INSTALL file
for details).
