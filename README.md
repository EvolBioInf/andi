# About

This is the `andi` program for estimating genome diversity.

# License

This program is released under the GNU Public License 3.

# Build Instructions

## Prerequisites

This program has the following dependencies: kseq.h, RMQ_succinct, libdivsufsort and gsl. The former two are included in the `lib` directory. Please make sure, you have the other two installed before attempting a build. If you did get the source not as a tarball but straight from the git repository you will also need the autoconf tools.

## Compiling

Assuming you have installed all prerequisites, building is as easy as follows.

	$ ./configure
	$ make
	$ make install

If the build was successful you can get the usage instructions via `--help`.

	$ andi --help

Code documentation is provided via doxygen.

	$ make docs

## Known Bugs

There is a slight memory leak in the RMQ lib used. But it is relatively small (some KB) compared to the overall memory consumption (a few GB).

