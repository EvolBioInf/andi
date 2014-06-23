# About

This is the `andi` program for estimating genome diversity.

# Build Instructions

For the latest [stable release](https://github.com/EvolBioInf/andi/releases) of `andi` download the tar ball. If you'd like to contribute to this software, feel free to create a fork of our [git repository](https://github.com/EvolBioInf/andi) and send pull requests.

## Prerequisites

This program has the following dependencies: [kseq.h](http://lh3lh3.users.sourceforge.net/kseq.shtml), RMQ_succinct, [libdivsufsort](https://code.google.com/p/libdivsufsort/) and the [GNU Scientific Library](http://www.gnu.org/software/gsl/). The former two are included in the `lib` directory. Please make sure, you have the other two installed before attempting a build. If you did get the source, not as a tarball, but straight from the git repository, you will also need autoconf, automake and the libtool.

## Compiling

Assuming you have installed all prerequisites, building is as easy as follows.

	$ ./configure
	$ make
	$ make install

Excessive build instructions are located in `INSTALL`. If the build was successful you can get the usage instructions via `--help` or the man page.

	$ andi --help
	$ man andi

Code documentation is provided via doxygen.

	$ make docs

# Links and Additional Resources

The release of this software is accompanied by a paper from Haubold et al. (in prep.). It explains the used *anchor distance* strategy in great detail. As soon as the various data sets used for validation are available online we will link to them here. The `maf2phy.awk` script used in the validation process is provided in the `scripts` directory.

# License

This program is released under the GNU Public License 3.

## Known Bugs

There is a slight memory leak in the RMQ lib used. But it is relatively small (some KB) compared to the overall memory consumption (a few GB).

