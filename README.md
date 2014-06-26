# About

This is the `andi` program for estimating genome diversity.

# Build Instructions

For the latest [stable release](https://github.com/EvolBioInf/andi/releases) of `andi` download the tar ball. If you'd like to contribute to this software, feel free to create a fork of our [git repository](https://github.com/EvolBioInf/andi) and send pull requests.

## Prerequisites

This program has the following dependencies: [kseq.h](http://lh3lh3.users.sourceforge.net/kseq.shtml), RMQ_succinct and [libdivsufsort](https://code.google.com/p/libdivsufsort/). The former two are included in the `lib` directory. Please make sure you installed the latter before attempting a build. If you did get the source, not as a tarball, but straight from the git repository, you will also need autoconf, automake and the libtool.

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

The release of this software is accompanied by a paper from Haubold et al. (in prep.). It explains the used *anchor distance* strategy in great detail. The `maf2phy.awk` script used in the validation process is located under `scripts`. Simulations were done using our own [simK](http://guanine.evolbio.mpg.de/bioBox/) tool.

## Data Sets

1. 29. E. coli strains: [data](http://guanine.evolbio.mpg.de/andi/eco29.fasta.gz)
2. 109 E. coli ST131 strains ([paper](http://www.pnas.org/content/early/2014/03/28/1322678111.abstract)): 
	* [99 newly sequenced strains](https://github.com/BeatsonLab-MicrobialGenomics/ST131_99)
	* [10 previously published strains](http://guanine.evolbio.mpg.de/andi/st131_extra.tgz)
3. 3085 Streptococcus pneumoniae strains ([paper](http://www.nature.com/ng/journal/v46/n3/full/ng.2895.html)): ftp://ftp.sanger.ac.uk/pub/pathogens/Streptococcus/pneumoniae/Maela_assemblies.tgz

# License

Copyright © 2014 Fabian Klötzl  
License GPLv3+: GNU GPL version 3 or later.

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at <http://gnu.org/licenses/gpl.html>.

## Known Bugs

There is a slight memory leak in the RMQ lib used. But it is relatively small (some KB) compared to the overall memory consumption (a few GB).
