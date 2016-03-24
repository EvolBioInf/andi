[![Build Status](https://travis-ci.org/EvolBioInf/andi.svg?branch=master)](https://travis-ci.org/EvolBioInf/andi) [![Coverage Status](https://coveralls.io/repos/EvolBioInf/andi/badge.svg?branch=master)](https://coveralls.io/r/EvolBioInf/andi?branch=master)

# About

This is the `andi` program for estimating the evolutionary distance between closely related genomes. These distances can be used to rapidly infer phylogenies for big sets of genomes. Because `andi` does not compute full alignments, it is so efficient that it scales even up to thousands of bacterial genomes.

This readme covers all necessary instructions for the impatient to get `andi` up and running. For extensive instructions please consult the [manual](andi-manual.pdf).


# Installation and Usage

Stable versions of `andi` are available via package managers. For manual installation see below.

For Debian and Ubuntu (starting 16.04):

    sudo apt-get install andi

For OS X with Homebrew:

    brew install science/andi

For ArchLinux with aura:

    sudo aura -A andi
    
With a successful installation you can get the usage instructions via `--help` or the man page.

    $ andi --help
    $ man andi

You can simply use `andi` with your genomes in `FASTA` format.

    $ andi S1.fasta S2.fasta
    2
    S1     0.0  0.1
    s2     0.1  0.0

From this distance matrix the phylogeny can be inferred via neighbor-joining. Check the [manual](andi-manual.pdf) for a more thorough description.


## Manual installation

If your system does not support one of the above package managers you have to manually build the latest [stable release](https://github.com/EvolBioInf/andi/releases) from a tarball. See the [manual](andi-manual.pdf) for extensive building instructions.

This program has the following external dependencies: [libdivsufsort](https://github.com/y-256/libdivsufsort) and the [GSL](https://www.gnu.org/software/gsl/). Please make sure you installed both before attempting a build. If you did get the source, not as a tarball, but straight from the git repository, you will also need the autotools.

Assuming you have installed all prerequisites, building is as easy as follows.

    $ autoreconf -fi -Im4  # optional when build from tarball
    $ ./configure
    $ make
    $ make install

Excessive build instructions are located in `INSTALL`. 

# Links and Additional Resources

The release of this software is accompanied by a paper from [Haubold et al.](http://bioinformatics.oxfordjournals.org/content/31/8/1169). It explains the used *anchor distance* strategy in great detail. The `maf2phy.awk` script used in the validation process is located under `scripts`. Simulations were done using our own [simK](http://guanine.evolbio.mpg.de/bioBox/) tool. For a demo visualising the internals of andi visit our [GitHub pages](http://evolbioinf.github.io/andi/).

## Data Sets

1. 29 E. coli and Shigella strains: [data](http://guanine.evolbio.mpg.de/andi/eco29.fasta.gz)
2. 109 E. coli ST131 strains ([paper](http://www.pnas.org/content/early/2014/03/28/1322678111.abstract)): 
	* [99 newly sequenced strains](https://github.com/BeatsonLab-MicrobialGenomics/ST131_99)
	* [10 previously published strains](http://guanine.evolbio.mpg.de/andi/st131_extra.tgz)
3. 3085 Streptococcus pneumoniae strains ([paper](http://www.nature.com/ng/journal/v46/n3/full/ng.2895.html)): ftp://ftp.sanger.ac.uk/pub/pathogens/Streptococcus/pneumoniae/Maela_assemblies.tgz

## License

Copyright © 2014 - 2016 Fabian Klötzl  
License GPLv3+: GNU GPL version 3 or later.

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at <http://gnu.org/licenses/gpl.html>.

Some files may be licensed differently.

## Contact

In case of bugs or unexpected errors don't hesitate to send me a mail: kloetzl@evolbio.mpg.de
