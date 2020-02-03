
# Cointerfaces
**Cointerfaces** is a protocol to build a joint alignment (or paired alignment) of a pair of interacting domains and detect inter-domain physical contacts using the joint alignment built. It takes as input two HMM profiles and return two FASTA files codifying the joint alignment as well as a file with coevolutioary scores for each pair of positions. This protocol was used in the publication (["Conservation of coevolving protein interfaces bridges prokaryote-eukaryote homologies in the twilight zone"](http://www.pnas.org/content/113/52/15018.full)).

Cointerfaces uses the following software:
- [fmpl](https://github.com/simomarsili/fmpl) implementation of plmDCA developed by Simone Marsili (@simomarsili) to compute the DCA coevolutionary model.
- [HMMER](http://hmmer.org/) sequence analysis software developed in Eddy lab.
- [HHsuite](https://github.com/soedinglab/hh-suite) sequence analysis software developed in Soeding lab.
- [cdhit](http://weizhongli-lab.org/cd-hit/) developed by Dr. Weizhong Li at Dr. Adam Godzik's Lab.

## Installation

The following instructions are given for its installation in Ubuntu. The installation in any linux distribution should be very similar, with the use of a different package manager.

Cointerfaces run on linux requires perl, python and fortran compiler:
> sudo apt-get install perl python3 gfortran

The following python packages are needed
> pip install bio numpy

And for perl:
> cpan -i Getopt::Long FindBin Cwd File::Basename List::Util

Extract the data
> cd data/ensembl_bacteria/
> tar xvzf selection_100_genomes.tar.gz
> cd ../..

A selection of 100 genomes is included for testing. For real case scenarios, please download the whole collection (ensembl bacteria release 23) following the instructions in file data_download_instructions.

**Modify the project_root variable in src/localpahts file to match the path of the project directory**

Compile fmpl program
> cd src/protocol/dca_model/mpl/
> make

To recompile, use ``make clean; make``


## Test

The included test should take around 5-10 minutes
> cd test
> ./run_test.py

The test will quantify if there are relevant differences between the computed contact prediction against a stored one. The files *old in test/ correspond to the stored ones which can be used for comparison (e.g. the compare the joint alignment)


## Miscelaneous

By default, cointerfaces uses two strategies for pairing interacting sequences:
1. Genomic proximity: Homologous sequences to the HMM profiles have been found in the same gene or a pair genes close in the their genome (by default, 300 base pairs)
2. Uniquines in genomes: Only one hit is found for each of the two HMM profiles in a given genomes. This option can be set with a parameter, it should be turn off if there are no external evidence that the two HMM profiles correspond to a known pair of interacting protein domain families.

The protocol can be ran from a pair of protein sequences. This can be achieved by making HMM profiles using jackhmmer, but it's needed to use the same version of HMMER (e.g. use jakchmmer in third_party_software/hmmer/binaries)

It is possible to use a newer version of the data or software (modifying localpaths accordingly). But this will probably require some minor modifications in the protocol to work properly.




