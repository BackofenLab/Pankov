Pankov: Probabilistic Sankoff-like simultaneous alignment and folding of RNAs inspired by Markov chains

This is the Pankov repository for structural alignment of RNAs, which is an extended branch of LocARNA package.


## Installation 
This branch is planned to be merged into the main LocARNA repository, on a long term. Please refer to README.txt and INSTALL files to follow the detailed installation instruction.

### Quick guide
1. Make sure autotools dependencies are installed and run: `libtoolize && aclocal &&  autoheader&& automake -a -c && autoconf`

2. Make sure Vienna RNA tools are installed and run the make configuration: `./configure --with-vrna=PATH_TO_VIENNA RNA`

3. Then run `make`

## Usage
After installation you can run Pankov binary file `pankov` from the bin directory in the installation path or the src directory.
Please always pass these parameters to use Pankov which uses the conditional probably `--track-closing-bp  --use-cond`, the parameters are mandatory for the proper run of the program.

Call example:
`./src/pankov seqA.fa seqB.fa --write-structure --track-closing-bp  --use-cond`

Further info can be found using `--help` option as well as the LocARNA package documentation.
