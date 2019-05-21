 Pankov: Probabilistic Sankoff-like RNA alignment and folding using Markovproperty with a near-exact energy model

This is the Pankov repository for structural alignment of RNAs, which is an extended branch of LocARNA package.


## Installation 
This branch is planned to be merged into the main LocARNA repository, on a long term. Please refer to README.txt and INSTALL files and follow the detailed installation instruction.

After installation you can run Pankov binary file `sparse_n` from the bin directory in the installation path.
Please always pass these parameters to use Pankov which uses the conditional probably `--track-closing-bp  --use-cond`, the parameters are mandatory for the proper run of the program.

Call example:
`./bin/sparse_n seqA.fa seqB.fa --write-structure --track-closing-bp  --use-cond`
