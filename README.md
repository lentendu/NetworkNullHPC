NetworkNullHPC re-implements in R and bash the OTU co-occurrence network inferrence method developped by Connor, Barberàn and Clauset (2017) and use the SLURM job scheduler to parallelize the computations.
The method is generalized to negative correlations in order to assess both co-occurrence and co-exclusion.

INSTALLATION
------------

To install, just clone the repository

	git clone http://github.com/lentendu/NetworkNullHPC.git/

and add the newly created directory into your PATH variable.

NetworkNullHPC is only available for Linux server using the [SLURM job scheduler](https://slurm.schedmd.com/).

It request [GNU Parallel](https://www.gnu.org/software/parallel/) and R loaded in the environment as a module.

In R, the following packages are requested:
 - plyr
 - tidyr
 - dplyr
 - vegan
 - Hmisc
 - igraph
 - foreach
 - doParallel
 - rhdf5

USAGE
-----

Usage informations will be displayed by invoking:

	NetworkNullHPC.sh -h

The expected input OTU table format contains sample as rows and OTU as columns, with the first column containing one field less to be properly read by the R read.table command.

Three outputs are produced:
 - a co-occurrence edge list
 - a co-exclusion edge list
 - a text summary of the option used

REFERENCE
---------

Connor, N., Barberán, A., & Clauset, A. (2017). Using null models to infer microbial co-occurrence networks. PLOS ONE, 12(5), e0176751. doi:10.1371/journal.pone.0176751
