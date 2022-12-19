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
 - hdf5r
 - compositions

USAGE
-----

Usage informations will be displayed by invoking:

	NetworkNullHPC.sh -h

The expected input OTU table format is a TAB or space separated file containing samples as rows and OTUs as columns, with the first column containing one field less, so that the first row and the first column could be used as the OTU and sample names, respectively. If the OTU table was edited under Windows OS, pay attention to use Unix compliant end of line (\n).

Three outputs are produced:
 - a text summary of the option used and the output network sizes
 - a co-occurrence edge list (if significant edges found)
 - a co-exclusion edge list (if significant edges found)
  The edge lists contain three columns: the two OTUs linked by the edge and the median Spearman's rank correlation value calculated over the bootstraped random noise addition matrices.

REFERENCES
----------

Connor, N., Barberán, A., & Clauset, A. (2017). Using null models to infer microbial co-occurrence networks. PLOS ONE, 12(5), e0176751. doi:10.1371/journal.pone.0176751

Lentendu, G., & Dunthorn, M. (2021). Phylogenetic relatedness drives protist assembly in marine and terrestrial environments. Global Ecology and Biogeography, 30(7), 1532–1544. doi: https://doi.org/10.1111/geb.13317
