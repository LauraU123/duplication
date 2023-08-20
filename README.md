# Duplication mutation analysis

This workflow analyses the evolution of the  duplication in the RSV G gene, which is 72nt in length in RSV-A and 60nt in RSV-B. 
It is based on the without-G workflow:[RSV without G](https://github.com/LauraU123/without_G_workflow) workflow, which constructs RSV phylogenetic analyses which exclude G from the treebuilding steps.



First, a tree based on the without_G workflow is constructed. Following this step, 
the sequences of all branches are reconstructed based on the tree and the root sequence. 
These branches are then aligned in a series of steps, and the duplicated regions in the G gene are then cut out.

Mutations in each copy of the duplication (preduplication and postduplication copies 1 and 2) are  found by comparing the mutations in each branch from tree tip to root. 


Finally, the mutations are scaled by total tree length in which they occur. This is graphed separately for synonymous and nonsynonymous mutations.
Two sample KS statstics and mutation rates based on the MLE of the Poisson distribution are calculated for each copy of the duplication. 

 
The duplicated region in the G gene is extracted from the phylogenetic tree and aligned in a series of steps, and subsequently mutations in each duplication copy (preduplication, postduplication copy 1 and copy 2) are analysed. 

### Input Data

Input metadata and sequences for RSV-A and RSV-B are available on data.nextstrain.org

* [RSV-A sequences](https://data.nextstrain.org/files/workflows/rsv/a/sequences.fasta.xz)
* [RSV-A metadata](https://data.nextstrain.org/files/workflows/rsv/a/metadata.tsv.gz)

* [RSV-B sequences](https://data.nextstrain.org/files/workflows/rsv/b/sequences.fasta.xz)
* [RSV-B metadata](https://data.nextstrain.org/files/workflows/rsv/a/metadata.tsv.gz)


These data are generously shared by labs around the world and deposited in NCBI genbank by the authors.
Please contact these labs first if you plan to publish using these data.


## Installation

Follow the standard [installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.


## Output
 
 The workflow generates the following outputs:

	- graph of each mutation copy location - synonymous (png)
	- graph of each mutation copy location - nonsynonymous (png)
	- cumulative distribution for each mutation copy (png)
	- KS goodness of fit test for each distribution (csv)
	- mutation rate of each copy (csv)

## Running the workflow

To run the workflow, Snakmake is required, as well as Nextstrain
Once installed, run:


``` snakemake --cores all ```  



