# Duplication mutation analysis

This workflow analyses the evolution of the  duplication in the RSV G gene.
It is based on the without-G workflow:[RSV without G](https://github.com/LauraU123/without_G_workflow) workflow, which constructs RSV phylogenetic analyses which exclude G from the treebuilding steps.

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



