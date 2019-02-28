# Reproducing Ebola infection analysis with MSF

Open the link 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84188

and download GSE84188_norm_matrix_corrected.txt.gz file.

Extract the columns
 
`Donor9.EBOV 6H, 1D, 2D with Donor9.Mock 6H, 1D, 2D`

`Donor10.EBOV 6H, 1D, 2D with Donor10.Mock 6H, 1D, 2D`

`Donor11.EBOV 6H, 1D, 2D with Donor11.Mock 6H, 1D, 2D`

these are the 3 replicates for each time-point. Using these replicates Differential gene expression analysis was carried out with EdgeR (version 3.4.2) using an upper-quartile normalization.

# Modulated Sub-graph Finder

Modulated Sub-graph Finder (**MSF**) is used to find the significantly dis-regulated sub-graphs or cluster of genes from the host cell signaling network, giving these sub-graphs an overall significance of modulation by combining the individual p-values of the genes derived from differential genes expression analysis. 

## Prerequisites

* Java version 8
* Jdk 1.8

## Download Linux && Windows

Download the ModulatedSubgraphFinder jar file.


## Input Data Preparation

To find the Modulated sub-graphs you need two files. One file is tab-separated text file containing the output of EdgeR analysis. Second file is a directed adjacency list with 3 columns, first two columns with interacting genes and the third column direction of the interaction. The gene identifiers should be same in both the input files. How the interaction file was created is given below.

**MSF** has seven argument parameters 

* `-p`	The path to differential gene expression analysis file 
* `-i`	The path to network file (Interaction file)
* `-t`	Software used for differential  gene expression analysis (DEseq2 or EdgeR)
* `-e`	The extension limit (1 to 3 genes extension)
* `-m`	The merging limit (1 to 3 genes merging)
* `-k`	Output extra files (InitialGraphs, ExtendedGraphs and MergedGraphs)
* `-o`	The path to output folder

Navigate to folder containing the jar file, EdgeR result and Reactome interaction file and run command

`java -jar ModulatedSubPathFinder.jar -p EdgeR_Results.txt -i Interactions.txt -m 2 -e 2 -t EdgeR -k yes -o /home/Documents/`

The default extension limit is 2 and merging limit is 1.

## Reactome Interaction File

The interaction file used to build and test **MSF** was downloaded from reactome https://reactome.org/download-data under Functional interactions (FIs) version 2016. This file was pre-filtered before use. All direct activation and inhibition interactions were filtered to be used by **MSF** testing.

### Filtering

First interactions that had no directions were removed from the file. Then indirect interactions were removed followed by interactions like expression regulation and catalyze. The remaining direct interactions were checked manually.






## Built With

Eclipse Neon

## Version

Version 1.0
DOI: 10.5281/zenodo.1400242 

## Authors

**Mariam Farman** 

## License

This project is licensed under the MIT license.





