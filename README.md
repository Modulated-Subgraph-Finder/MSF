# Modulated Sub-graph Finder

Modulated Sub-graph Finder (**MSF**) used to find the significantly dis-regulated sub-graphs or cluster of genes from the host cell signaling network, giving these sub-graphs an overall significance of modulation by combining the individual p-values of the genes derived from differential genes expression analysis. 

## Prerequisites

* Java version 8
* Jdk 1.8

## Download Linux & Windows

Download the ModulatedSubgraphFinder.jar file and the example files for tutorial.


## Input Data Preparation

To find the Modulated sub-graphs you need two files. One file is tab-separated text file containing the output of DESeq2/EdgeR analysis. Second file is a directed adjacency list with 3 columns, first two columns with interacting genes and the third column direction of the interaction. The gene identifiers should be same in both the input files. Example files included in Docs folder.

## Tutorial

**MSF** has seven argument parameters 

* `-p`	The path to differential gene expression analysis file 
* `-i`	The path to network file (Interaction file)
* `-t`	Software used for differential  gene expression analysis (DEseq2 or EdgeR)
* `-e`	The extension limit (1 to 3 genes extension)
* `-m`	The merging limit (1 to 3 genes merging)
* `-k`	Output extra files (InitialGraphs, ExtendedGraphs and MergedGraphs)
* `-o`	The path to output folder

Navigate to folder containing the jar file and the example files and run command

`java -jar ModulatedSubgraphFinder.jar -p ExampleDEGAnalysis_DEseq2.csv -i ExampleInteractions.csv -t DEseq2 -k yes -o /home/Documents/`

The default extension and merging limit is 2. If you do not want to output extra files skip the -k argument.

### Output Files

#### InitialGraphs

This is a text file that contains the initial sub-graphs that are found by combining the individual P-values of the gene. The snapshot shows there core sub-graphs were found, written in a line with the combined p-value of the sub-graph at the end.

`[ppp2r1a, tgfbr2, tgfb2, ppp2ca] 1.1103580381200124E-18`

`[smad2, skp1a, smad3] 1.6569054494446366E-7`
 
`[nog, bmp2, bmp5, bmpr2, gdf7, bmp6] 3.482452728268697E-4`

#### ExtendedGraphs

This file shows if any sub-graphs were extended by adding genes beyond its immediate neighborhood. For example graph is extended by 4 genes.

`[ppp2r1a, tgfbr2, tgfb2, ppp2ca] 1.1103580381200124E-18`

`[smad2, skp1a, smad3, smad7, ifng, acvr1, smad5] 1.9580322709452766E-10`
 
`[nog, bmp2, bmp5, bmpr2, gdf7, bmp6] 3.482452728268697E-4`

#### MergedGraphs

This is the text file showing if any or all the (extended/unextended) sub-graphs merge with each other or not. In this case the graphs did not merge.


`[ppp2r1a, tgfbr2, tgfb2, ppp2ca] 1.1103580381200124E-18`

`[smad2, skp1a, smad3, smad7, ifng, acvr1, smad5] 1.9580322709452766E-10`
 
`[nog, bmp2, bmp5, bmpr2, gdf7, bmp6] 3.482452728268697E-4`

#### SourcesAndSinks

This output file gives details about the genes found in the sub-graphs. It shows graph number followed by genes in the graph. Then each gene from the graph, its fold change, individual p-value and in the last if it was identified as a source, intermediate or sink. The sources have impact score against it.

`[Graph 1]`
 
`[ppp2r1a, tgfbr2, tgfb2, ppp2ca]`

`[ppp2r1a, 343.680151982419, 1.4E-7, Source, 3/4 = 75.000]`

`[tgfbr2, 154.225593438929, 0.03294361100258, Intermediate]`

`[tgfb2, 505.093869916889, 4.0E-5, Sink]`

`[ppp2ca, 525.820339004501, 6.44660834994E-4, Source, 3/4 = 75.000]`

`[Graph 2]`

`[smad2, skp1a, smad3, smad7, ifng, acvr1, smad5]`

`[smad2, 28.1552749566574, 1.51528460481E-4, Intermediate]`

`[skp1a, 843.231746773892, 0.014525956708056, Sink]`

`[smad3, 176.442049495414, 0.146036813182214, Source, 2/7 = 28.571]`

`[smad7, 606.619477359416, 0.562318600418486, Source, 4/7 = 57.143]`

`[ifng, 135.622073938663, 2.12566137193E-4, Sink]`

`[acvr1, 100.463207564239, 0.017364084706018, Source, 2/7 = 28.571]`

`[smad5, 675.586240911088, 0.005281642821289, Sink]`

`[Graph 3]`

`[nog, bmp2, bmp5, bmpr2, gdf7, bmp6]`

`[nog, 48.0012337399973, 0.005498079382534, Intermediate]`

`[bmp2, 186.03559350895, 0.237098851922717, Sink]`

`[bmp5, 57.9279489451268, 0.419057181261883, Sink]`

`[bmpr2, 814.343493775222, 0.017218841862483, Source, 2/6 = 33.333]`

`[gdf7, 142.559614133974, 0.011845850200211, Source, 5/6 = 83.333]`

`[bmp6, 442.598787102063, 0.403941986198972, Sink]`


#### NetworkFile

A text file with directed adjacency list for **MSF** identified modulated sub-graphs. This file could further be used to visualize the sub-graphs in other tools for example in Cytoscape. The last two columns are used as edge attributes to be imported in Cytoscape.

`tgfbr2 ppp2r1a |- 0  1 `

`tgfb2 tgfbr2 <- 0  1 `

`tgfbr2 ppp2ca |- 0  1 `

`skp1a smad2 <- 0  1 `

`skp1a smad3 <- 0  1 `

`nog bmp2 -> 1  0 `

`nog bmp5 -> 1  0 `

`bmp2 bmpr2 <- 0  1 `

`nog gdf7 <- 0  1 `

`nog bmp6 -> 1  0 `

`smad7 smad2 -> 1  0 `

`ifng smad7 <- 0  1 `

`acvr1 smad3 <-> 1  1 `

`acvr1 smad5 -> 1  0 `

#### SourceWeight

This file has node attributes for each gene to be imported into Cytoscape. First column with gene name, second column the node size depending on the source weight and the last column is LogFoldchange to show type of regulation.

`ppp2r1a  750.000  -1.15250657211695` 
 
`tgfbr2  10.000  0.357347052992014`

`tgfb2  10.000  -0.273710254909285`

`ppp2ca  750.000  -0.225269184801659`

`smad2  10.000  1.20822323638527`

`skp1a  10.000  0.47396118683819`

`smad3  285.714  -0.193102843027487`

`smad7  571.429  0.100940948274761`

`ifng  10.000  -0.411195572576727`

`acvr1  285.714  -0.360351622623896`

`smad5  10.000  0.43197314616271`

`nog  10.000  -0.397877024648883`

`bmp2  10.000  -0.135915332971465`

`bmp5  10.000  -0.21104227566159`

`bmpr2  333.333  0.255224975603306`

`gdf7  833.333  0.334306828205764`

`bmp6  10.000  0.079033613236002`


## Tutorial MSF to StringApp
### Prerequisites

* Cytoscpae
* StringApp (Cytoscape Plugin)
* enhanceGraphics (Cytoscape Plugin)

### Getting Started

#### Importing MSF networks

`File-> Import-> Network-> File-> NetworkFile.text-> Advanced Options-> Tick Delimiter **SPACE**-> Untick Use first line as column names`

Click column 1 and select as source node
Click on column 2 and set as target node 
Click on column 3 and set as interaction type.
Click ok.

The MSF network has been imported to Cytoscape, to add directional follow the next steps

Click on Style
Select Edge attributes
Drop down Source Arrow Shape

`Column -> Column 6`

 `Mapping Type -> Discrete Mapping`
 
 `For 1 select the arrow ->`
 
Drop down Target Arrow Shape

`Column -> Column 4`

 `Mapping Type -> Discrete Mapping`
 
 `For 1 select the arrow ->`
 
 To add node attributes follow steps
 
`File-> Import-> Table-> File-> SourceWeight.text 
Where to import Table Data -> To selected netowkrs only
Network List -> NetworkFile.text
Import Data as -> Node Table Columns
Key Column for Networks -> Shared name
Advanced Options-> Tick Delimiter **SPACE**-> Untick Use first line as column names`

Click on Style
Select Node attributes
Drop down size 

`Column -> Column 3`

 `Mapping Type -> Passthrough Mapping`
 
 Drop down Border Paint 
 
 `Column -> Column 5`
 
` Mapping Type -> Continous Mapping`

 `Set Border width to 25`

#### Importing StringApp network
After successful download of Cytoscape and StringApp and enhanceGraphics inside cytoscape. Follow the steps

`File->Import->Network->Public Databases`

A new window would open then select `Data source ` as `String : protein query`. Next select the organism as Homo sapiens. Then enter the complete list of gene in MSF sub-graphs and click `import`. The string network from String database would be downloaded with the user input gene list but interactions from String database. Since we would have our own interaction we can remove the String interactions by

`select->Edges->Select all edges` The right click on any edge in the network and click `Edit->cut`. All the string database edges would be removed.

#### Merging networks

We would overlay our network on the String network with our edges only. To do so 

`Tools->Merge->Networks` a new window would open, using the `Union` option select the string network and **MSF** identified network. From the advance options  select the appropriate  names of genes to be merged For example `SharedName` for **MSF** network and `display name` for String network. Then click `merge`.

Now you would see **MSF** identified network overlapped  with string network.

#### Enrichment Analysis

To perform gene enrichment analysis on **MSF** network got to `Apps->STRING->Set as STRING network`. Next go to `Apps->STRING Enrichment->Retrieve Functional Enrichment`. This would perform functional enrichment for GO terms and KEGG pathways. Below the network a table with enrichment analysis would open. You can filter the table to a particular type of enrichment for example KEGG pathway using the filter symbol on top left of table. Using the settings option on top right you could make single genes to pie charts to show for example in which KEGG pathways is the gene present with different colors for different KEGG pathways.

#### Exporting Network

You can either save the network as a png or a session file.

## Reactome Interaction File

The interaction file used to build and test **MSF** was downloaded from reactome https://reactome.org/download-data under Functional interactions (FIs) version 2016. This file was pre-filtered before use. All direct activation and inhibition interactions were filtered to be used by **MSF** testing.

### Filtering

First interactions that had no directions were removed from the file. Then indirect interactions were removed followed by interactions like expression regulation and catalyze. The remaining interaction were checked manually.






## Built With

Eclipse Neon

## Version

Version 1.0 of the tool.

## Authors

**Mariam Farman** 
farman@tbi.univie.ac.at

## License

This project is licensed under the Creative Commons Attribution 4.0 International License.




