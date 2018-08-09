# Modulated Sub-graph Finder

Modulated Sub-graph Finder (**MSF**) used to find the significantly dis-regulated sub-graphs or cluster of genes from the host cell signaling network, giving these sub-graphs an overall significance of modulation by combining the individual p-values of the genes derived from differential genes expression analysis. 

## Prerequisites

* Java version 8
* Jdk 1.8

## Download Linux && Windows

Download the ModulatedSubgraphFinder jar file and the example files for tutorial.


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

`java -jar ModulatedSubPathFinder.jar -p ExampleDEAnalysis.txt -i ExampleInteractions.txt -t DEseq2 -k -o /home/Documents/`

The default extension and merging limit is 2.

### Output Files

#### InitialGraphs

This is a text file that contains the initial sub-graphs that are found by combining the individual P-values of the gene. The snapshot shows there core sub-graphs were found, written in a line with the combined p-value of the sub-graph at the end.

`[ppp2r1a, tgfbr2, tgfb2, ppp2ca] 1.1103580381200124E-18`

`[smad2, skp1a, smad3] 1.6569054494446366E-7`
 
`[nog, bmp2, bmp5, bmpr2, gdf7, bmp6] 3.482452728268697E-4`

#### ExtendedGraphs

This file shows if any sub-graphs were extended by adding genes beyond its immediate neighborhood.

`[ppp2r1a, tgfbr2, tgfb2, ppp2ca] 1.1103580381200124E-18`

`[smad2, skp1a, smad3, smad7, ifng, acvr1, smad5] 1.9580322709452766E-10`
 
`[nog, bmp2, bmp5, bmpr2, gdf7, bmp6] 3.482452728268697E-4`

#### MergedGraphs

This is the text file showing if any or all the (extended/unextended) sub-graphs merge with each other or not.

#### SourcesAndSinks

This output file gives details about the genes found in the sub-graphs. It shows graph number followed by genes in the graph. Then each gene from the graph, its interactions in the graph, its fold change, individual p-value and in the last if it was identified as a source, intermediate or sink.

`[Graph 1]`

`[ppp2r1a, tgfbr2, tgfb2, ppp2ca]`

`[ppp2r1a, -tgfbr2, -1.15250657211695, 1.4E-7, Sink]`

`[tgfbr2, -ppp2ca-ppp2r1a-tgfb2, 0.357347052992014, 0.03294361100258, Intermediate]`

`[tgfb2, -tgfbr2, -0.273710254909285, 4.0E-5, Source]`

`[ppp2ca, -tgfbr2, -0.225269184801659, 6.44660834994E-4, Sink]`

`[Graph 2]`

`[smad2, skp1a, smad3, smad7, ifng, acvr1, smad5]`

`[smad2, -smad7-acvr1-skp1a, 1.20822323638527, 1.51528460481E-4, Sink]

 [skp1a, -smad2-smad3, 0.47396118683819, 0.014525956708056, Source]`
 
`[smad3, -smad7-skp1a, -0.193102843027487, 0.146036813182214, Sink]

 [smad7, -smad2-smad3-smad5-ifng, 0.100940948274761, 0.562318600418486, Intermediate]`
 
`[ifng, -smad7, -0.411195572576727, 2.12566137193E-4, Source]`

`[acvr1, -smad5-smad2, -0.360351622623896, 0.017364084706018, Source]`

`[smad5, -smad7-acvr1, 0.43197314616271, 0.005281642821289, Sink]`

`[Graph 3]`

`[nog, bmp2, bmp5, bmpr2, gdf7, bmp6]`

`[nog, -bmp2-bmp5-bmp6-gdf7, -0.397877024648883, 0.005498079382534, Source]`

`[bmp2, -bmpr2-nog, -0.135915332971465, 0.237098851922717, Intermediate]`

`[bmp5, -bmpr2-nog, -0.21104227566159, 0.419057181261883, Intermediate]`

`[bmpr2, -bmp2-bmp5-bmp6-gdf7, 0.255224975603306, 0.017218841862483, Sink]`

`[gdf7, -bmpr2-nog, 0.334306828205764, 0.011845850200211, Intermediate]`

`[bmp6, -bmpr2-nog, 0.079033613236002, 0.403941986198972, Intermediate]`


#### NetworkFile

A text file with directed adjacency list for **MSF** identified modulated sub-graphs. This file could further be used to visualize the sub-graphs in other tools for example in Cytoscape.

`ppp2r1a tgfbr2 |-`

`tgfbr2 tgfb2 <-`

`ppp2ca tgfbr2 |-`

`smad2 skp1a <-`

`smad3 skp1a <-`

`nog bmp2 ->`

`nog bmp5 ->`

`bmpr2 bmp2 <-`

`gdf7 nog <-`

`nog bmp6 ->`

`smad7 smad2 ->`

`smad7 ifng <-`

`acvr1 smad3 <->`

`acvr1 smad5 ->`

## Tutorial MSF to StringApp
### Prerequisites

* Cytoscpae
* StringApp (Cytoscape Plugin)
* enhanceGraphics (Cytoscape Plugin)

### Getting Started
#### Importing networks

After successful  download of Cytoscape and StringApp and enhanceGraphics inside cytoscape. Follow the steps

`File->Import->Network->Public Databases`

A new window would open then select `Data source ` as `String : protein query`. Next select the organism as Homo sapiens. Then enter the gene list and click `import`. The string network from String database would be downloaded with the user input gene list but interactions from String database. Since we would have our own interaction we can remove the String interactions by

`select->Edges->Select all edges` The right click on any edge in the network and click `Edit->cut`. All the string database edges would be removed. To import our network (**MSF** identified NetworkFile)


`File->Import->Network->File` and select the networkFile.

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

## License

This project is licensed under the Creative Commons Attribution 4.0 International License.




