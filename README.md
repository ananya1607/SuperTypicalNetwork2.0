# SuperTypicalNetwork2.0

Construct and characterise enhancer-promoter interaction networks from chromatin loop data, distinguishing SE constituents from typical enhancers.

## What the pipeline does
Given a set of chromatin loops and coordinates for super enhancer constituents, typical enhancers and promoters, the pipeline:
1. Reads and harmonises the inputs, optionally relaxing the loop anchors by a fixed number of base pairs to improve overlap detection
2. Intersects the loop anchors with SEc, TE and promoter coordinates
3. Classifies each anchor as P, SE or TE. In the case where an anchor overlaps more than one element type, it applies a fixed priority: *P > SEc > TE*
4. Builds a node and edge list. chooses all those loops that resolve to a relevant element (P, SEc, TE) to generate edge list and node list, making sure to remove self edges, and duplicated edges are collapsed, and the maximum edge weight is chosen.
5. Constructs an undirected graph
6. Computes per-node centrality within each component (degree, strength, betweenness, closeness, eigenvector, harmonic)
7. Saves the outputs  


## Packages Required:
- data.table
- dplyr
- tidyr
- purrr
- tibble
- igraph
- GenomicRanges
- IRanges
- InteractionSet

## Requirements: 
- Loops file containing loop anchors, and interaction frequency (normalised)
- Super enhancer constituents' coordinates
- Typical enhancer coordinates
- Promoter coordinates

## Inputs and Input formats
### For creation of enhancer promoter interaction network
#### 1. config_file.csv
A .csv file containing the following information
- cell_line = column giving information about the sample
- condition =  column giving information about whether the sample is baseline or differentiated/treated
- loops_file = absolute path to the loops file containing the information about loops called
- se_file = absolute path to super enhancer constituent file
- te_file = absolute path to Typical enhancer file
- prom_file = absolute path to promoter file
- out_dir = absolute path to the output subdirectory where the results should be pasted
- relax_by =  numeric value in bp, using which the loop anchors will be relaxed (can be set to 0 if relaxation is not wanted)


sample config file:
```
"cell_line","condition","loops_file","se_file","te_file","prom_file","out_dir","relax_by"
"foreskin_keratinocyte","day0","~/Documents/data/foreskin_keratinocyte/loops_raw_day0.csv","~/Documents/data/foreskin_keratinocyte/se_coords_foreskin_keratinocyte.csv","~/Documents/data/foreskin_keratinocyte/te_coords_foreskin_keratinocyte.csv","~/Documents/data/foreskin_keratinocyte/pr_coords_foreskin_keratinocyte.csv","~/Documents/results/foreskin_keratinocyte/day0",0 

```

#### 2. loops_file.csv 
A .csv file containing the following information
- chr1,start1,end1: chromosome, start and end position of the first anchor of the loop
- chr2,start2,end2: chromosome, start and end position of the second anchor of the loop
- value: normalised interaction frequency of the loop (given by the loop caller)
- loop_length: abs(start1-end2), determining the length of the loop extruded

sample:
```
chr1,start1,end1,chr2,start2,end2,value,loop_length
11,44285000,44290000,11,44610000,44615000,2.707059,325000
11,19935000,19940000,11,20615000,20620000,3.5252247,680000
11,31335000,31340000,11,32180000,32185000,4.7458797,845000

```

#### 3. se_file.csv, te_file.csv, prom_file.csv
A .csv file containing the following information
- chrom, start, end: chromosome, start and end position of the genomic element
- id: identifier of the genomic element

sample:
```
chrom,start,end,id
1,109383139,109384593,SE_01_003900656_1
1,109385296,109386726,SE_01_003900656_2
1,109387828,109393871,SE_01_003900656_3
1,109393909,109395142,SE_01_003900656_4
1,110337431,110342543,SE_01_003900218_1
```

### For calculation of partner conservation index
#### 1. output.csv
Output.csv file after running 2_execution.R script

sample:

```
"name","LA1_type","Cluster","degree","strength","betweenness","closeness","eigenvector","harmonic","cluster_type","cell_line","condition"
"SE_03_000100405_1","SE",1,2,16.719615,0,0.379222385438733,0.261716130055924,39.0539323779481,"P-SE-TE","foreskin_keratinocyte","day0"
"TE_03_000100755","TE",1,10,103.015875,25,0.476838058858646,1,64.2753427094221,"P-SE-TE","foreskin_keratinocyte","day0"
"TE_03_000101088","TE",1,5,18.32796,13,0.338459268546109,0.113019058990157,23.1568205917296,"P-SE-TE","foreskin_keratinocyte","day0"
"TE_03_000101830","TE",1,5,67.10218,0,0.425964389354371,0.942143939718321,56.7961570097734,"P-SE-TE","foreskin_keratinocyte","day0"
```

#### 2. edges.csv
The edges.csv file is a file output after running the 2_execution.R script in the condition subdirectory

```
LA1,LA2,value,interaction_type
ENSR00000449872,ENSR00001318017,2.7851398,P-P
ENSR00000449872,TE_01_003909714,4.5675254,P-TE
ENSR00000458419,ENSR00001319928,3.3427083,P-P
ENSR00000706021,TE_01_003909018,5.180407,P-TE
ENSR00000862678,TE_01_003901800,3.7546432,P-TE
```



## Main Functions
##### <ins>Function name: process_loops_se_te_pr()</ins>

Input:
- loops_file -> path to the loops file (as csv)
- prom_file -> path to the promoter coordinates file (as csv)
- se_file -> path to the super enhancer constituent coordinates file (as csv)
- te_file -> path to the typical enhancer coordinates file (as csv)
- relax_by -> number of bases to relax the loop coordinates by, for better intersection (in bp)

Output: 
- se, te and pr <- dataframe of super enhancers constituents, typical enhancers and promoters respectively, having consistent column names and formatted in order
- loops <- dataframe of loops

##### <ins>Function name: give_truth_table()</ins>

Input:
loops, te, se, pr <- dataframe of loops, typical enhancers, super enhancer constituents and promoters, respectively, created using process_loops_se_te_pr() function 

Output:
- loops_gi <- GInteraction object of loop coordinates
- promoter_gr, super_enhancers_gr, typical_enhancers_gr <- GRanges object of P, SEc and TE coordinates respectively
- hits1_p, hits2_p <- S4Vector object containing the hits between P and anchor1 and anchor2 of loops respectively
- hits1_se, hits2_se <- S4Vector object containing the hits between SEc and anchor1 and anchor2 of loops respectively
- hits1_te, hits2_te <- S4Vector object containing the hits between TE and anchor1 and anchor2 of loops respectively

##### <ins>Function name: giveloop_ids()</ins>

Input:
- loops_gi <- GInteraction object of loop coordinates
- pr_gr, te_gr, se_gr <- granges object of super enhancer constituents, typical enhancers and promoters respectively, having consistent column names and formatted in order
- hits1_p, hits2_p, hits1_se, hits2_se, hits1_te, hits2_te <- S4Vector object containing the hits between the genomic elements (P, SEc, TE) and the loop anchors

Output
loops_gi <- GInteraction object containing the IDs of the genomic elements that have hits against the loop anchors

##### <ins>Function name: classify_interaction(), make_loop_interaction_categories()</ins>

Input:
- loops_gi <- GInteraction object containing the ids of intersecting genomic elements of each loop anchor in metadata

Output:
- loops_gi <- GInteraction object containing the interaction type of each loop observed
*Also prints the number of loops in each interaction category as a bar plot*
*priority chosen (P>SEc>TE)*

##### <ins>Function name: give_total_edges_nodes()</ins>

Input:
- loops_g1 <- GRanges object of Loops containing the interaction type
- loops <- loops data frame 

What it does: for each of the loops that are observed, it removes redundancy in interaction type, 
               makes it consistent and fetches the relevant IDs based on the interaction type and 
               generates pairwise combinations between the nodes to make the edge dataframe and the nodes dataframe
               generates some graphs

Output:
- total_nodes <- df containing all the nodes and the type (SEc, TE or P)
- total_edges <- df containing all the edges (IDs), their interaction frequency and interaction type

##### <ins>Function name: calculate_pci()</ins>
Input:
- cell_line <- sample name as mentioned in the output csv file
- csv_path <- absolute path to network properties .csv file
- ed1, ed2, ed3 <- absolute paths to sample_edges.csv pertaining to each conditon, ed1 corresponds to the baseline condition, ed3 corresponds to most differentiated condition
- condition_levels <- a list containing condition levels from least differentiated to most differentiated condition

What it does: for a given sample with three conditions, and network properties of nodes from all conditions, this function calculates the partner conservation index for the Super enhancer constituents and Typical enhancers present in each condition. 

## Execution
Run the 2_execution.R code in RStudio, making sure that the inputs are formatted properly and that absolute paths are given.

*In case the user wants to build a bipartite network, then set is_bipartite = TRUE in 2_execution.R code*

## Outputs
This code returns results across conditions, as well as combined results
1. all_enhancer_properties.csv
   .csv file containing network properties, sample identity and cluster information of the nodes across all the given input samples
2. all_results.rds
   RDS file containing the analysed file that includes additional information about the loop object used, nodes and edges, cluster classification, granges loop object (based on which the enhancer promoter network is created), network igraph object, enhancer properties, cell line and condition (just as an easier way to access information)
3.  pipeline.log
   ile containing the execution of code (if the pipeline doesn't run properly, refer to this file to look at errors)
4. outputs in condition subdirectories:
   1. igraph_object.rds
      RDS file containing the network in igraph format
   2. nodes.csv, edges.csv
      csv files containing the edges and nodes for a particular condition
   3. enhancer_properties.csv
      csv file containing the network properties for nodes belongng to a particular condition 


