# Global_Proteome_Analysis

Manual for Python script for the analysis of large-scale temporal proteomics data

Data input

Changes in protein abundances between different time points using the TMT quantification option provided by the PEAKS Q software tool. Each TMT-experiment was analyzed separately. Auto normalization was used for quantitative analysis of the proteins, in which the global ratio was calculated from the total intensity of all labels in all quantifiable peptides. Quantitative analysis was performed using protein identifications containing at least 2 unique peptides, which peptide identifications were filtered against 1% FDR. However, the number of unique peptides for a positive protein quantification could be altered here. The significance method for evaluating the observed abundance changes was set to ANOVA and the significance score was expressed as the -10xlog10(p), where p is the significance testing p-value. The p-value represents the likelihood that the observed change is caused by random chance. Results from PEAKS Q were exported to “proteins.csv”, containing the quantified proteins for each analysis separately. This file was used as data input for the proteomics data analysis pipeline in Python. These quantification files need to be placed in the folder “crude_data” before the Python analysis. Furthermore, to access the simplified Central Carbon Metabolism KEGG database, an Excel file “Simplified Enzymes Central Carbon Metabolic Pathway.xlsx” should be located in the folder “KEGG_GO_Uniprot”.

Script

This script consists of several consecutive sections, of which the content and functions will be explained below. This script is dependent on the functions in the “all_proteomics_functions.py” file. To this end, the “all_proteomics_functions.py” is run first to import all functions. 


0.	Import the required module: biopython

1.	**Parameters.** Fill in the required parameters, tailored to your research. Parameters include the TMT labelling scheme (label_scheme), the technical replicates, (TMT-labelled) groups, different experiments, biological replicates. Furthermore, the different growth phases and the time in hours (relative to glucose depletion) at which the samples were taken are required (xdata and xaxis parameters). Finally, output directories need to be specified.
2.	**Generate gene names (and optionally GO terms) for each yeast genome accession number.** For interpretation of the data it is required to convert the accession numbers to gene names, and for functional interpretation, a GO term analysis can be performed.
3.	**Import and modify the proteomics datafile.** The data from all TMT experiments are imported into a Python dictionary using the “import_prot_data” function, and thereafter the data are normalized towards the bridging samples using the “adjust_ratios” function. Proteins that have a number of identified unique peptides that is below the selected threshold (2 by default, can be adjusted in section 1), are removed. Subsequently, all TMT data are concatenated into one large dataframe using the “concatenate_data” function. Lastly, the data are divided into the different biological experiments and replicates using the “separate_data” function. This concatenation and separation is required because one TMT experiment could contain the data from different biological experiments.
4.	**KEGG pathways. **Import a selected set of proteins from the central carbon metabolism (CCM) from an Excel file “Simplified Enzymes Central Carbon Metabolic Pathway.xlsx” and define the selected pathways by annotating these in the parameter “pathways”. Furthermore, import all listed KEGG pathways for the selected organism (S. cerevisiae by default), using the “KEGG_tool” function. 
5.	**Make clustermap to check reproducibility of technical and biological replicates.** This is done using the “make_clustermap” function, modified from the Seaborn “clustermap” function, in which Euclidean distances metric and the average linkage method were used. Here, only proteins detected in all three biological are visualized.  
6.	**Average biological replicates.** After the reproducibility of the experiments is analysed, calculate the average intensities and ratios of the biological replicates, to use for all further analyses. To this end, the ratio/fold change of each protein in a specific condition was calculated relative to the bridging sample (bridging sample = 1). Furthermore, calculate the standard deviations and note the number of replicates that a protein was identified in for subsequent statistical evaluation. 
7.	**FC plot marker proteins.** For data validation, a set of marker proteins are defined which are analysed using fold change (FC) plots. These could be altered in the parameter “genelist”. 
8.	**Make clustermap of CCM of averaged biological replicates.** The central carbon metabolism pathways and proteins are plotted on a heatmap to roughly examine overall temporal trends of the different biological experiments, using the “make_clustermap” function. Proteins that are not detected in any of the experiments are not visualized, solely printed. No significance threshold is 
9.	**Calculate p-values between selected experiments. **In order to assess the statistical significance of observed differences between two desired biological experiments, p-values are calculated using the means, standard deviations, and number of observations. These parameters are saved in the data frame of each experiment. This analysis can either be done to compare two biological experiments (here: MG vs WT or AN vs O2), or to compare different timepoints within one biological experiment (here: transition into stationary). If you want to compare specific growth phases, this could be defined in ‘xaxis’ (in AN vs O2).
10.	**Make fold change (FC) plots including significance values for single proteins of interest.** For example, all proteins in a certain pathway can be selected. For a good comparison of the separate experiments (e.g. comparison of conditions or strains), the fold changes remained normalised to the bridging samples. The p-values of significant differences in means between two experiments were annotated with ***, **, *, or ˆ for p-values below 0.001, 0.01, 0.05, and 0.1, respectively. In the last section of the Python script (section “A. FC plots additional proteins of interest”), examples can be found of different ways to retrieve and plot proteins of interest trough the 
11.	**STRING functional enrichment analysis and volcano plot comparing desired experiments or timepoints.** To determine whether specific GO terms or KEGG pathways are enriched within the significantly up- or down- regulated proteins under a particular condition compared to another, a functional enrichment analysis was done using the STRING database via the “STRING” function. For example, the growth phases that you would like to compare, could be defined in the parameter “xaxis”. In addition, the experiment could also be selected in the parameter “exp”. The “STRING” function requires a dictionary containing the up- and downregulated proteins (with "upregulated" and "downregulated" as keys, obtained via the function “up_down_regulated”, the species identifier (4932 for S. cerevisiae), the functional categories that should be assessed ("PMID", "Component", "Process", "Function", "KEGG", and/or "Interpro"), the FDR threshold (<0.05 in this study), and an optional set of "background genes" with as alternative background the whole species proteome. A function named “backgroundgene_2_string” was written that retrieves the protein-specific string identifiers for the background genes/proteins, which are all the proteins that were detected in the experiments that were compared. The STRING function gives a printed tabular output containing for each significantly enriched or depleted term/pathway: the functional category, the term-id/pathway-id, the FDR, the p-value, the number of up-/downregulated proteins which can be annotated to that term, the total number of proteins in the background set that can be annotated to that term, and the term description. This table is also exported as a CSV file. 
Volcano plots are used to visualize global proteome differences between two experiments. These plots were generated using the “gene_exp.volcano” function; a modified version of the “GeneExpression.volcano” function from the "Bioinfokit.visuz" module in Python. The log2 fold change threshold can be altered, which is defined by “lfc_thr”. Furthermore, proteins that are quantified in less than 2 biological replicates are not visualized. 
