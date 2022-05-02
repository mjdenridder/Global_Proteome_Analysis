
"""
Script for processing large scale proteomics data. 

Author: Wiebeke van den Brandeler

December 2021
"""

#%% 0. Import required modules 

import all_proteomics_functions as apf
import re
import numpy as np
import pandas as pd
import os
from   matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import math
import scipy.stats as st
import matplotlib.patches as mpatches


#%% 1. data specific parameters and directory specifications

print("\n\033[4m1. Loading parameters\033[0m")

Unique_thr   = 2 #threshold for the minimum amount of found unique peptides
num_brep     = 3 #amount of biological replicates

label_scheme = {'TMT01': {'Group 1': 'MG_AN_16_ME',
                          'Group 2': 'MG_AN_16_LE',
                          'Group 3': 'MG_AN_16_ES',
                          'Group 4': 'MG_AN_16_MS',
                          'Group 5': 'bridge_TMT01'},
                'TMT02': {'Group 1': 'MG_AN_17_ME',
                          'Group 2': 'MG_AN_17_LE',
                          'Group 3': 'MG_AN_17_ES',
                          'Group 4': 'MG_AN_17_MS',
                          'Group 5': 'bridge_TMT02'},
                'TMT03': {'Group 1': 'MG_AN_18_ME',
                          'Group 2': 'MG_AN_18_LE',
                          'Group 3': 'MG_AN_18_ES',
                          'Group 4': 'MG_AN_18_MS',
                          'Group 5': 'bridge_TMT03'},
                'TMT04': {'Group 1': 'MG_O2_16_ME',
                          'Group 2': 'MG_O2_16_ED',
                          'Group 3': 'MG_O2_16_MD',
                          'Group 4': 'MG_O2_16_MS',
                          'Group 5': 'bridge_TMT04'},
                'TMT05': {'Group 1': 'MG_O2_17_ME',
                          'Group 2': 'MG_O2_17_ED',
                          'Group 3': 'MG_O2_17_MD',
                          'Group 4': 'MG_O2_17_MS',
                          'Group 5': 'bridge_TMT05'},
                'TMT06': {'Group 1': 'MG_O2_18_ME',
                          'Group 2': 'MG_O2_18_ED',
                          'Group 3': 'MG_O2_18_MD',
                          'Group 4': 'MG_O2_18_MS',
                          'Group 5': 'bridge_TMT06'},
                'TMT07': {'Group 1': 'MG_O2_16_LE',
                          'Group 2': 'MG_O2_17_LE',
                          'Group 3': 'MG_O2_18_LE',
                          'Group 4': 'bridge_TMT07'},
                'TMT08': {'Group 1': 'WT_AN_19_ME',
                          'Group 2': 'WT_AN_19_LE',
                          'Group 3': 'WT_AN_19_ES',
                          'Group 4': 'WT_AN_19_MS',
                          'Group 5': 'bridge_TMT08'},
                'TMT09': {'Group 1': 'WT_AN_20_ME',
                          'Group 2': 'WT_AN_20_LE',
                          'Group 3': 'WT_AN_20_ES',
                          'Group 4': 'WT_AN_20_MS',
                          'Group 5': 'bridge_TMT09'},
                'TMT10': {'Group 1': 'WT_AN_21_ME',
                          'Group 2': 'WT_AN_21_LE',
                          'Group 3': 'WT_AN_21_ES',
                          'Group 4': 'WT_AN_21_MS',
                          'Group 5': 'bridge_TMT10'},
                'TMT11': {'Group 1': 'WT_O2_28_ME',
                          'Group 2': 'WT_O2_28_ED',
                          'Group 3': 'WT_O2_28_MD',
                          'Group 4': 'WT_O2_28_MS',
                          'Group 5': 'bridge_TMT11'},
                'TMT12': {'Group 1': 'WT_O2_29_ME',
                          'Group 2': 'WT_O2_29_ED',
                          'Group 3': 'WT_O2_29_MD',
                          'Group 4': 'WT_O2_29_MS',
                          'Group 5': 'bridge_TMT12'},
                'TMT13': {'Group 1': 'WT_O2_30_ME',
                          'Group 2': 'WT_O2_30_ED',
                          'Group 3': 'WT_O2_30_MD',
                          'Group 4': 'WT_O2_30_MS',
                          'Group 5': 'bridge_TMT13'},
                'TMT14': {'Group 1': 'WT_O2_28_LE',
                          'Group 2': 'WT_O2_29_LE',
                          'Group 3': 'WT_O2_30_LE',
                          'Group 4': 'bridge_TMT14'}}

#define the technical and biological replicates
TR1             = ['126', '127C', '128C', '129C', '130C']
TR2             = ['127N', '128N', '129N', '130N', '131']
TR              = [TR1, TR2]

groups          = ['Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5']
experiments     = ['MG_AN', 'MG_O2', 'WT_AN', 'WT_O2']
TMTs            = [['TMT01', 'TMT02', 'TMT03'], ['TMT04', 'TMT05', 'TMT06', 'TMT07'], 
                   ['TMT08', 'TMT09', 'TMT10'], ['TMT11', 'TMT12', 'TMT13', 'TMT14']]
breps           = [['16', '17', '18'], ['16', '17', '18'], 
                   ['19', '20', '21'], ['28', '29', '30']]

#define phases and hours at which samples were taken
xdata_AN        = [7.5, 10.5, 13.5, 16.5]
xdata_AN        = [x - 13.5 for x in xdata_AN]
xaxis_AN        = ['ME', 'LE', 'ES', 'MS']
xdata_O2        = [6, 9, 12, 16.5, 27]
xdata_O2        = [x - 12 for x in xdata_O2]
xaxis_O2        = ['ME', 'LE', 'ED', 'MD', 'MS']
gluc_data       = [92.4, 83.5, 71.7, 44.3, 11.1, 0.06, 0.08, 0.07,   0,  0]
xdata_gluc      = [-7.5,   -6, -4.5,   -3, -1.5,    0,  1.5,    3, 4.5, 15]

#minor paralogs that were removed in the MG strain
minor_paralogs  = ['HXK1', 'GLK1', 'TDH1', 'TDH2', 'GPM2', 'GPM3', 'ENO1', 'PYK2', 
                   'PDC5', 'PDC6', 'ADH2', 'ADH4', 'ADH5']

#file and folder specifications
folder_prot_data  = 'crude_data'                        # name of data folder
folder_output     = 'Proteomics results'                # name of output folder
folderspec_output = os.getcwd() + "/" + folder_output   # path of output folder

# loop to check if the folder already exists in that path and to make the folder if it does not exist yet
if os.path.isdir(folderspec_output):
    print ("\nDirectory %s already exists in this folder." %folder_output)
else: 
    try:
        os.mkdir(folderspec_output)
    except OSError:
        print ("\nCreation of the directory %s failed" % folderspec_output)
    else:
        print ("\nSuccessfully created the directory %s " % folderspec_output)

#other output folders
folder_clustermap   = folderspec_output + '/Clustermaps_reproducibility'
folder_FC           = folderspec_output + '/FC_plots'
folder_heatmap      = folderspec_output + '/Heatmaps_CCM'
folder_STRING       = folderspec_output + '/STRING_enrichment'
folder_volcano      = folderspec_output + '/Volcano_plots'

folders = [folder_clustermap, folder_FC, folder_heatmap, folder_STRING, folder_volcano]

# loop to check if the folder already exists in that path and to make the folder if it does not exist yet
for folder in folders:
    if not os.path.isdir(folder):
        try:
            os.mkdir(folder)
        except OSError:
            print ("\nCreation of the directory %s failed" % folder)
        else:
            print ("\nSuccessfully created the directory %s " % folder)


#%% 2. Uniprot accession modifications

print("\n\033[4m2. Retrieving Uniprot gene names and/or associated GO terms\033[0m")     

#generate df containing uniprot accessions for yeast and corresponding gene names
Uniprot2Gene = apf.Uniprot_2_Gene(filename = 'Uniprot2Gene.csv', folder_spec = os.getcwd() + '/KEGG_GO_Uniprot')

#generate dictionary containing all uniprot accessions for yeast, corresponding gene names, and annotated GO terms
Uni_GO = apf.Uniprot2GO(Uniprot2Gene, filename = 'Uni2GO_yeast.txt', folder_spec = os.getcwd() + '/KEGG_GO_Uniprot')

#%% 3. import and modify the proteomics datafile

print("\n\033[4m3. Import and modify the proteomics datafile\033[0m")

# import the data and store in dictionary with their TMT number
# it is important that the prot data files are named like 'proteins_TMTX.csv'
backup      = apf.import_prot_data(folder_prot_data, label_scheme, TR, groups)

#adjust the ratios to make sure the bridging samples have a ratio of 1
backup      = apf.adjust_ratios(backup)

# Find data below the unique peptides threshold and remove it
for TMT in backup:
    RowstoDelete2 = []

    keep_cols    = [col for col in backup[TMT].columns if 'Unique' in col]
    temp         = backup[TMT][keep_cols].astype(float) 
    
    for i, row in temp.iterrows():
        if row.min() < 3:
            RowstoDelete2.append(i)
        
    backup[TMT] = backup[TMT].drop(RowstoDelete2) 
    
#merge all the proteomics data in one large dataframe
all_combined= apf.concatenate_data(backup)

# assess gene names to accessions in the df
accessions  = [re.split('\|',x)[0] for x in all_combined['Accession'].tolist() if x == x]
gene_names  = []

for accession in accessions:
    if (Uniprot2Gene['Accession'] == accession).any(0):
        k   = np.where(Uniprot2Gene['Accession'] == accession)[0][0]
        gene_names.append(Uniprot2Gene['Gene name'][k])
    else:
        print('\nNo gene name found for accession number %s' %accession)
        gene_names.append("") #make empty cell
    
all_combined.insert(2, 'Gene name', gene_names)
all_combined['Accession'] = accessions 

#make separate dataframes of the separate biological replicates add to dictionary
#also calculate the normalised intensities 
sep_data = apf.separate_data(all_combined, experiments, TMTs, breps)

#%% 4. KEGG pathways

print("\n\033[4m4. Import the KEGG pathways\033[0m")

# load the manually created pathways for the CCM from the excel file
pathway_file = os.getcwd() + '/KEGG_GO_Uniprot/Simplified Enzymes Central Carbon Metabolic Pathway.xlsx'
read_file    = pd.read_excel(pathway_file)
read_file    = read_file[[col for col in read_file.columns if not 'Unnamed' in col]]
pathways     = ['Glycolysis', 'TCA cycle', 'Anaplerotic/Gluconeogenic reactions', 
                'Alcoholic fermentation', 'Acetate metabolism', 'Glycerol metabolism', 
                'Glycogen metabolism', 'Pentose Phosphate Pathway']
KEGG_data    = {}

for pathway in pathways:
    keep_cols = [col for col in read_file.columns if pathway.lower() in col.lower()]
    
    if len(keep_cols) == 0:
        continue
    
    KEGG_data[pathway] = {}
    
    for col in keep_cols:
        if 'gene' in col:
            KEGG_data[pathway]['Gene name'] = [x for x in list(read_file[col]) if str(x) != 'nan']
            
        elif 'accession' in col:
            KEGG_data[pathway]['Accession'] = [x for x in list(read_file[col]) if str(x) != 'nan']


#open and load the saved uniprot data from the KEGG tool 
KEGG_data_all = apf.KEGG_tool(folder_spec = os.getcwd() + "/KEGG_GO_Uniprot", filename = 'KEGG_2_Uniprot_all_pathways.txt')

# find unidentified proteins in the CCM
pathway_dict = {}
for key in KEGG_data:
    pathway_dict[key] = KEGG_data[key]['Gene name']

plot_pathways = ['Glycolysis', 'TCA cycle', 'Anaplerotic/Gluconeogenic reactions', 
                 'Alcoholic fermentation', 'Acetate metabolism', 'Glycerol metabolism', 
                 'Glycogen metabolism', 'Pentose Phosphate Pathway']

gene_list = []
no_data   = {}

#find genes that are not present in any of the experiments and remove them
for exp in sep_data:
    for brep in sep_data[exp]:
        exp_data = sep_data[exp][brep]
                
        #loop over all pathways to find the selected pathways
        for pathway in plot_pathways:
            gene_list = pathway_dict[pathway]
            Out_path, no_data = apf.select_pathway(exp_data, pathway_dict[pathway], no_data)

remove_genes = []

for key, value in no_data.items():
    if value >= len(sep_data): #if the gene is absent in all the measurements, remove it
        remove_genes.append(key)

print('\nProteins not found in any of the experiments:')
for pathway in plot_pathways:
    for gene in remove_genes:
        if gene in pathway_dict[pathway]:
            print('\t* '+ pathway, gene)
            indx = pathway_dict[pathway].index(gene)
            del pathway_dict[pathway][indx]

#define the list of genes you want to compare the data on
glycolysis_list = KEGG_data['Glycolysis']['Gene name']+KEGG_data['Alcoholic fermentation']['Gene name']
MG_genenames    = [gene for gene in glycolysis_list if gene not in minor_paralogs]

#%% 5. make clustermap to check reproducibility of technical and biological replicates

print("\n\033[4m5. Make clustermap to check reproducibility of technical and biological replicates\033[0m")

#iterate over all the separated biological replicate data
for exp in sep_data:
    
    #concatenate the data from all the biological replicates under the same conditions
    for brep in sep_data[exp]: 
        concat_data = sep_data[exp][brep]
        break
    
    for brep in sep_data[exp]:
        keep_cols = sep_data[exp][brep].columns.difference(concat_data.columns)
        concat_data = pd.concat([concat_data, sep_data[exp][brep][keep_cols]], axis = 1)
                            
    concat_data = concat_data.dropna()  #drop NaN values, for better overview in the image
    
    #for better visualisation of the data, normalise to the avaerage of the three ME phases
    concat_data.loc[:, [col for col in concat_data.columns if 'Ratio' in col and 'TR' in col]] = concat_data[[col for col in concat_data.columns if 'Ratio' in col and 'TR' in col]].div(concat_data[[col for col in concat_data.columns if 'ME' in col and 'Ratio' in col and 'TR' in col]].mean(axis=1), axis =0)
    concat_data = concat_data.drop([col for col in concat_data.columns if 'Ratio' in col and 'TR' not in col], axis=1)
    
    #plot clustermaps
    xaxis  = [col for col in concat_data.columns if 'Ratio' in col and 'TR' in col]
    xaxis  = [re.split('Ratio', item)[1] for item in xaxis]
    xaxis  = [re.split(exp+'_', item)[1] for item in xaxis]
    folder = folder_clustermap
            
    figname = 'Clustermap_techreps_' + exp +'.pdf'
    
    if 'MG' in exp:
        title1 = 'MG '
    else:
        title1 = 'WT '
    if 'AN' in exp:
        title2 = 'anaerobic'
    else:
        title2 = 'aerobic'
    figtitle = 'Clustermap technical and biological replicates ' + title1 + title2
    
    c = ["lightcoral", "red","darkred","black", "darkgreen","green", "palegreen"]
    v = [0,.2, .35, .5,0.65, 0.8, 1.]
    l = list(zip(v,c))
    cmap=LinearSegmentedColormap.from_list('rg',l, N=51)
    
    clustergrid = apf.make_clustermap(concat_data, folder=folder, figname=figname, gene_labels = 'Gene name',
                                      figtitle=None, xaxis=xaxis, 
                                      show=True, savefig=True, col_cluster=True, 
                                      row_cluster=False, adjust_ratio = False, barplot = False, 
                                      figsize=(9,12), yticklabels=False, cmap=cmap, center=0, 
                                      dendrogram_ratio=(.1, .13), vmin = -2, vmax =3)



#%% 6. average the biological replicates

print("\n\033[4m6. Calculate the averages and standard deviations of the biological replicates\033[0m")

avg_data = {}

for exp in sep_data:
    
    #concatenate the data from all the biological replicates under the same conditions
    for brep in sep_data[exp]: 
        sep_data[exp][brep].index = sep_data[exp][brep]['Gene name']
        concat_data = sep_data[exp][brep]
        break
    
    for brep in sep_data[exp]:
        sep_data[exp][brep].index = sep_data[exp][brep]['Gene name']
        keep_cols = sep_data[exp][brep].columns.difference(concat_data.columns)
        concat_data = pd.concat([concat_data, sep_data[exp][brep][keep_cols]], axis = 1)
    
    # select the correct sample times and phases for the conditions
    if 'AN' in exp:
        xdata = xdata_AN
        xaxis = xaxis_AN
    
    elif 'O2' in exp:
        xdata = xdata_O2
        xaxis = xaxis_O2
    
    new_vals = pd.DataFrame()
    
    # calculate intensity and ratio means, standard deviations and number of identificantions
    for phase in xaxis:
        keep_cols           = [col for col in concat_data.columns if 'Corrected Intensity' in col and phase in col]
        
        corrected_intensity = concat_data[keep_cols].mean(axis=1).to_frame()
        corrected_intensity.columns = ['Corrected Intensity ' + exp +'_'+phase + ' brep avg']
        
        keep_cols           = [col for col in concat_data.columns if 'Ratio' in col and 'avg' in col and phase in col]
        new_ratios          = concat_data[keep_cols].mean(axis=1).to_frame()
        new_ratios.columns  = ['Ratio ' + exp +'_'+phase + ' brep avg']
        sd                  = concat_data[keep_cols].std(axis=1).to_frame()
        sd.columns          = ['sd ' + exp +'_'+phase + ' brep avg']
        
        keep_cols           = [col for col in concat_data.columns if 'Ratio' in col and phase in col and 'avg' in col and not 'brep' in col]
        num                 = (3 - concat_data[keep_cols].isna().sum(axis=1)).to_frame()
        num.columns         = ['num ' + exp +'_'+phase + ' brep avg']
        new_vals            = pd.concat([new_vals, corrected_intensity, new_ratios, sd, num], axis=1)
    
    new_genenames = []
    insig         = []
    
    for gene in list(concat_data['Gene name']):
        if gene in list(concat_data['Gene name'][concat_data[keep_cols].isnull().sum(axis=1)==2]):
            insig.append(gene)
            new_genenames.append(gene+'**')
        elif gene in list(concat_data['Gene name'][concat_data[keep_cols].isnull().sum(axis=1)==1]):
            new_genenames.append(gene+'*')
        else:
            new_genenames.append(gene)        
    
    # clean up de dataframes, remove unnecessary columns
    drop_cols                       = [col for col in concat_data.columns if 'Ratio' in col or ('RI' in col) or 'FC' in col or 'TR' in col or 'Intensity' in col]
    concat_data                     = concat_data.drop(drop_cols, axis=1)
    
    # append the new columns to the dataframe
    concat_data                     = pd.concat([concat_data, new_vals], axis =1)
    concat_data['Gene name sign']   = new_genenames
    
    # make the index of the dataframe equal to the gene names
    concat_data.index               = concat_data['Gene name']
    
    #store the data in a dict
    avg_data[exp]       = concat_data
    

#%% 7. make clustermap of average breps

print("\n\033[4m8. Make clustermap of the CCM proteins of the averaged biological replicates to analyse the overall trends\033[0m")

for exp in avg_data:
    
    if 'AN' in exp:
        xdata = xdata_AN
        xaxis = xaxis_AN
    elif 'O2' in exp:
        xdata = xdata_O2
        xaxis = xaxis_O2

    #give printed line as output to visualise progress of script
    print('\nGenerating heatmaps for %s...' %exp)
    
    exp_data          = avg_data[exp]
    Out_path_combined = pd.DataFrame()
    
    #loop over all pathways and find the proteomics data for each protein in these pathways
    for pathway in plot_pathways:
        
        #retrieve genes in pathway
        gene_list           = pathway_dict[pathway] 
        
        #retrieve prot data for these genes and annotate the pathway name in the df
        Out_path, no_data   = apf.select_pathway(exp_data, pathway_dict[pathway], no_data)
        Out_path['Pathway'] = pathway
        
        #if there are no recognised proteins in the pathway, skip to nextpathway
        if len(Out_path) == 1:
            continue
        
        #concatenate the dataframe to create large dataframe with all data to be plotted
        Out_path_combined   = pd.concat([Out_path_combined, Out_path])
    
    #select the data needed for the heatmap
    keep_cols    = [col for col in Out_path_combined.columns if 'Ratio' in col and 'avg' in col or 'Gene name' in col or 'Pathway' in col]
    cluster_data = Out_path_combined[keep_cols]
    
    #print CCM proteins that were not detected in those experiments
    print('\nProteins not found in the %s experiments:' %exp)
    for i, item in enumerate(list(cluster_data['Gene name'][cluster_data.isnull().sum(axis=1)>1])):
        print('\t* '+ str(cluster_data['Pathway'][cluster_data.isnull().sum(axis=1)>1].iloc[i]), item)
    
    #for a better overview in the heatmap, remove the undetected proteins
    cluster_data.index  = cluster_data['Gene name']
    drop_data           = cluster_data[cluster_data.isnull().sum(axis=1)>1]
        
    Out_path_combined.index = Out_path_combined['Gene name']            
    Out_path_combined       = Out_path_combined.drop(drop_data.index)
            
    #plot clustermaps
    folder   = folder_heatmap                
    figname  = 'Clustermap_avg_CCM_' + exp +'.pdf'
    
    if 'MG' in exp:
        title1 = 'MG '
    else:
        title1 = 'WT '
    if 'AN' in exp:
        title2 = 'anaerobic'
    else:
        title2 = 'aerobic'
    
    figtitle = 'Heatmap and intensity barplot CCM ' + title1 + title2

    
    c    = ["lightcoral", "red","darkred","black", "darkgreen","green", "palegreen"]
    v    = [0,.2, .35, .5,0.65, 0.8, 1.]
    l    = list(zip(v,c))
    cmap = LinearSegmentedColormap.from_list('rg',l, N=51)

    clustergrid = apf.make_clustermap(Out_path_combined, folder=folder, figname=figname, 
                                      figtitle=figtitle, xaxis=xaxis, 
                                      show=True, savefig=True, col_cluster=False, 
                                      row_cluster=False, adjust_ratio = True, barplot = True, 
                                      figsize=(17,20), yticklabels=1, cmap=cmap, center=0, 
                                      row_colors='row_colors', dendrogram_ratio=(.1, .1), vmin= -3, vmax = 3)

#%% 8. calculate p-values between selected experiments (here: MG vs WT)

print("\n\033[4m9a. Calculate the p-values for the difference in means between the MG and WT data\033[0m")

exp_data_AN = pd.DataFrame()
exp_data_O2 = pd.DataFrame()

for exp in avg_data:
    
    #if by mistake, minor isozymes were detected in the MG-strain, remove these (only if first PEAKS is manually analysed!!)
    if 'MG' in exp:
        for gene in minor_paralogs:
            if gene in avg_data[exp].index:
                avg_data[exp].loc[gene, [col for col in avg_data[exp].columns if 'Accession' not in col and 'Gene name' not in col]] = np.nan

    exp_data = avg_data[exp]
    
    #combine the anaerobic MG and WT experiments and the aerobic MG and WT experiments
    if 'AN' in exp:
        keep_cols   = avg_data[exp].columns.difference(exp_data_AN.columns)
        exp_data_AN = pd.concat([exp_data_AN, exp_data[keep_cols]], axis =1)
    
    elif 'O2' in exp:
        keep_cols   = avg_data[exp].columns.difference(exp_data_O2.columns)
        exp_data_O2 = pd.concat([exp_data_O2, exp_data[keep_cols]], axis =1)

data = {'AN': exp_data_AN,
        'O2': exp_data_O2}

for exp in data:
    
    #select the correct sample times for the condition
    if 'AN' in exp:
        xaxis = xaxis_AN        
    elif 'O2' in exp:
        xaxis = xaxis_O2
    
    exp_data   = data[exp]
    
    #compare WT and MG for each phase and calculate the p-value for the difference in means        
    for phase in xaxis:
           
        mean1       = [col for col in exp_data.columns if 'Ratio' in col and phase in col and 'MG' in col][0]
        mean2       = [col for col in exp_data.columns if 'Ratio' in col and phase in col and 'WT' in col][0]
        sd1         = [col for col in exp_data.columns if 'sd' in col and phase in col and 'MG' in col][0]
        sd2         = [col for col in exp_data.columns if 'sd' in col and phase in col and 'WT' in col][0]
        n1          = [col for col in exp_data.columns if 'num' in col and phase in col and 'MG' in col][0]
        n2          = [col for col in exp_data.columns if 'num' in col and phase in col and 'WT' in col][0]
        
        p_val = st.ttest_ind_from_stats(exp_data[mean1], exp_data[sd1], exp_data[n1], exp_data[mean2], exp_data[sd2], exp_data[n2], equal_var=False, alternative='two-sided')[1]
        p_val = pd.DataFrame(p_val)
        p_val.columns = ['p-value MGvsWT '+exp+'_'+phase]
        p_val.index = exp_data.index
        
        #add the p-value column to the correct dataframes
        for key in avg_data:
            if exp in key:
                avg_data[key]['p-value MGvsWT '+exp+'_'+phase] = p_val


#%% 9. FC plots glycolysis with significance score MG vs WT

print("\n\033[4m10a. Make fold change plots of the glycolytic proteins to compare MG and WT\033[0m")

#define the list of genes you want to compare the data on
gene_list           = glycolysis_list

#remove minor paralogs
genelist            = [ elem for elem in gene_list if elem not in minor_paralogs ]

Out_gene_combined   = {}

fig, axs = plt.subplots(math.ceil(len(genelist)/3),3, figsize = (7*3, 6*math.ceil(len(genelist)/3)))

i = 0
j = 0

marker      = ['o', 'o', 's', 's']
fillstyle   = ['full', 'none', 'full', 'none']
color       = plt.rcParams['axes.prop_cycle'].by_key()['color']

#iterate over the genes in the selected gene list
for k, gene in enumerate(genelist):
    
    #define the position of the subplot to fill (j=column, i=row)
    if j  == 3:
        j  = 0
        i += 1
    
    Out_gene_combined[gene] = {}
    Out_gene                = pd.DataFrame()
    color_i                 = 0
    leg                     = []
    
    #iterate over the different experiments (MG/WT (an)aerobic)
    for exp in avg_data:
        if 'AN' in exp:
            xaxis    = xaxis_AN
            xdata    = xdata_AN
        
        elif 'O2' in exp:
            xaxis    = xaxis_O2
            xdata    = xdata_O2
        
        #check if the gene was identified in in the TMT exp
        if gene in list(avg_data[exp]['Gene name']):
            
            #retrieve the ratios for the gene
            out      = avg_data[exp].loc[avg_data[exp].index == gene]
            out_plot = out[[col for col in out.columns if 'Ratio' in col]]
            Out_gene = pd.concat([Out_gene, out])
            #retrieve the standard deviations for the error bars
            yerr     = out[[col for col in out.columns if 'sd' in col]]
            yerr     = [i[0] for i in yerr.T.to_numpy()]
            
            #different plotting lines if the avg data is based on less data points
            if '**' in str(out['Gene name sign']):
                fmt  = marker[color_i]+':'
            elif '*' in str(out['Gene name sign']):
                fmt  = marker[color_i]+'-.'
            else: 
                fmt  = marker[color_i]+'-'
            
            #only plot if there are more than 3 data points for that gene
            if out_plot.isnull().sum().sum() < 3:
                
                #make a list of amount of detected unique peptides to show in graph legend
                unique = list(out[[col for col in out.columns if 'Unique' in col]].iloc[0])
                unique = [int(un) for un in unique if un==un]
                
                label  = exp
                lns    = axs[i, j].errorbar(xdata, list(out_plot.iloc[0]), fmt = fmt, fillstyle = fillstyle[color_i], yerr = yerr, color = color[color_i], linewidth = 3, markersize =11, label=exp+ ' ' + str(unique), capsize = 5, elinewidth = 1.2, alpha=0.7, zorder=2)
                leg.append(lns)
                
                #annotate the significance level of the difference between MG and WT with 1, 2 or 3 *'s
                for l, phase in enumerate(xdata):
                    annot = list(out_plot.iloc[0])[l]
                    if float(avg_data[exp][[col for col in avg_data[exp].columns if 'MGvsWT' in col 
                                            and xaxis[l] in col]].loc[gene]) < 0.001:
                        axs[i, j].annotate('***', xy = (xdata[l], annot), xytext = (xdata[l]+0.3, annot), 
                                           color = 'k', fontsize = 15)
                    elif float(avg_data[exp][[col for col in avg_data[exp].columns if 'MGvsWT' in col 
                                              and xaxis[l] in col]].loc[gene]) < 0.01:
                        axs[i, j].annotate('**', xy = (xdata[l], annot), xytext = (xdata[l]+0.3, annot), 
                                           color = 'k', fontsize = 15)
                    elif float(avg_data[exp][[col for col in avg_data[exp].columns if 'MGvsWT' in col 
                                              and xaxis[l] in col]].loc[gene]) < 0.05:
                        axs[i, j].annotate('*', xy = (xdata[l], annot), xytext = (xdata[l]+0.3, annot), 
                                           color = 'k', fontsize = 15)
                    elif float(avg_data[exp][[col for col in avg_data[exp].columns if 'MGvsWT' in col 
                                              and xaxis[l] in col]].loc[gene]) < 0.1:
                        axs[i, j].annotate('^', xy = (xdata[l], annot), xytext = (xdata[l]+0.3, annot), 
                                           color = 'k', fontsize = 12)
            
            color_i += 1 #move to the next color
    
    #plot the glucose concentration on the secondary y-axis        
    ax2 = axs[i,j].twinx()
    lns = ax2.errorbar(xdata_gluc, gluc_data, fmt = '--.', markersize = 7, color='grey', zorder=1, alpha=0.3, label= 'glucose (mM)', linewidth=5)
    ax2.set_ylabel('Glucose concentration (mM)', fontsize = 18)
    ax2.tick_params(axis='y', colors='grey', labelsize=16)
    ax2.yaxis.label.set_color('grey')
    ax2.set_ylim(bottom=-1)
    
    #retrieve the labels for in the legend
    leg.append(lns)
    labs = [l.get_label() for l in leg]
    
    #adjust characteristics of the subplot        
    axs[i, j].set_title(gene.title(), fontsize = 20, fontweight = 'bold')
    axs[i, j].legend(leg, labs, fontsize = 16)
    axs[i, j].set_xticks([int(item) for item in xdata_AN]+xdata_O2)
    axs[i, j].set_xticklabels([int(item) for item in xdata_AN]+xdata_O2, fontsize = 16)
    axs[i, j].tick_params(axis='y', labelsize = 16)
    axs[i, j].set_xlabel('Time (h)', fontsize = 18)
    axs[i, j].set_ylabel('Fold change', fontsize = 18)
    axs[i, j].grid(True, linewidth=0.2, color='lightgrey', alpha=0.6)
    axs[i, j].set_xlim(right=20, left=-7)
    axs[i, j].set_ylim(top=3.8, bottom=0.5)
    
    #move to the next subplot
    j += 1
            
    Out_gene_combined[gene] = Out_gene

#save the figure
fig.tight_layout(w_pad=2)
folder  = folder_FC  
figname = 'FC_plots_glycolysis_significance_WT_MG.pdf'
fig.savefig(folder + '/' + figname, bbox_inches = 'tight')   

#%% 10. STRING enrichment and volcano plot compare WT to MG

print("\n\033[4m11a. STRING functional enrichment analysis comparing the MG and WT strain\033[0m")

# find accessions of up- and downregulated genes to link to GO enrichment analysis
all_Changes     = {}
all_Changes_g   = {}
folder          = folder_STRING
dic_up          = {}
dic_down        = {}

#specify the log2 fold change threshold you want to apply
lfc_thr         = (0.32, 0.32)

exp_data_AN = pd.DataFrame()
exp_data_O2 = pd.DataFrame()

#concatenate the WT and MG data of the (an)aerobic experiments for normalisation
for exp in avg_data:

    exp_data = avg_data[exp]
    
    if 'AN' in exp:
        keep_cols   = avg_data[exp].columns.difference(exp_data_AN.columns)
        exp_data_AN = pd.concat([exp_data_AN, exp_data[keep_cols]], axis =1)
    
    elif 'O2' in exp:
        keep_cols   = avg_data[exp].columns.difference(exp_data_O2.columns)
        exp_data_O2 = pd.concat([exp_data_O2, exp_data[keep_cols]], axis =1)

data = {'AN': exp_data_AN,
        'O2': exp_data_O2}

for exp in data:
    
    #create placeholders
    Changes   = {}
    Changes_g = {}
    new_cols  = []
    
    #select correct sample times
    if 'AN' in exp:
        xaxis = xaxis_AN
    elif 'O2' in exp:
        xaxis = xaxis_O2
        
    #retrieve a set of background genes to give as input to the STRING enrichment tool
    exp_data              = data[exp]
    string_background     = exp_data.copy()
    string_background     = string_background.dropna(subset=[col for col in exp_data.columns if 'Unique' in col][:3], how = 'all')
    string_background     = string_background.dropna(subset=[col for col in exp_data.columns if 'Unique' in col][3:], how = 'all')
    all_string_identifier = apf.backgroundgene_2_string(list(string_background.index))

    #adjust ratios
    keep_cols  = [col for col in exp_data.columns if 'Ratio' in col and 'avg' in col]
    new_ratios = exp_data[keep_cols]
    
    #normalize each phase to the WT data     
    for phase in xaxis:
        col_ratio   = [col for col in keep_cols if phase in col]
        new_ratios.loc[:,col_ratio] = new_ratios[col_ratio].div(exp_data[[col for col in keep_cols if phase in col and 'WT' in col][0]], axis=0)

    #save recalculated ratios and the log2 of the ratios to the dataframe
    exp_data[keep_cols] = new_ratios
    log2FC              = np.log2(exp_data[keep_cols].astype(float))    
    for col in log2FC.columns:
        new_cols.append(col.replace('Ratio', 'Log2 FC'))    
    exp_data[new_cols] = log2FC
        
    # determine sign up- and downregulated genes for each phase, and make volcano plots
    for phase in xaxis:
        lfc = [col for col in exp_data.columns if 'Log2 FC' in col and 'MG' in col and phase in col][0]
        pv  = [col for col in exp_data.columns if 'p-value' in col and 'MGvsWT' in col and phase in col][0]
        equal, upregulated, downregulated = apf.up_down_regulated(exp_data, lfc = lfc, pv=pv, lfc_thr=lfc_thr)
        upregulated_g     =   upregulated['Gene name'].tolist()
        downregulated_g   = downregulated['Gene name'].tolist()
        upregulated       =   upregulated['Accession'].tolist()
        downregulated     = downregulated['Accession'].tolist()
        Changes[phase]    = {'upregulated': upregulated, 'downregulated': downregulated}
        Changes_g[phase]  = {'upregulated': upregulated_g, 'downregulated': downregulated_g}
        
        #remove genes detected in less than 2 of the bio replicates        
        exp_data2 = exp_data.copy()
        exp_data2 = exp_data2.dropna(subset=[pv]) 
        
        #give printed output
        print('\nMG vs WT upregulated genes under ' +exp+ ' conditions in phase '+phase+': '+ str(upregulated_g))
        print('MG vs WT downregulated genes under ' +exp+ ' conditions in phase '+phase+': '+ str(downregulated_g))
        print('\n'+'\t'.join(['phase', '#prot', '#sign', '#upreg', '#downreg']))
        print('\t\t'.join([phase, str(len(exp_data2)), str(len(equal)+len(upregulated)+len(downregulated)), str(len(upregulated)), str(len(downregulated))]))

        # make and save volcano plot
        figname        = 'Volcano_plot_MGvsWT_' + exp +'_'+ phase+'.pdf' 
        figname        = folder_volcano+'/'+figname
        
        if exp == 'O2':
            title1 = 'aerobic'
        else:
            title1 = 'anaerobic'
            
        figtitle = title1 + ': ' + phase
        apf.gene_exp.volcano(df=exp_data2, lfc = lfc, pv = pv, lfc_thr=lfc_thr, figname=figname, title = figtitle, show=True, geneid = ['Gene name'], genenames='deg', gstyle=1, gfont=13, dotsize = 15, lfc_thr2=(0.41,0.32), dim=(6.5,6), axlabelfontsize=17, plotlegend=True) 
        
        #find the occurrences of the up- and downregulated genes
        for gene in upregulated_g:
                if gene not in dic_up.keys():
                    dic_up[gene] = 1
                else:
                    dic_up[gene] += 1
        
        for gene in downregulated_g:   
                if gene not in dic_down.keys():
                    dic_down[gene] = 1
                else:
                    dic_down[gene] += 1
        
    all_Changes[exp]    = Changes
    all_Changes_g[exp]  = Changes_g
    filename            = folder + '/string_MGvsWT_' + exp
    
    print1              = 'MG ' + exp
    print2              = 'WT ' + exp
    
    # find and print enriched and depleted GO terms for each phase
    apf.STRING(Changes, xaxis, print1, print2, filename, background=all_string_identifier)

#%% A. FC plots additional proteins of interest

print("\n\033[4mA. Make fold change plots for additional proteins of interest\033[0m")

#input the gene names of interest in one of the following three ways:
    
#via KEGG pathways eg:
genelist = KEGG_data['Glycogen metabolism']['Gene name']
genelist = KEGG_data_all['sce00240']['Gene name']

#via certain keywords eg:
genelist = [ elem for elem in avg_data['WT_AN'].index if 'HSP' in elem]

#manually give gene names as input eg the ROS proteins from bisschops:
genelist = ["SOD1", 'SOD2', 'CTA1']


