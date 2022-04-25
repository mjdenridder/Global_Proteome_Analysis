"""
Various functions for processing large scale proteomics data. 

Author: Wiebeke van den Brandeler

December 2021
"""

#%% Function to print progress of a loop, script or function
def progressBar(name, value, endvalue, bar_length = 33, width = 20):
    
    #import required modules
    import sys 
    
    percent = float(value) / endvalue
    arrow = '-' * int(round(percent*bar_length) - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write("\r{0: <{1}} : [{2}]{3}%".format(\
                     name, width, arrow + spaces, float("{:.1f}".format(percent*100))))
        
    sys.stdout.flush()
    if value == endvalue:     
          sys.stdout.write('\n\n')

#%% Uniprot accession to gene name converter
import os

def Uniprot_2_Gene(url = 'https://www.uniprot.org/docs/yeast.txt', firstname = 'AAC1', filename = 'Uniprot2Gene.csv', folder_spec = os.getcwd() + '/KEGG_GO_Uniprot')    :
    """
    

    Parameters
    ----------
    url : STRING. URL where the full genome of the selected species can be found with
        their corresponding accession numbers. Default = S. cerevisiae :
        url = 'https://www.uniprot.org/docs/yeast.txt'
        
    firstname : STRING. The first gene name that is named on the website, so 
        the start position of the data is known. The default is 'AAC1'.
    
    filename : STRING. How to save the Uniprot 2 gene dataframe. The default is 
    'Uniprot2Gene.csv'.
    
    folder_spec: STRING. Where to save the Uniprot 2 gene dataframe. The default
    is os.getcwd() + '/KEGG_GO_Uniprot'

    Returns
    -------
    df : dataframe containing gene names and corresponding accession numbers
        of S. cerevisiae

    """
    
    #import required modules 
    import urllib
    import pandas as pd
    import os
    
    print('\nInitializing the Uniprot 2 Gene tool...\n')
    
    # loop to check if the folder already exists in that path and to make the folder if it does not exist yet
    if os.path.isdir(folder_spec):
        print ("Directory %s already exists." %folder_spec)
    else: 
        try:
            os.mkdir(folder_spec)
        except OSError:
            print ("Creation of the directory %s failed" % folder_spec)
        else:
            print ("Successfully created the directory %s " % folder_spec)
    
    filename = folder_spec + '/' +filename
    
    # check if the file name already exists and if so, import the existing uniprot to GO data 
    if os.path.isfile(filename):
        print("File %s exists" %filename)
        df = pd.read_csv(filename)
    
    else:        
        item = urllib.request.urlopen(url).read()                   #read url
        gene_accession = []                                         #create placeholder
        
        # decode the url and split into readable output
        read_output = item.decode("utf-8").split('\n')                 
        for i, item in enumerate(read_output):
            if firstname in item:
                read_output = read_output[i:]
        
        for i, item in enumerate(read_output):
            if item == '':
                read_output = read_output[:i]
        
        for i, item in enumerate(read_output):
            
            #print progress
            progressBar("Annotating gene names to all accessions",i, len(read_output) - 1)
            
            split_output = item.split(' ')
            split_output = [i for i in split_output if i != '']
            rem_item     = []
            for j, item2 in enumerate(split_output):
                if ';' in item2:
                    split_output[j] = split_output[j].replace(";", "")
                    rem_item.append(j+1)
                else:
                    break
            if len(rem_item )==1:
                del split_output[rem_item[0]]
            elif len(rem_item) > 1:
                del split_output[rem_item[0]:(rem_item[-1]+1)]
                
            #append gene name to placeholder    
            gene_accession.append([split_output[0], split_output[2]])
        
        # store data in dataframe and save to .csv file
        df = pd.DataFrame(gene_accession, columns=['Gene name', 'Accession'])
        df.to_csv(filename, index = False)
    
    return df    

#%% Annotate Uniprot accession numbers with all related GO terms

import os

def Uniprot2GO(df, filename = 'Uni2GO_yeast.txt', folder_spec = os.getcwd() + '/KEGG_GO_Uniprot', col_Uniprot = ['go(biological%20process)', 'go(molecular%20function)', 'go(cellular%20component)'], colname = ['biological process', 'molecular function', 'cellular component'], print_nogo = False):
    """
    This function takes all the uniprot accession numbers from the proteins 
    identified in the experiments and associates them with GO-terms. 
    
    Due to technical errors occurring by repeatedly accessing the Uniprot site,
    the script may stop and return the error: 'Connection refused'. When this
    happens, run the script over and over again untill 100% completion.
    
    Another problem is that from time to time, by accident a different
    accession number is taken from the website than the accession number that
    was given as input. To make sure the correct accession number was processed,
    containing the correct GO-terms, another loop was built in.

    Parameters
    ----------
    df : datafrane containing the gene names and accession numbers of the whole 
        yeast genome
    filename : STRING, optional
        name of the Uni2Go txt file. default = 'Uni2GO_yeast.txt'
    folder_spec : STRING, optional
        specify the directory. default = os.getcwd() + '/Uni2GO'
    col_Uniprot : LIST OF STRINGS, optional
        specify the type of GO terms you want to analyse, to give as input in the 
        URL. default = ['go(biological%20process)', 'go(molecular%20function)', 
                        'go(cellular%20component)']
    colname : LIST OF STRINGS, optional
        specify the names of the type of GO terms you want to analyse.
        default = ['biological process', 'molecular function', 'cellular component']
    print_nogo : BOOL, optional
        choose True if you want a printed list of accession numbers that did not
        have any found annotated GO terms. default = False

    Returns
    -------
    Uni_GO : dictionary containing the accession numbers 
        with the GO-ids and GO-names.

    """
    
    #import required modules
    import json
    import urllib
    import os
    
    print('\nInitializing the Uniprot2GO tool...\n')
    
    # create empty dictionary
    Uni_GO = {}
    # make list of the Uniprot accessions of all proteins found in the experiments
    Uniprot_accessions = list(df['Accession'])
    
    # loop to check if the folder already exists in that path and to make the folder if it does not exist yet
    if os.path.isdir(folder_spec):
        print ("Directory %s already exists." %folder_spec)
    else: 
        try:
            os.mkdir(folder_spec)
        except OSError:
            print ("Creation of the directory %s failed" % folder_spec)
        else:
            print ("Successfully created the directory %s " % folder_spec)
    
    filename = folder_spec + '/' +filename
        
    # check if the file name already exists and if so, import the existing uniprot to GO data 
    if os.path.isfile(filename):
        print ("File %s already exists, import the data\n" %filename)
        with open(filename, "r") as fp: 
            Uni_GO = json.load(fp)
            
            
            """
            Uncomment the section below if the Uniprot to GO data is not fully
            complete yet. This is commented out because some accessions have no
            associated GO terms and are not included in the dict, so the output 
            data is never the same length as the input data. This causes the 
            script to run from approx 96% for each new csv file, which is not a 
            problem but this makes the script longer.
            """
            
            # ----------------------------------- start comment below this line
            
            # # check if the complete data is imported
            # if len(Uni_GO) == len(Uniprot_accessions):
                
            #     # a check to see whether all the uniprot accessions from the experiments were associated with GO terms        
            #     for item in Uniprot_accessions:
            #         if item not in Uni_GO.keys():
            #             print('%s is not in the dict. reason: no go-term' %item)
            
            # # if data is not complete, keep loading in the data
            # else:
            #     start = len(Uni_GO)-1                                                             # the start position to continue loading in the data is the length of the already existing data
                
            #     # import the GO terms for the remaining Uniprot accessions starting from the start position
            #     for j, item in enumerate(Uniprot_accessions[start:]):
                
            #         progressBar("Annotating GO terms to all accessions",j+1+start, len(Uniprot_accessions) - 1)
                    
            #         time.sleep(0.2) #use this to try not overload the uniprot server
                    
            #         GOname = []
            #         GOid = []
                    
            #         for i, col in enumerate(col_Uniprot):
            #             # define urls 
            #             url_left='https://www.uniprot.org/uniprot/?query='
            #             url_right='&sort=score&columns=id,' + col + ',&format=tab'
                        
            #             #Retrieve Uniprot data
            #             url_full = url_left+item+url_right                                      # full URL to extract info from example: https://www.uniprot.org/uniprot/?query=YGL253W&sort=score&columns=id,genes(PREFERRED),&format=tab
            #             item2 = urllib.request.urlopen(url_full).read()                         # read url
            #             read_output = item2.decode("utf-8").split('\n')                         # decode the url to make string from bytes and split at enters ('\n')
            #             split_output = read_output[1].split('\t')                               # split the accession number from the GO terms by splitting at tabs ('\t')
                        
            #             # check whether the correct accession number was retrieved from the site and else, iterate over the website lines to find the correct one
            #             if split_output[0] != item:
            #                 for item3 in read_output:
            #                     if item3.split('\t')[0] == item:
            #                         #print('Changing output line of %s' %item)
            #                         split_output = item3.split('\t')
            #                         break     
                        
            #             #if there are no GO terms associated with this accession number, continue to the next accession
            #             if split_output[1] == '':
            #                 continue
                        
            #             GO_term = split_output[1].split('; ')                           # split the GO terms in separate strings by splitting at ;
                        
            #             #append the GO-ids and GO-names to a list 
            #             GO_name = []
            #             GO_id   = []
                        
            #             for term in GO_term:
            #                 GO_name.append(term.split(' [')[0])
            #                 GO_id.append(term.split(' [')[1].split(']')[0])
                        
                        
            #             GOname.append([colname[i], GO_name])
            #             GOid.append([colname[i], GO_id])
                        
            #             # save the information in a dictionary
            #             Uni_GO[split_output[0]] = { 'Gene name' : df['Gene name'].loc[j],
            #                                         colname[i]: {'GO-name': GO_name, 
            #                                                     'GO-id'  : GO_id  }}
                        
            #             # save and publish the Uniprot_data list to a .txt file, for later access
            #             with open(filename, "w") as fp:                                         
            #                 json.dump(Uni_GO, fp)
                
            
            # ----------------------------- stop the commenting out from here
            
            
            # a check to see whether all the uniprot accessions from the experiments were associated with GO terms        
            if print_nogo:
                for item in Uniprot_accessions:
                    if item not in Uni_GO.keys():
                        print('%s is not in the dict. reason: no go-term' %item)
                
    # if the file name does not exist yet in the current directory, start from scratch
    else:
        print ("File %s does not exist, creating the file may take a while (up to 1,5 hours)...\n" %filename)    
        
        # import the GO terms for the Uniprot accessions
        for j, item in enumerate(Uniprot_accessions):
            
            progressBar("Annotating GO terms to all accessions",j, len(Uniprot_accessions) - 1)
            
            #time.sleep(0.2) #use this to try not overload the uniprot server
            
            GOname = []
            GOid = []
            
            for i, col in enumerate(col_Uniprot):
                # define urls 
                url_left='https://www.uniprot.org/uniprot/?query='
                url_right='&sort=score&columns=id,' + col + ',&format=tab'
                
                #Retrieve Uniprot data
                url_full = url_left+item+url_right                                      # full URL to extract info from example: https://www.uniprot.org/uniprot/?query=YGL253W&sort=score&columns=id,genes(PREFERRED),&format=tab
                item2 = urllib.request.urlopen(url_full).read()                         # read url
                read_output = item2.decode("utf-8").split('\n')                         # decode the url to make string from bytes and split at enters ('\n')
                split_output = read_output[1].split('\t')                               # split the accession number from the GO terms by splitting at tabs ('\t')
                
                # check whether the correct accession number was retrieved from the site and else, iterate over the website lines to find the correct one
                if split_output[0] != item:
                    for item3 in read_output:
                        if item3.split('\t')[0] == item:
                            #print('Changing output line of %s' %item)
                            split_output = item3.split('\t')
                            break     
                
                #if there are no GO terms associated with this accession number, continue to the next accession
                if split_output[1] == '':
                    continue
                
                GO_term = split_output[1].split('; ')                           # split the GO terms in separate strings by splitting at ;
                
                #append the GO-ids and GO-names to a list 
                GO_name = []
                GO_id   = []
                
                for term in GO_term:
                    GO_name.append(term.split(' [')[0])
                    GO_id.append(term.split(' [')[1].split(']')[0])
                
                
                GOname.append([colname[i], GO_name])
                GOid.append([colname[i], GO_id])
            
            Uni_GO[split_output[0]] = { 'Gene name' : df['Gene name'].loc[j]}
            
            for i, col in enumerate(GOname):
                
                # save the information in a dictionary
                Uni_GO[split_output[0]][GOname[i][0]] = {'GO-name': GOname[i][1], 
                                                         'GO-id'  : GOid[i][1]  }
            
            # save and publish the Uniprot_data list to a .txt file, for later access
            with open(filename, "w") as fp:                                         
                json.dump(Uni_GO, fp)
        
        # a check to see whether all the uniprot accessions from the experiments were associated with GO terms        
        if print_nogo:
            for item in Uniprot_accessions:
                if item not in Uni_GO.keys():
                    print('%s is not in the dict. reason: no go-term' %item)

    return Uni_GO

#%% KEGG tool

import os

def KEGG_tool(folder_spec = os.getcwd() + "/KEGG_GO_Uniprot", filename = 'KEGG_2_Uniprot_all_pathways.txt', organism = 'sce'):
    """
    

    Parameters
    ----------
    folder_spec : STRING, optional
        DESCRIPTION. The default is os.getcwd() + "/KEGG_GO_Uniprot".
    filename : STRING, optional
        DESCRIPTION. The default is 'KEGG_2_Uniprot_all_pathways.txt'.
    organism : STRING, optional
        organism you want to list the pathways for. The default is 'sce'.

    Returns
    -------
    DICTIONARY
        Dict containing the KEGG pathway ID, the path name, annotated proteins,
        the uniprot accession numbers of the proteins, and the kegg query names.

    """
    

    # import required modules
    import io
    import os
    import re
    import urllib
    import json
    import pandas as pd
    from Bio.KEGG import REST
    
    # Some code to return a Pandas dataframe, given tabular text
    def to_df(result):
        return pd.read_table(io.StringIO(result), header=None)
    
    print('\nInitializing KEGG tool...\n')
    
    # Define directory and filename to save the results to
    
    # loop to check if the folder already exists in that path and to make the folder if it does not exist yet
    if os.path.isdir(folder_spec):
        print ("\nDirectory %s already exists." %folder_spec)
    else: 
        try:
            os.mkdir(folder_spec)
        except OSError:
            print ("\nCreation of the directory %s failed" % folder_spec)
        else:
            print ("\nSuccessfully created the directory %s " % folder_spec)
    
    
    filename = folder_spec + '/'+filename 
    
    #  Retrieve all KEGG pathways and the annotated gene nemaes
    
    #retrieve all the yeast pathways
    result = REST.kegg_list("pathway", organism).read()
    result = to_df(result)
    
    result.columns = ['path number', 'path name']
    all_pathways = []
    
    for i, row in result.iterrows():
        all_pathways.append([re.split(':', row['path number'])[1], re.split(' - ', row['path name'])[0]])
    
    all_pathways = pd.DataFrame(all_pathways, columns =['Path number', 'Path name'])
    
    
    #retrieve the annotates genes for each pathway, and the corresponding Uniprot accessions
    
    col = 'genes(PREFERRED)'  #Uniprot query
    
    #if the file already (partly) exists, import the file
    if os.path.isfile(filename):
        print ("File %s exists\n" %filename)
        with open(filename, "r") as fp: 
            all_pathway_data = json.load(fp)
    else:
        all_pathway_data = {}
        print("File %s does not exist yet, importing KEGG data may take a while (15 mins)...\n" %filename)
    
    #iterate over each pathway and find annotated gene names     
    for i, row in all_pathways.iterrows():
        progressBar("Importing all KEGG pathways",i, len(all_pathways) - 1)
        pathway= row['Path number']
        pathname = row['Path name']
        
        if pathway in all_pathway_data.keys():
            continue
        
        result = REST.kegg_get(pathway).read()
        df = to_df(result)
        df.columns = ['A']
        
        df = df[df.A.str.contains(';')].iloc[1:]
        
        temp = [re.search('    (.*);', x).group(1) for x in df['A'].tolist()]
        temp = [re.split(' ', x) for x in temp]
        
        query = []
        gene_name = []
        accession = []
        remove_indx = []
        
        for i, x in enumerate(temp):
            if len(list(filter(len, x) )) != 2:
                continue
            query.append(list(filter(len, x))[0])
            gene_name.append(list(filter(len, x))[1])
            
        if len(query) != 0:    
            for j, item in enumerate(query):
                #print('{:.2f}% of q completed'.format(j/len(query)*100))
                # define urls 
                url_left='https://www.uniprot.org/uniprot/?query='
                url_right='&sort=score&columns=id,' + col + ',&format=tab'
                
                #Retrieve Uniprot data
                url_full = url_left+item+url_right                                      # full URL to extract info from example: https://www.uniprot.org/uniprot/?query=YGL253W&sort=score&columns=id,genes(PREFERRED),&format=tab
                item2 = urllib.request.urlopen(url_full).read()                         # read url
                read_output = item2.decode("utf-8").split('\n')                         # decode the url to make string from bytes and split at enters ('\n')
                if len(read_output) < 2:
                    remove_indx.append(j)
                    continue
                split_output = read_output[1].split('\t')                               # split the accession number from the GO terms by splitting at tabs ('\t')
                accession.append(split_output[0])
            
            for index in sorted(remove_indx, reverse=True):
                del query[index]
                del gene_name[index]
                
            all_pathway_data[pathway] = {'Path name'    : pathname,
                                         'Query'        : query,
                                         'Gene name'    : gene_name,
                                         'Accession'    : accession}  
    
        # save and publish the Uniprot_data list to a .txt file, for later access
        with open(filename, "w") as fp:                                         
            json.dump(all_pathway_data, fp)  
    
    #import MG pathway 
    pathway = 'sce00000'
    pathname = 'Minimal Glycolysis'   
    
    query = []
    accession = []
    
    if pathway not in all_pathway_data.keys():
        gene_name = ["HXK2", "PGI1", "PFK2", "PFK1", "FBA1", "TPI1", "TDH3", "PGK1", "GPM1", "ENO2", "CDC19", "PDC1", "ADH1", "ADH3"]
        
        for name in gene_name:
            if name in all_pathway_data['sce00010']['Gene name']:
                idx = all_pathway_data['sce00010']['Gene name'].index(name)
                query.append(all_pathway_data['sce00010']['Query'][idx])
                accession.append(all_pathway_data['sce00010']['Accession'][idx])
        
        all_pathway_data[pathway] = {'Path name'    : pathname,
                                     'Query'        : query,
                                     'Gene name'    : gene_name,
                                     'Accession'    : accession}  
    
        # save and publish the Uniprot_data list to a .txt file, for later access
        with open(filename, "w") as fp:                                         
            json.dump(all_pathway_data, fp)  
    
    return all_pathway_data

#%% Import proteomics data

def import_prot_data(folder, label_scheme, TR, groups, to_excel = False, excel_name = 'alldata'):
    """
    

    Parameters
    ----------
    folder : STRING. 
        Name of directory where the output files from PEAKS are saved
    label_scheme : DICTIONARY
        labelling scheme used for the different TMT experiments. Should be adjusted
        for each new set of experiments, otherwise the code does not work
    TR : NESTED LIST OF STRINGS
        e.g. [['126', '127C', '128C', '129C', '130C'],
              ['127N', '128N', '129N', '130N', '131']].
    groups : LIST OF STRINGS
        e.g. ['Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5'].
    to_excel : BOOL, optional
        If you want to save the data in one large excel file with the TMT experiments as sheet 
        names, choose True. The default is False.
    excel_name : STRING, optional
        As what do you want to save the excel file? The default is 'alldata'.

    Returns
    -------
    Backup: DICTIONARY.
        Dict containing the TMT exp numbers as keys and the data as values.

    """
    
    #import required modules
    import pandas as pd
    import re
    import os
    
    if to_excel:
        
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer = pd.ExcelWriter(os.getcwd()+'/'+folder+'/' + excel_name + '.xlsx', 
                                engine='xlsxwriter')
    
    #specify directory; list files in folder
    filenames = os.listdir(folder)  
    # create list of .csv filenames from folder                                                    
    filenames = [ filename for filename in filenames if filename.endswith(".csv") ]      
    filenames.sort()
    
    backup = {}
    
    #iterate over the different PEAKS output files
    for filename in filenames: 
        df          = pd.read_csv(folder+'/'+filename)
        df          = df.sort_values(by='Protein ID')
        sheetname   = re.split('.csv',filename)[0]
        sheetname   = re.split('_', sheetname)[1]
        
        if to_excel:
            # save all the crude data to an XlsxWriter Excel object.
            df.to_excel(writer, sheet_name=sheetname, index=False)
        
        #TMT07 and TMT14 have one less datapoint, adjust for that
        if 'TMT07' in filename or 'TMT14' in filename:
            index   = 3
            a       = TR[0][:index] + TR[0][index+1 :]
            b       = TR[1][:index] + TR[1][index+1 :]
            TR_2    = [a, b]
        else:
            TR_2    = TR
        
        # manually calculate average intensities, as PEAKS does this differently
        intensities = [col for col in df.columns if 'Intensity' in col]
        
        for i, tr in enumerate(TR_2[0]):
            average = df[[col for col in intensities if (TR_2[0][i] in col or TR_2[1][i] in col) and 'Group' not in col]].mean(axis=1)
            df[[col for col in intensities if groups[i] in col]] = average.to_frame()
            
        #find indexes/columns to keep
        keep_cols = [col for col in df.columns 
                      if 'Intensity' in col or 'Ratio' in col
                      or ('Protein ID' in col) 
                      or ('Significance' in col) 
                      or ('-10lgP' in col) 
                      or ('Unique' in col)
                      or ('Accession' in col)
                      or ('Description' in col)]
        
        df          = df[keep_cols]
        drop_cols   = []
        
        #adjust the column names of the dataframe, include the TMT exp number
        for i, col in enumerate(keep_cols):
            if 'Significance' in col or ('-10lgP' in col) or ('Unique' in col):
                keep_cols[i] = col+ ' ' + sheetname
                continue
            elif 'Intensity Group' in col:
                group = re.search('y (.*)\(TMT', col).group(1)
                if 'no data' in label_scheme[sheetname][group]:
                    df = df.drop(col, axis =1)
                    drop_cols.append(i)
                else:
                    keep_cols[i] = re.split(' ',col)[0] + ' ' + label_scheme[sheetname][group] + ' avg'
                continue
            elif 'Ratio Group' in col:
                group = re.search('o (.*)\(TMT', col).group(1)
                if 'no data' in label_scheme[sheetname][group]:
                    df = df.drop(col, axis =1)
                    drop_cols.append(i)
                else:
                    keep_cols[i] = re.split(' ',col)[0] + ' ' + label_scheme[sheetname][group] + ' avg'
                continue
            
            
            #adjust the technical replicate names and find empty columns
            for j, tr in enumerate(TR_2[0]):
                if tr in col:
                    if df[col].loc[0] == '0' or df[col].loc[0] == '-':
                        df = df.drop(col, axis =1)
                        drop_cols.append(i)
                        break
                    group = groups[j]
                    keep_cols[i] = re.split(' ',col)[0] + ' ' + label_scheme[sheetname][group] + ' TR1'
                    break
            for j, tr in enumerate(TR_2[1]):
                if tr in col:
                    if df[col].loc[0] == '0' or df[col].loc[0] == '-':
                        df = df.drop(col, axis =1)
                        drop_cols.append(i)
                        break
                    group = groups[j]
                    keep_cols[i] = re.split(' ',col)[0] + ' ' + label_scheme[sheetname][group] + ' TR2'
                    break
        
        # remove the empty columnbs
        for index in sorted(drop_cols, reverse=True):
            del keep_cols[index]    
            
        # store the data in a dictionary with the TMT exp number as keys
        df.columns          = keep_cols
        backup[sheetname]   = df
    
    if to_excel:
        writer.save()
    
    return backup

#%% adjust ratios so that the bridging sample is '1' in every case

def adjust_ratios(df):
    """
    

    Parameters
    ----------
    df : DICTIONARY
        TMT experiments as keys, df's containing the data from 1 TMT experiment as values.

    Returns
    -------
    backup : DICTIONARY
        Dict containing the data from the TMT experiments just as in the input, except 
        the ratios are normalised to the bridging samples.

    """
    
    backup = df
    
    for key in backup.keys():
        
        # find the columns containing the bridging sample intensities
        col_bridge_TR1      = [col for col in backup[key].columns 
                               if 'Intensity bridge' in col and 'TR1' in col]
        col_bridge_TR2      = [col for col in backup[key].columns 
                               if 'Intensity bridge' in col and 'TR2' in col]
        col_bridge_avg      = [col for col in backup[key].columns 
                               if 'Intensity bridge' in col and 'avg' in col]
        
        # find all the columns containing intensities
        col_intensity_TR1   = [col for col in backup[key].columns 
                               if 'Intensity' in col and 'TR1' in col]
        col_intensity_TR2   = [col for col in backup[key].columns 
                               if 'Intensity' in col and 'TR2' in col]
        col_intensity_avg   = [col for col in backup[key].columns 
                               if 'Intensity' in col and 'avg' in col]
        
        # find the column names in which the ratios should be stored
        col_ratio_TR1       = [col for col in backup[key].columns 
                               if 'Ratio' in col and 'TR1' in col]
        col_ratio_TR2       = [col for col in backup[key].columns 
                               if 'Ratio' in col and 'TR2' in col]
        col_ratio_avg       = [col for col in backup[key].columns 
                               if 'Ratio' in col and 'avg' in col]
        
        # calculate the normalised ratios by dividing the intensities by the bridge
        ratios_TR1          = backup[key][col_intensity_TR1].div(backup[key][col_bridge_TR1[0]], axis = 0) 
        ratios_TR2          = backup[key][col_intensity_TR2].div(backup[key][col_bridge_TR2[0]], axis = 0)
        ratios_avg          = backup[key][col_intensity_avg].div(backup[key][col_bridge_avg[0]], axis = 0)
     
        # replace the old ratios with the new, normalised ratios
        backup[key][col_ratio_TR1] = ratios_TR1
        backup[key][col_ratio_TR2] = ratios_TR2
        backup[key][col_ratio_avg] = ratios_avg
        
    return backup

#%% concatenate all data in one large dataframe

def concatenate_data(df, overlap = False):
    """
    

    Parameters
    ----------
    df : DICTIONARY
        Dict containing all the data from the TMT experiments.
    overlap : BOOL, optional
        Select True if you only want to keep genes/proteins that were detected in all 
        of the TMT experiments. The default is False.

    Returns
    -------
    all_combined : DATAFRAME
        A large df containing all the data from all the TMT experiments.

    """
    
    #import required modules
    import pandas as pd
    
    #save only certain columns, and paste all the data in one large dataframe
    for key in df.keys():
        
        # set the protein ID as index, so proteins detected in multiple TMT experiments
        # can be placed in the correct row
        df[key] = df[key].set_index('Protein ID')
        
        # import all the columns of the first TMT experiment data
        if key == next(iter(df)):   
            all_combined = df[key]
            continue
        
        # only keep certain columns of the other TMT exp data, to prevent duplicate columns
        keep_cols = [col for col in df[key].columns if (key in col) or ('Intensity' in col) 
                     or ('Ratio' in col) or ('Protein ID' in col)]
        
        
        if not overlap:
            all_combined = pd.concat([all_combined, df[key][keep_cols]], axis=1)
            
            #make sure the accession is noted and is not 'nan'
            all_combined.loc[all_combined.Accession != all_combined.Accession, 'Accession'] = df[key].loc[all_combined.Accession != all_combined.Accession, 'Accession']
        
        if overlap:
            #select only proteins that are present in all measurements
            all_combined = pd.concat([all_combined, df[key][keep_cols]], axis=1, join="inner")

    #remove 'no data' columns
    all_combined = all_combined.drop([col for col in all_combined.columns if 'no data' in col], axis=1)      
    
    return all_combined

#%% separate data into their biological replicate data

def separate_data(df, experiments, TMTs, breps):
    """
    

    Parameters
    ----------
    df : DATAFRAME
        Dataframe containing all the experimental data you want to separate.
    experiments : LIST OF STRINGS
        e.g. ['MG_AN', 'MG_O2', 'WT_AN', 'WT_O2'].
    TMTs : NESTED LIST OF STRINGS
        e.g. [['TMT01', 'TMT02', 'TMT03'], ['TMT04', 'TMT05', 'TMT06', 'TMT07'], 
              ['TMT08', 'TMT09', 'TMT10'], ['TMT11', 'TMT12', 'TMT13', 'TMT14']].
    breps : NESTED LIST OF STRINGS
        e.g. [['16', '17', '18'], ['16', '17', '18'], 
              ['19', '20', '21'], ['28', '29', '30']].

    Returns
    -------
    sep_data : DICTIONARY
        Contains the biological replicate description as keys and the separated 
        dataframes as values.

    """
    
    #import required modules
    import pandas as pd
    
    sep_data = {}                                                              # create placeholder for the output
    
    basic_cols = [col for col in df.columns
                  if 'Protein Group' in col
                  or 'Accession' in col
                  or 'Gene name' in col]
    
    for i, exp in enumerate(experiments):
        sep_data[exp] = {}
        for j, num in enumerate(breps[i]):
            
            # separate data from one biological replicate
            keep_cols   = [col for col in df.columns
                           if num in col and exp in col]        
            exp_data    = df[keep_cols]

            
            # calculate the normalised total intensity 
            keep_cols   = [col for col in exp_data.columns if 'Intensity' in col and 'avg' in col]
            df2         = exp_data[keep_cols]
            TI          = df2.sum(axis = 0, skipna = True)
            
            #correct TI's to make sure it's always 2500000000, for bar plots
            factor      = 2500000000/TI
            
            for col in TI.to_frame().T.columns:   
                corrected_intensity = (exp_data[col]*factor[col]).to_frame()
                exp_data.loc[:, corrected_intensity.columns[0].replace('Intensity', 'Corrected Intensity')] = corrected_intensity
            
            # extract Unique peptides
            unique = df['#Unique ' + TMTs[i][j]].to_frame()
            
            # concatenate all the data in one large dataframe
            exp_data = pd.concat([df[basic_cols], unique, exp_data], axis =1)
            exp_data = exp_data.dropna(subset = ['Gene name'])
            
            # save the separated data for each biological replicate
            sep_data[exp][exp+'_'+num] = exp_data 
    
    return sep_data


#%% clustermap

def make_clustermap(df, folder, figname, figtitle, xaxis, gene_labels = 'Gene name sign', pv_thr = 0.05, savefig=True, show=True, adjust_ratio = False, barplot=False, **kwargs):
    """
    

    Parameters
    ----------
    df : DATAFRAME
        Dataframe containing the ratios you want to plot.
    folder : STRING
        Directory where you want to save the figure.
    figname : STRING
        Name as what you want to save the figure.
    xaxis : LIST of strings
        Growth phases at which the data was measured. E.g. ['ME', 'LE', 'ES', 'MS']
    gene_labels : STRING, optional
        Column name where the gene names are listed. The default is 'Gene name sign'.
    pv_thr : FLOAT, optional
        p-value threshold. The default is 0.05.
    savefig : BOOL, optional
        Do you want to save the figure to your computer? The default is True.
    show : BOOL, optional
        Do you want to show the plot in Python? The default is True.
    adjust_ratio : BOOL, optional
        If you want to normalise the ratio to the first phase of the xaxis, fill in True. 
        If that was already done, fill in False. The default is False.
    barplot : BOOL, optional
        If True, plot two bargraphs next to the heatmap, showing the absolute abundances 
        of the proteins. The default is False.
    **kwargs : DICTIONARY

    Returns
    -------
    clustergrid

    """   
    
    #import required modules
    import numpy as np
    import re
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.pyplot import gcf

                
    keep_cols = [col for col in df.columns if 'Ratio' in col]
    cluster_data = df[keep_cols]
    
    if adjust_ratio:
        col_ratio = [col for col in cluster_data.columns if xaxis[0] in col]
        cluster_data = cluster_data.div(cluster_data[col_ratio[0]], axis = 0)
        
    cluster_data = np.log2(cluster_data.astype(float))

    keep_cols = []
    new_xaxis = []
    if len(xaxis) == 9:
        for item in xaxis:
            for col in cluster_data.columns:
                item2 = re.split('_', item)
                if item2[0] in col and item2[1] in col:
                    keep_cols.append(col)
                    new_xaxis.append(item)
    else:
        for item in xaxis:
            for col in cluster_data.columns:
                if item in col:
                    keep_cols.append(col)
                    new_xaxis.append(item)
    
    xaxis = new_xaxis
    cluster_data = cluster_data[keep_cols]
    cluster_data.columns = xaxis
    cluster_data.index = [item.title() for item in df[gene_labels].tolist()]
    cluster_data = cluster_data.replace([np.inf, -np.inf], np.nan)
    
    if 'row_colors' in kwargs.keys():
        colorscheme = dict(zip(df['Pathway'].unique(), plt.cm.rainbow(np.linspace(0,1,len(df['Pathway'].unique())))))
        row_colors = df['Pathway'].map(colorscheme)
        row_colors.index = cluster_data.index.tolist()
        row_colors = row_colors.rename('Pathway')
        kwargs['row_colors'] = row_colors
    
    if 'row_cluster' in kwargs.keys():
        if kwargs['row_cluster'] == False:
            mask = cluster_data.isnull()
            kwargs['mask'] = mask
    
    # extract the indexes from the clustering analysis
    clustergrid = sns.clustermap(cluster_data, **kwargs)
    clustergrid.ax_heatmap.set_yticklabels(clustergrid.ax_heatmap.get_ymajorticklabels(), fontsize = 13)
    clustergrid.ax_heatmap.set_xticklabels(clustergrid.ax_heatmap.get_xmajorticklabels(), fontsize = 14)
    
    if barplot:
        rowcol_legendpos = 0.23
    else:
        rowcol_legendpos = 0.5
        
    if 'row_colors' in kwargs.keys():
        plt.setp(clustergrid.ax_row_colors.get_xticklabels(), fontsize = 13)
        for label in df['Pathway'].unique():
            clustergrid.ax_col_dendrogram.bar(0, 0, color=colorscheme[label], label=label, linewidth=0);
            clustergrid.ax_col_dendrogram.legend(title='$\\bf{Pathway}$', title_fontsize=14, loc="center", ncol=2, fontsize = 14, bbox_to_anchor=(rowcol_legendpos, -0.045), bbox_transform=gcf().transFigure)    


    clustergrid.ax_cbar.set_title('$Log_2$  FoldChange', fontsize = 13)
    clustergrid.ax_cbar.tick_params(axis='y', labelsize =12)
    
    if barplot:
        clustergrid.fig.subplots_adjust(right=0.4)
        ax = clustergrid.fig.add_axes([0.5, -0.0, 0.24, 0.94])
        
        keep_cols = [col for col in df.columns if 'Intensity' in col and 'avg' in col]
        avg_intensity = df[keep_cols].mean(axis = 1)
        ax.barh(np.arange(len(df)), avg_intensity, log=True, zorder = 3)
        ax.set_yticks(np.arange(len(df)))
        ax.invert_yaxis()  # labels read top-to-bottom
        ax.set_yticklabels(cluster_data.index.tolist(), fontsize=13)
        ax.tick_params(axis='x', labelsize=13)
        ax.set_xlabel('Average intensity ($Log_{10} scale$', fontsize = 14)
        ax.grid(True, which = 'major', linewidth=1, axis = 'x', alpha = 0.6, zorder = 0)
        ax.grid(True, which = 'minor', linewidth=0.5, axis = 'x', alpha =0.4, zorder = 0)
        ax.set_xlim(10000, 180000000)
        
        ax2 = clustergrid.fig.add_axes([0.79, -0.0, 0.24, 0.94])
        
        keep_cols = [col for col in df.columns if 'Intensity' in col and 'avg' in col]
        avg_intensity = df[keep_cols].mean(axis = 1)
        ax2.barh(np.arange(len(df)), avg_intensity, log=False, zorder = 3)
        ax2.set_yticks(np.arange(len(df)))
        ax2.invert_yaxis()  # labels read top-to-bottom
        ax2.set_yticklabels(cluster_data.index.tolist(), fontsize=13)
        ax2.tick_params(axis='x', labelsize=13)
        ax2.set_xlabel('Average intensity', fontsize = 14)
        ax2.grid(True, which = 'major', linewidth=1, axis = 'x', alpha = 0.6, zorder = 0)
        ax2.grid(True, which = 'minor', linewidth=0.5, axis = 'x', alpha =0.4, zorder =0)
        ax2.set_xlim(0, 180000000)
    
    if savefig:
        clustergrid.savefig(folder+ '/'+ figname)
    
    if not show:
        plt.close()
    
    return clustergrid 

#%% select pathways of interest
def select_pathway(df, pathway_genes, no_data):
    """
    Function that selects proteins in proteomics data that belong to a specific pathway

    Parameters
    ----------
    df : dataframe containing the data of 1 biological replicate (which is the
        average of one or multiple technical replicates)
    pathway_genes : list containing the uniprot accession numbers and corresponding
        gene names
    no_data: dictionary where you want to store the undidentified proteins/genes

    Returns
    -------
    Out_path : dataframe containing only the experimental data of a certain
        pathway

    """
    
    #import required modukes
    import pandas as pd
    import numpy as np
    
    Out_path    = pd.DataFrame()
    gene_list   = pathway_genes
    empty       = pd.DataFrame(df.iloc[0]).T.applymap(lambda x: np.nan)
    out         = []
    yes_data    = {}
    
    for gene in gene_list:
        if gene in list(df['Gene name']):
            if gene not in yes_data.keys():
                out = (pd.DataFrame(df.iloc[np.where(df['Gene name'].str.contains(pat = gene))[0][0]]).T)
                if out.shape[1] == 2:
                    out = out.iloc[0]
                Out_path = pd.concat([Out_path, out])
                yes_data[gene] = 1
        else:
            empty['Gene name'] = gene
            out = empty
            Out_path = pd.concat([Out_path, out])
            if gene not in no_data.keys():
                no_data[gene] = 1
            else:
                no_data[gene] += 1
    
    return Out_path, no_data


#%% find up and downregulated genes
def up_down_regulated(df, lfc='Log2 Fold Change', pv='Avg p-value', lfc_thr=(0.58, 0.58), pv_thr=(0.05, 0.05)):
    
    equal           = df.loc[(df[lfc] >= -lfc_thr[1]) & (df[lfc] <= lfc_thr[0]) & (df[pv] < pv_thr[1])] #equal
    upregulated     = df.loc[(df[lfc] >= lfc_thr[0]) & (df[pv] < pv_thr[0])]  # upregulated
    downregulated   = df.loc[(df[lfc] <= -lfc_thr[1]) & (df[pv] < pv_thr[1])]  # downregulated
    
    return equal, upregulated, downregulated

# %% create volcano-, relative intensity- & fold change plot for each biological replicate


__all__ = ['gene_exp', 'general']

class gene_exp:
    
    """
    This class is used for visualising gene or protein expression differences via a volcano plot
    """

    def __init__(self):
        pass

    def geneplot(d, geneid, lfc, lfc_thr2, pv_thr, genenames, gfont, pv, gstyle):
        
        #import reuired modules
        import matplotlib.pyplot as plt
        import sys
        
        lfc_thr = lfc_thr2
        if genenames is not None and genenames == "deg":
            for i in d.index:
                if (d.loc[d.index == i, lfc].iloc[0] >= lfc_thr[0] and d.loc[d.index == i, pv].iloc[0] < pv_thr[0]) or \
                    (d.loc[d.index == i, lfc].iloc[0] <= -lfc_thr[1] and d.loc[d.index == i, pv].iloc[0] < pv_thr[1]):
                    if gstyle==1:
                        plt.text(d.loc[d.index == i, lfc].iloc[0], d.loc[d.index == i, 'logpv_add_axy'].iloc[0], i,
                                      fontsize=gfont)
                    elif gstyle==2:
                        plt.annotate(i, xy=(d.loc[d.index == i, lfc].iloc[0], d.loc[d.index == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(8, 8), textcoords='offset points', size=8,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
        elif genenames is not None and type(genenames) is tuple:
            for i in d.index:
                if i in genenames:
                    if gstyle==1:
                        plt.text(d.loc[d.index == i, lfc].iloc[0], d.loc[d.index == i, 'logpv_add_axy'].iloc[0], i,
                                      fontsize=gfont)
                    elif gstyle==2:
                        plt.annotate(i, xy=(d.loc[d.index == i, lfc].iloc[0], d.loc[d.index == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=4,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
        elif genenames is not None and type(genenames) is dict:
            for i in d.index:
                if i in genenames:
                    if gstyle==1:
                        plt.text(d.loc[d.index == i, lfc].iloc[0], d.loc[d.index == i, 'logpv_add_axy'].iloc[0],
                                      genenames[i], fontsize=gfont)
                    elif gstyle == 2:
                        plt.annotate(genenames[i], xy=(d.loc[d.index == i, lfc].iloc[0], d.loc[d.index == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)

    def volcano(df, lfc='Log2 Fold Change', pv='Avg p-value', lfc_thr=(0.58, 0.58), pv_thr=(0.05, 0.05), lfc_thr2=(0.58, 0.58), color=("green", "grey", "red", "orange"),
                valpha=1, geneid=None, genenames=None, gfont=12, dim=(5, 5), r=300, ar=0, dotsize=11, markerdot="o",
                sign_line=True, gstyle=1, show=False, axtickfontsize=14,
                axtickfontname="Arial", axlabelfontsize=15, axlabelfontname="Arial", axxlabel=None,
                axylabel=None, xlm=None, ylm=None, plotlegend=True, legendpos='upper left',
                figname='volcano', legendanchor=(1.05,1), plottitle = True, title='Volcano plot',
                legendlabels=['significantly upregulated', 'below significance threshold', 'significantly downregulated', 'significant'], theme=None, plotMG = False, df2="dataframe", geneid2=None, genenames2=None, Out_path_MG = None, MG_genenames = None):
        
        #import required modules
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.colors import ListedColormap
        
        _x = r'$ log_{2}$(Fold difference)'
        _y = r'$ -log_{10}$(P-value)'
        color = color
        # check if dataframe contains any non-numeric character
        assert general.check_for_nonnumeric(df[lfc]) == 0, 'dataframe contains non-numeric values in lfc column'
        assert general.check_for_nonnumeric(df[pv]) == 0, 'dataframe contains non-numeric values in pv column'
        # this is important to check if color or logpv exists and drop them as if you run multiple times same command
        # it may update old instance of df
        df = df.drop(['color_add_axy', 'logpv_add_axy'], axis=1, errors='ignore')
        assert len(set(color)) == 4, 'unique color must be size of 3'
        df.loc[(df[lfc] >= -lfc_thr[1]) & (df[pv] < pv_thr[1]), 'color_add_axy'] = color[3] #equal
        df.loc[(df[lfc] >= lfc_thr[0]) & (df[pv] < pv_thr[0]), 'color_add_axy'] = color[0]  # upregulated
        df.loc[(df[lfc] <= -lfc_thr[1]) & (df[pv] < pv_thr[1]), 'color_add_axy'] = color[2]  # downregulated
        df['color_add_axy'] = df['color_add_axy'].fillna(color[1])  # intermediate
        df['logpv_add_axy'] = -(np.log10(df[pv]))
        
        # plot
        assign_values = {col: i for i, col in enumerate(color)}
        color_result_num = [assign_values[i] for i in df['color_add_axy']]
        
        
        new_color = [col for i, col in enumerate(color) if i in color_result_num]
        color = new_color
        legendlabels = [leg for i, leg in enumerate(legendlabels) if i in color_result_num]
        assign_values = {col: i for i, col in enumerate(color)}
        color_result_num = [assign_values[i] for i in df['color_add_axy']]

        if theme == 'dark':
            general.dark_bg()
        plt.subplots(figsize=dim)
        if plotlegend:
            s = plt.scatter(df[lfc], df['logpv_add_axy'], c=color_result_num, cmap=ListedColormap(color), alpha=valpha,
                            s=dotsize, marker=markerdot)
            plt.legend(handles=s.legend_elements()[0], labels=legendlabels, loc=legendpos, fontsize=14, bbox_to_anchor=legendanchor)
        else:
            plt.scatter(df[lfc], df['logpv_add_axy'], c=color_result_num, cmap=ListedColormap(color), alpha=valpha,
                        s=dotsize, marker=markerdot)
        if sign_line:
            plt.axhline(y=-np.log10(pv_thr[0]), linestyle='--', color='#7d7d7d', linewidth=1, alpha=0.2)
            plt.axvline(x=lfc_thr[0], linestyle='--', color='#7d7d7d', linewidth=1, alpha=0.2)
            plt.axvline(x=-lfc_thr[1], linestyle='--', color='#7d7d7d', linewidth=1, alpha=0.2)
        gene_exp.geneplot(df, geneid, lfc, lfc_thr2, pv_thr, genenames, gfont, pv, gstyle)
        
        if plottitle:
            plt.title(title, fontsize = 18, fontweight = 'bold')
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
            
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        general.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)
        general.get_figure(show, r, figname, theme)
        
        return df
        
class general:
    
    """
    This class is used for plotting
    """
    
    def __init__(self):
        pass

    rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                   '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                   '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                   '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')

    @staticmethod
    def get_figure(show, r, fig_name, theme):
        
        #import required modules
        import matplotlib.pyplot as plt
        
        if show:
            plt.show()
        else:
            plt.savefig(fig_name, bbox_inches='tight', dpi=r)
        if theme == 'dark':
            plt.style.use('default')
        plt.clf()
        plt.close()


    @staticmethod
    def axis_labels(x, y, axlabelfontsize="small", axlabelfontname=None):
        
        #import required modules
        import matplotlib.pyplot as plt
        
        plt.xlabel(x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.ylabel(y, fontsize=axlabelfontsize, fontname=axlabelfontname)

    @staticmethod
    def axis_ticks(xlm=None, ylm=None, axtickfontsize="small", axtickfontname=None, ar=None):
        
        #import required modules
        import matplotlib.pyplot as plt
        import numpy as np
        
        if xlm:
            plt.xlim(left=xlm[0], right=xlm[1])
            plt.xticks(np.arange(xlm[0], xlm[1], xlm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.xticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.yticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

    @staticmethod
    def depr_mes(func_name):
        print("This function is deprecated. Please use", func_name )
        print("Read docs at https://reneshbedre.github.io/blog/howtoinstall.html")

    @staticmethod
    def check_for_nonnumeric(pd_series=None):
        
        #import required modules
        import pandas as pd
        
        if pd.to_numeric(pd_series, errors='coerce').isna().sum() == 0:
            return 0
        else:
            return 1

    @staticmethod
    def pvalue_symbol(pv=None, symbol=None):
        if 0.05 >= pv > 0.01:
            return symbol
        elif 0.01 >= pv > 0.001:
            return 2 * symbol
        elif pv <= 0.001:
            return 3 * symbol
        else:
            return None

    @staticmethod
    def get_file_from_gd(url=None):
        
        #import required modules
        import pandas as pd
        
        get_path = 'https://drive.google.com/uc?export=download&id=' + url.split('/')[-2]
        return pd.read_csv(get_path, comment='#')

    @staticmethod
    def dark_bg():
        
        #import required modules
        import matplotlib.pyplot as plt
        
        plt.style.use('dark_background')


#%% STRING background gene generator

def backgroundgene_2_string(gene_list):
    """
    

    Parameters
    ----------
    gene_list : LIST of strings
        List of background genes.

    Returns
    -------
    all_string_identifier : LIST
        Converted genes to string identifiers.

    """
    
    #import required modules
    import requests
    
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"
    
    # Set parameters
    params = {
    
        "identifiers"   : "\r".join(gene_list), # your protein list
        "species"       : 4932,                 # species NCBI identifier 
        "limit"         : 1,                    # only one (best) identifier per input protein
        "echo_query"    : 1,                    # see your input identifiers in the output
    }
    
    # Construct URL
    request_url = "/".join([string_api_url, output_format, method])
    
    # Call STRING
    results = requests.post(request_url, data=params)
    
    # Read and parse the result
    all_string_identifier = []
    for line in results.text.strip().split("\n"):
        l = line.split("\t")
        string_identifier = l[2]
        all_string_identifier.append(string_identifier)
    
    return all_string_identifier
    
#%% STRING

def STRING(changes, xaxis, print1, print2, filename, sp = 4932, categories = ['Function', 'Process','Component', 'KEGG'], fdr_thr = 0.05, output_format = "json", method = "enrichment", **kwargs):
    """
    

    Parameters
    ----------
    changes : DICTIONARY
        Dictionary containing the up- and downregulated genes for each growt phase.
    xaxis : LIST
        Growth phases.
    print1 : STRING
        First experiment you are comparing, eg 'MG AN'.
    print2 : STRING
        Second experiment you are comparing, eg 'WT AN'.
    filename : STRING
        Specify the name of the output csv file.
    sp : INT, optional
        Species number. The default is 4932, for S. cerevisae.
    categories : LIST of strings, optional
        The functional catagories you want to compare. The default is ['Function', 'Process','Component', 'KEGG'].
    fdr_thr : FLOAT, optional
        False discovery rate threshold. The default is 0.05.
    output_format : STRING, optional
        Desired output format from the STRING website. The default is "json".
    **kwargs : DICTIONARY
        Provide a list of background genes for example, use string identifiers for this.

    Returns
    -------
    save_data : DICTIONARY
        Save the output in the form of dataframes, to a dictionary with the growth phases as keys.

    """
    
    #import required modules
    import requests 
    import json
    import pandas as pd

    request_url = "/".join(["https://string-db.org/api", output_format, "enrichment"])
    
    save_data = {}
    
    for x in xaxis:
        if x not in changes:
            print("\nPhase %s not found in the dictionary" %x)
            continue
        if len(changes[x]['upregulated']) == 0:
            print('\nNo upregulated genes found in growth phase %s, skipping GO enrichment analysis' %x)
        else: 
             
            my_genes = changes[x]['upregulated']
            
            if "background" in kwargs.keys():
                params = {"identifiers" : "%0d".join(my_genes), # your genes
                          "species" : sp, # species NCBI identifier 
                          "background_string_identifiers" : "%0d".join(kwargs["background"])}
            else:
                params = {"identifiers" : "%0d".join(my_genes), # your genes
                          "species" : sp, # species NCBI identifier 
                          }
            
            response = requests.post(request_url, data=params)
            data = json.loads(response.text)
            
            df_cat, df_term, df_fdr, df_pvalue, df_genes, df_backgenes, df_description = [], [], [], [], [], [], []
            
            if len(data) == 0:
                print("\nNo enriched GO-terms found for %s phase." %x)
            
            else:
                print('\nEnriched GO terms in %s %s phase, compared to the %s are:' % (print1, x, print2))
                print("\t".join(["category","GO-term", "    fdr", "p-value", "#genes", '#genes_background' "    description"]))
                
                for row in data:
                    if 'error' in row.keys():
                        print(params["identifiers"])
                        continue
                    term = row["term"]
                    fdr = float(row["fdr"])
                    p_value = float(row["p_value"])
                    description = row["description"]
                    category = row["category"]
                    number_of_genes = int(row["number_of_genes"])
                    number_of_genes_in_background = int(row["number_of_genes_in_background"])
                
                    for cat in categories:
                        
                        if category == cat and fdr < fdr_thr: #'PMID','Component','Process','Function','KEGG','Interpro'
                            ##print significant GO Process annotations
                            df_cat.append(category)
                            df_term.append(term)
                            df_fdr.append(fdr)
                            df_pvalue.append(p_value)
                            df_genes.append(number_of_genes)
                            df_backgenes.append(number_of_genes_in_background)
                            df_description.append(description)
                            
                            ##print significant GO Process annotations
                            if category == "Component" or category == 'Function':
                                fdr = str('\t'+str(fdr))
                                description = str('\t' + description)
                                print("*"+("\t".join([category,term, str(fdr), str(p_value), str(number_of_genes), str(number_of_genes_in_background), description])))
                            elif category == "KEGG":
                                fdr = str('\t'+str(fdr))
                                print("*"+("\t".join([category,term, str(fdr), str(p_value), str(number_of_genes), str(number_of_genes_in_background), description])))
                            else:
                                print('*'+("\t".join([category,term, str(fdr), str(p_value), str(number_of_genes), str(number_of_genes_in_background), description]))) #for printing everything
                
                df = pd.DataFrame()
                df['Category'] = df_cat
                df['GO term'] = df_term
                df['FDR'] = df_fdr
                df['p-value'] = df_pvalue
                df['#genes'] = df_genes
                df['#genes background'] = df_backgenes
                df['Description'] = df_description
                
                df.to_csv(filename+'_enriched_'+x+'.csv', index=False)
                
                save_data[x] = {'upregulated': df_description}
                
        if len(changes[x]['downregulated']) == 0:
            print('\nNo downregulated genes found in growth phase %s, skipping GO enrichment analysis' %x)
        else: 
            my_genes = changes[x]['downregulated']
            
            if "background" in kwargs.keys():
                params = {"identifiers" : "%0d".join(my_genes), # your genes
                          "species" : sp, # species NCBI identifier 
                          "background_string_identifiers" : "%0d".join(kwargs["background"])}
            else:
                params = {"identifiers" : "%0d".join(my_genes), # your genes
                          "species" : sp, # species NCBI identifier 
                          }
            
            response = requests.post(request_url, data=params)
            data = json.loads(response.text)
            
            df_cat, df_term, df_fdr, df_pvalue, df_genes, df_backgenes, df_description = [], [], [], [], [], [], []
        
            
            if len(data) == 0:
                print("\nNo depleted GO-terms found for %s phase." %x)
            
            else:

                print('\nDepleted GO terms in %s %s phase, compared to the %s are:' % (print1, x, print2))
                print("\t".join(["category","GO-term", "    fdr", "p-value", "#genes", '#genes_background' "    description"]))
    
                for row in data:
                    if 'error' in row.keys():
                        print(params["identifiers"])
                        continue
                    term = row["term"]
                    fdr = float(row["fdr"])
                    p_value = float(row["p_value"])
                    description = row["description"]
                    category = row["category"]
                    number_of_genes = int(row["number_of_genes"])
                    number_of_genes_in_background = int(row["number_of_genes_in_background"])
                
                    for cat in categories:
                        
                        if category == cat and fdr < fdr_thr: #'PMID','Component','Process','Function','KEGG','Interpro'
                            df_cat.append(category)
                            df_term.append(term)
                            df_fdr.append(fdr)
                            df_pvalue.append(p_value)
                            df_genes.append(number_of_genes)
                            df_backgenes.append(number_of_genes_in_background)
                            df_description.append(description)
                            ##print significant GO Process annotations
                            if category == "Component" or category == 'Function':
                                fdr = str('\t'+str(fdr))
                                description = str('\t' + description)
                                print("*"+("\t".join([category,term, str(fdr), str(p_value), str(number_of_genes), str(number_of_genes_in_background), description])))
                            elif category == "KEGG":
                                fdr = str('\t'+str(fdr))
                                print("*"+("\t".join([category,term, str(fdr), str(p_value), str(number_of_genes), str(number_of_genes_in_background), description])))
                            else:
                                print('*'+("\t".join([category,term, str(fdr), str(p_value), str(number_of_genes), str(number_of_genes_in_background), description]))) #for printing everything
                
                if len(df_cat) != 0:
                
                    df = pd.DataFrame()
                    df['Category'] = df_cat
                    df['GO term'] = df_term
                    df['FDR'] = df_fdr
                    df['p-value'] = df_pvalue
                    df['#genes'] = df_genes
                    df['#genes background'] = df_backgenes
                    df['Description'] = df_description
                    
                    
                    df.to_csv(filename+'_depleted_'+x+'.csv', index=False)
                    
                    save_data[x] = {'downregulated': df_description}
                    
    return save_data