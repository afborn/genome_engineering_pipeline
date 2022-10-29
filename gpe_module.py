#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# gpe module
#import for retrieving predicted gRNAs
from pybiomart import Server
import pandas as pd
import urllib.request
import os
# import nest_asyncio
# nest_asyncio.apply()
from selenium import webdriver
import time
# pd.set_option('display.max_rows', None)
# pd.options.display.max_columns = None
# pd.options.display.max_colwidth = 200
#imports for primer design
from Bio.Seq import Seq
import primer3
import numpy as np
import requests
import sys
import warnings
warnings.filterwarnings('ignore')

# define functions

def check_csv_file(file):    
    df_from_file = pd.read_csv(file)    
    try:
        df_from_file.columns = df_from_file.columns.str.upper()
    except:
        raise SystemExit("Problem with csv file")        
    try:
        if 'GENE ID' in df_from_file.columns:
            df_from_file['GENE ID'] = df_from_file['GENE ID'].str.upper()
    except:
        raise SystemExit("Please provide a csv file containing column named 'GENE ID'")
    try:
        number_of_ensg_entries = df_from_file['GENE ID'].str.contains('ENSG').sum()
        if not len(df_from_file) == (number_of_ensg_entries):
            raise SystemExit("Please check that supplied csv file contains a column named 'GENE ID'")
    except:
        raise SystemExit("Program aborted - please retry")
    
    print("Successfully processed csv file")
    
    return df_from_file


def extract_gene_name_from_id(csv_file_ens_gene_id, dataset='hsapiens_gene_ensembl'):
    '''Function to extract gene names using gene ids.
    Function accepts a dataframe that is supplied by user and that 
    contains ENSEMBL gene ids in a column labeled "Gene id".
    Function will test for existence of column "Gene id", prior to extracting gene names 
    (and gene ids) using Biomart Ensembl server.
    Function returns a dataframe containing ENSEMBL gene id's and gene names.
    '''
    #initiate dict to store gene sequences for ENSEMBL names
    gene_id = []
    gene_symbol = []
    
    #check if supplied csv correct as df 
    df_ens_gene_name = check_csv_file(csv_file_ens_gene_id)   
    
    try:
        #connect to Biomart Server
        server = Server(host='http://www.ensembl.org')
        # use #dataset.list_filters() to see available filters for list_filters method
        #dataset.list_filters()
        #generate dataset (homo sapiens)
        dataset = (server.marts['ENSEMBL_MART_ENSEMBL']
                   .datasets[dataset])
        
        for entry in range(len(df_ens_gene_name)):
            gene_ens_id = df_ens_gene_name["GENE ID"].iloc[entry]
            #parse dataset for gene_id
            data_gene = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                          filters={'link_ensembl_gene_id': gene_ens_id})
            
            #if gene_id found and unique, return gene_id, gene_symbol
            if len(data_gene.index) == 1:
                gene_id.append(data_gene.iat[0,0])
                gene_symbol.append(data_gene.iat[0,1])
    except:
        raise SystemExit("Program aborted - unable to extract human gene names for supplied gene Id's")
        #abort
    
    print("Successfully identified gene names from IDs")
    
    return pd.DataFrame(zip(gene_id, gene_symbol), columns=('GENE ID', 'GENE NAME'))



def extract_gene_seq_from_ens_id(csv_file_ens_gene_id):        
    #initiate dict to store gene sequences for ENSEMBL IDs
    gene_sequences = []
    gene_id = []
        
    #check if supplied csv file is ok
    df_ens_gene_id = check_csv_file(csv_file_ens_gene_id)
    
    try:
        #REST API python3 Ensembl
        server = "https://rest.ensembl.org"
        #extract unique gene IDs in case of redundant ENSG ID entries in df
        unique_gene_IDs = df_ens_gene_id['GENE ID'].unique()
            #loop through unique_gene_IDs and extract gene sequence from ensembl
        for entry in unique_gene_IDs:
            ext = "/sequence/id/" + entry +"?"
            # retrieve plain text gene sequence
            r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
            # sanity check object r
            if not r.ok:
                r.raise_for_status()
                sys.exit("Unable to retrieve gene sequences from Ensembl - please try again later")
            gene_sequences.append(r.text)
            gene_id.append(entry)
    except:
        raise SystemExit("Unable to retrieve gene sequences from Ensembl")
    
    print("Successfully retrieved gene sequences from IDs")
    
    return pd.DataFrame(zip(gene_id, gene_sequences), columns=('GENE ID', 'GENE SEQUENCES'))



def construct_synthego_urls_gene_names_ids(df_ens_genes):
    '''Function constructs & returns a dataframe (df) of Synthego website URLs used to predict gRNAs.
    Function requires a df, comprising two columns labeled "GENE ID" and "GENE NAME"
    Function returns a dataframe containing a URL for all genes specified in input df, as well as 
    related gene id and gene name.
    '''
    try:
        #Capital letter column header
        df_ens_genes.columns = df_ens_genes.columns.str.upper()
        df_ens_genes["GENE ID"] = df_ens_genes["GENE ID"].str.upper()
        df_ens_genes["GENE NAME"] = df_ens_genes["GENE NAME"].str.upper()
    except:
        raise SystemExit("Program aborted - dataframe does no contain correct column labels'")
        #abort
    
    urls = []
    gene_id = []
    gene_name = []
    
    for entry in range(len(df_ens_genes)):
        url = 'https://design.synthego.com/#/design/results?genome=homo_sapiens_gencode_26_primary&nuclease=cas9&gene_id='         +df_ens_genes["GENE ID"].iloc[entry]+'&symbol='+df_ens_genes["GENE NAME"].iloc[entry]
        urls.append(url)
        gene_id.append(df_ens_genes["GENE ID"].iloc[entry])
        gene_name.append(df_ens_genes["GENE NAME"].iloc[entry])
            
    df_url = pd.DataFrame(zip(gene_id, gene_name, urls), columns=('GENE ID', 'GENE NAME', 'SYNTHEGO URL'))
    
    print("Successfully constructed URLs from gene names and gene id's for scraping predicted gRNAs")
    
    return df_url



def predict_gRNA_from_urls(df_genes_urls, path_to_gecko='D:\geckodriver\geckodriver.exe'):
    '''Function requires a dataframe containing ENSEMBL gene ids, gene names and Synthego URLs used for gRNA predicitons
    for respective genes. Required column labels = "GENE ID", "GENE NAME", "SYNTHEGO URL".
    Function returns a dataframe with ENSEMBL gene id's, gene names and 4 predicted gRNAs per gene.
    Using this function requires geckodriver and Firefox installed. 
    Function takes two arguments - 1. a dataframe that contains columns "GENE ID", "GENE NAME" and "SYNTHEGO URL", 2. absolute path to installed geckodriver (default value set)
    '''
    path_to_gecko = path_to_gecko
    counter = len(df_genes_urls)
    print("Retrieving gRNAs for " + str(counter) + " gene(s) \nWorking...")
   
    try:
        #Capital letter column header
        df_genes_urls.columns = df_genes_urls.columns.str.upper()
    except:
        raise SystemExit("Program aborted - dataframe does no contain correct column labels'")
        #abort
    
    #Store predictions and related gene ID's and gene names
    gRNAs_synth_array = []
    gene_id_array = []
    gene_name_array = []
    
    try:
        for entry in range(len(df_genes_urls)):
            url = df_genes_urls["SYNTHEGO URL"].iloc[entry]
            #suppress opening Firefox browser window
            os.environ['MOZ_HEADLESS'] = '1'
            # run firefox webdriver from executable path of your choice
            driver = webdriver.Firefox(executable_path = path_to_gecko)
            # get web page
            driver.get(url)
            #30 seconds wait to have page fully loaded - might need to increase time
            time.sleep(30)
            results = driver.find_elements("xpath", "//*[@class='ng-binding']")
            #extract gene symbol for respective entry
            gene_symbol = df_genes_urls["GENE NAME"].iloc[entry]
            #extract gene id for respective entry
            gene_id = df_genes_urls["GENE ID"].iloc[entry]
           
            gRNAs_synth = []
            gRNAs_synth_min = []
 
            # loop over results
            for result in results:
                gRNAs_synth.append(result.text)
                #extract minimal information for predicted gRNAs
                gRNAs_synth_min = gRNAs_synth[5:9]
            gRNAs_synth_array.append(gRNAs_synth_min)           
            
            #flatten array gRNAs_synth_array
            gRNAs_synth_array_final = [item for sublist in gRNAs_synth_array for item in sublist]
            
            #add gene id and gene name to arrays for later zipping of df
            for entry in range(4):
                gene_id_array.append(gene_id)
                gene_name_array.append(gene_symbol)
                
            counter -= 1
            if counter > 1:
                print("Remaining gRNA predictions: " + str(counter) + " genes \nWorking...")
            elif counter == 1:
                print("Retrieving gRNA predictions for last gene \n...almost there...")
                
    except:
        raise SystemExit("Failed to retrieve predicted gRNAs!") 
    
    #put lists together into df
    try:
        df_predictions = pd.DataFrame(zip(gene_id_array, gene_name_array, gRNAs_synth_array_final), columns=('GENE ID', 'GENE NAME', 'PREDICTED GRNA'))
    except:
        raise SystemExit("Unable to construct dataframe with gRNA predictions!")
    
    print("Predictions succesfully completed!")
    
    return df_predictions


def gRNA_hybridisation(df_gRNA_seq):
    '''Function to check if gRNA's sequence is identical to DNA coding strand, or non-coding-strand.
    
    In case, gRNA is identical in sequence to non-coding strand, the gRNA sequence will be reverse-complemented 
    to turn it into its coding-strand counterpart.
    
    The function takes 1 argument:
    1) Requirement: a dataframe containing columns, "GENE NAME", "PREDICTED GRNA" (containing gRNA sequences), "GENE ID" 
    (containing ENSEMBL Gene ID's that are predicted to be targeted by gRNA) and related DNA sequences ("GENE SEQUENCES)".
    Gene ID and Gene name according to ENSEMBL definitions.
    
    Function returns a dataframe containing the columns "GENE ID", "PREDICTED GRNA", "PREDICTED GRNA_T" & 
    "GRNA CODING STRAND" (sequence of gRNA on DNA coding strand), "GRNA LOCATION" (location of gRNA on DNA coding strand), 
    "GRNA HYBRIDISATION" (indicating which strand the "PREDICTED GRNA" is binding to - values either "coding strand" 
    or "non-coding strand")
    '''  
    
    try:
        #Capital letter column header
        df_gRNA_seq.columns = df_gRNA_seq.columns.str.upper()
        df_gRNA_seq['GENE ID'] = df_gRNA_seq['GENE ID'].str.upper()
        df_gRNA_seq['PREDICTED GRNA'] = df_gRNA_seq['PREDICTED GRNA'].str.upper()
        #replace U with T in gRNA sequences
        df_gRNA_seq["PREDICTED GRNA_T"] = df_gRNA_seq["PREDICTED GRNA"].str.replace("U",'T')  
    except:
        raise SystemExit("Dataframe does not contain correct column labels or column entries!")
        
    try:
        #array for storing orientation result
        for_rev_orientation_gRNA = []

        #using .find to match GRNA_T to gene sequence - if not found output will be -1.
        for row in range(len(df_gRNA_seq)):
            #get key from gRNA_predicted to retrieve value from gene_sequences dict
            #key = df_gRNA_seq["GENE ID"].iloc[row]
            value = df_gRNA_seq["GENE SEQUENCES"].iloc[row]
            orientation = value.find(df_gRNA_seq["PREDICTED GRNA_T"].iloc[row])
            for_rev_orientation_gRNA.append(orientation)

        #add for_rev_orientation_gRNA to df_gRNAs
        df_gRNA_seq["for_rev_gRNA"] = for_rev_orientation_gRNA

        #create column and pre-populate
        df_gRNA_seq["GRNA CODING STRAND"] = df_gRNA_seq["PREDICTED GRNA_T"]

        #In case of output -1 in for_rev_orientation_gRNA use reverse complement gRNA to search for match
        for entry in range(len(df_gRNA_seq)):
            #check if column "for_rev_gRNA" contains -1, if so reverse complement the PREDICTED GRNA_T sequence and
            #store it in GRNA ON CODING STRAND
            if df_gRNA_seq.iloc[entry, 5] == -1:
                df_gRNA_seq["GRNA CODING STRAND"].iloc[entry] = Seq(df_gRNA_seq["PREDICTED GRNA_T"].iloc[entry]).reverse_complement()
                df_gRNA_seq["GRNA CODING STRAND"].iloc[entry] = str(df_gRNA_seq["GRNA CODING STRAND"].iloc[entry])

        #determine gRNA position on coding strand
        forward_position = []
        #determine location of gRNA on coding strand
        for row in range(len(df_gRNA_seq)):
            #get key from gRNA_predicted to retrieve value from gene_sequences dict
            #key = df_gRNA_seq["GENE ID"].iloc[row]
            value = df_gRNA_seq["GENE SEQUENCES"].iloc[row]
            orientation = value.find(df_gRNA_seq["GRNA CODING STRAND"].iloc[row])
            forward_position.append(orientation)

        #add for_rev_orientation_gRNA to df_gRNAs
        df_gRNA_seq["GRNA LOCATION"] = forward_position    

        #create new column with empty strings
        df_gRNA_seq["GRNA HYBRIDISATION"] = ""

        #In case of output -1 in for_rev_gRNA, enter value 'to coding strand' to entry, otherwise 'to non-coding strand'
        for entry in range(len(df_gRNA_seq)):
            if df_gRNA_seq.iloc[entry, 5] == -1:
                df_gRNA_seq["GRNA HYBRIDISATION"].iloc[entry] = "to coding strand"
            else:
                df_gRNA_seq["GRNA HYBRIDISATION"].iloc[entry] = "to non-coding strand"

        df_gRNA_orientation = df_gRNA_seq[["GENE ID", "GENE NAME", "PREDICTED GRNA", "GRNA LOCATION", "GRNA HYBRIDISATION", "GRNA CODING STRAND",                                          "PREDICTED GRNA_T", "GENE SEQUENCES"]]
    except:
        raise SystemExit("Unable to create dataframe with information on gRNAs!")
    
    print("Starting gRNA analyses...")
    time.sleep(2)
    print("Analyses of gRNA completed successfully!")  
    
    return df_gRNA_orientation


def add_count_to_gRNA(df):
    df["GRNA NAME"] = ""
    gRNA_list = ['_gRNA1', '_gRNA2', '_gRNA3', '_gRNA4']
    tempList = list(gRNA_list)
    count = int((len(df)/4)-1)
    try:
        if (len(df)%4) == 0:
            for i in range(count):
                for element in tempList:
                    gRNA_list.append(element)
            df["REPEAT"]= gRNA_list
            df["GRNA NAME"] = df["GENE NAME"] + df["REPEAT"]
            df = df.drop(["REPEAT"], axis=1)
    except:
        raise SystemExit("Dataframe does not contain 4 gRNAs per gene")
    
    return df
        

def primer3_primer_around_gRNA(df_gRNAs_info, upstream_dist=250, downstream_dist=250):

    #list to hold sequences
    sequences_for_primer_design = []

    for row in range(len(df_gRNAs_info)):
        beginning = df_gRNAs_info["GRNA LOCATION"].iloc[row] - 250
        end = df_gRNAs_info["GRNA LOCATION"].iloc[row] + 250
        key = df_gRNAs_info["GRNA NAME"].iloc[row]
        value = df_gRNAs_info["GENE SEQUENCES"].iloc[row]
        seq_slice = value[beginning:end]
        sequences_for_primer_design.append(seq_slice)

    df_gRNAs_info['SEQ SLICE'] = sequences_for_primer_design


    primer_designed = {}
    for entry in range(len(df_gRNAs_info)):
        seq_dict = {
            'SEQUENCE_ID': df_gRNAs_info['PREDICTED GRNA'].iloc[entry],
            'SEQUENCE_TEMPLATE': df_gRNAs_info['SEQ SLICE'].iloc[entry],
        }
        primer_designed[df_gRNAs_info['GRNA NAME'].iloc[entry]] = primer3.designPrimers(seq_dict,    
            {
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_PICK_INTERNAL_OLIGO': 1,
                'PRIMER_INTERNAL_MAX_SELF_END': 8,
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 25,
                'PRIMER_OPT_TM': 60.0,
                'PRIMER_MIN_TM': 57.0,
                'PRIMER_MAX_TM': 63.0,
                'PRIMER_MIN_GC': 20.0,
                'PRIMER_MAX_GC': 80.0,
                'PRIMER_MAX_POLY_X': 100,
                'PRIMER_INTERNAL_MAX_POLY_X': 100,
                'PRIMER_SALT_MONOVALENT': 50.0,
                'PRIMER_DNA_CONC': 50.0,
                'PRIMER_MAX_NS_ACCEPTED': 0,
                'PRIMER_MAX_SELF_ANY': 12,
                'PRIMER_MAX_SELF_END': 8,
                'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                'PRIMER_PAIR_MAX_COMPL_END': 8,
                'PRIMER_PRODUCT_SIZE_RANGE': [[375, 500]],
            })

    primer_df = pd.DataFrame.from_dict(primer_designed, orient='index')
    primer_df["GRNA NAME"] = primer_df.index
    primer_df = primer_df[["GRNA NAME", "PRIMER_LEFT_0_SEQUENCE", "PRIMER_RIGHT_0_SEQUENCE", "PRIMER_PAIR_0_PRODUCT_SIZE", "PRIMER_LEFT_0", "PRIMER_RIGHT_0", "PRIMER_INTERNAL_0", "PRIMER_LEFT_0_TM", "PRIMER_RIGHT_0_TM"]]

    return primer_df

# gRNA_primer_df = gRNA_predicted.join(primer_df)
# gRNA_primer_df = gRNA_primer_df.drop(columns=["gRNA name", "gRNA reverse (-1 = True)"])
# gRNA_primer_df_wo_seq = gRNA_primer_df.drop(columns=["Seq for primer design"])
# gRNA_primer_df_wo_seq.head(10)

# gRNA_primer_df_wo_seq.to_csv(project_id + '_gRNA_primer.csv')

# #check if all genes have predicted gRNAs - if not, run predicted_gRNA function until all genes have predicted gRNAs)
# while len(gRNA_predicted[gRNA_predicted['gRNA sequence'].isnull()]) > 0:
#     print(gRNA_predicted[gRNA_predicted['gRNA sequence'].isnull()])
#     nan_rows = gRNA_predicted[gRNA_predicted['gRNA sequence'].isnull()]
#     nan_rows = nan_rows.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index='False')
#     nan_rows = nan_rows[['GENE ID', 'GENE NAME']]
#     nan_rows.columns = nan_rows.columns.str.upper()
#     predict_gRNA(nan_rows)


if __name__ == "__main__":
    print("This is just a module")

