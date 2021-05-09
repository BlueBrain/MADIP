#!/usr/bin/env python
# coding: utf-8

"""
MADIP: Molecular Atlas Data Integration Pipeline

This module provide helper functions for step_2_protein_ids_alignment.ipynb jupyter notebook


Copyright 2021 Blue Brain Project / EPFL 

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   
"""

import pandas as pd
import numpy as np

import re

import pickle as pkl

def process_uniprot_mapping_data():
    """
    Return two pandas df and two dicts with 
    ['Entry', 'Entry name', 'Status', 'Protein names', 'Gene names','Organism', 'Length', 'gene_id_entry_name'] 
    based on manually obtained Uniprot data from https://www.uniprot.org 
    with the queries as
    (taxonomy:"Mus musculus (Mouse) [10090]" OR taxonomy:"Rattus norvegicus (Rat) [10116]" OR taxonomy:"Homo sapiens (Human) [9606]") AND reviewed:yes
    (taxonomy:"Mus musculus (Mouse) [10090]" OR taxonomy:"Rattus norvegicus (Rat) [10116]" OR taxonomy:"Homo sapiens (Human) [9606]") AND reviewed:no
    """
    uniprot_rev = pd.read_csv('../data/uniprot_rev_taxonomyMRH_21july2020.tab', sep='\t')
    print("uniprot_rev: ",len(uniprot_rev))

    uniprot_unrev = pd.read_csv('../data/uniprot_taxonomyMRH_unreviewed_21july2020.gz', sep='\t')
    print("uniprot_unrev: ",len(uniprot_unrev))

    uniprot_rev['Gene names'] = uniprot_rev['Gene names'].str.upper()
    uniprot_rev['Gene names'] = uniprot_rev['Gene names'].str.split(' ')

    uniprot_unrev['Gene names'] = uniprot_unrev['Gene names'].str.upper()
    uniprot_unrev['Gene names'] = uniprot_unrev['Gene names'].str.split(' ')

    print("check == 1.0:",len(uniprot_rev['Entry'].unique())/len(uniprot_rev['Entry']))
    print("check == 1.0:",len(uniprot_unrev['Entry'].unique())/len(uniprot_unrev['Entry']))

    uniprot_rev = uniprot_rev.loc[~uniprot_rev['Gene names'].isna()]
    uniprot_unrev = uniprot_unrev.loc[~uniprot_unrev['Gene names'].isna()]

    uniprot_rev_dict = pd.Series(uniprot_rev['Gene names'].values,index=uniprot_rev['Entry']).to_dict() 
    uniprot_unrev_dict = pd.Series(uniprot_unrev['Gene names'].values,index=uniprot_unrev['Entry']).to_dict() 

    uniprot_rev_genes = list(set([item for sublist in uniprot_rev['Gene names'].tolist() for item in sublist]))
    uniprot_unrev_genes = list(set([item for sublist in uniprot_unrev['Gene names'].tolist() for item in sublist]))

    print("rev gene names",len(uniprot_rev_genes))
    print("unrev gene names",len(uniprot_unrev_genes))

    uniprot_rev['gene_id_entry_name'] = None
    for idx,row in uniprot_rev.iterrows():
        uniprot_rev.loc[idx, 'gene_id_entry_name'] = row['Entry name'].split('_')[0]

    uniprot_unrev['gene_id_entry_name'] = None
    for idx,row in uniprot_unrev.iterrows():
        uniprot_unrev.loc[idx, 'gene_id_entry_name'] = row['Entry name'].split('_')[0]

    print(len(uniprot_rev['gene_id_entry_name'].unique()))
    print(len(uniprot_unrev['gene_id_entry_name'].unique()))
    
    return uniprot_rev, uniprot_unrev, uniprot_rev_dict, uniprot_unrev_dict



def get_gene_unified(index, row):
    """
    Count number of gene name occurences, return candidate unified gene names
    """
    gn_list = row["gene_names"].replace(" ","").split(";")
    
    counts_occ = dict()
    for elem in gn_list:
        counts_occ[elem] = gn_study_dict.get(elem,0.0)
        
    # get keys with max value:
    max_value = max(counts_occ.values())
    
    if max_value >0:
        gns_mostFreq = [k for k,v in counts_occ.items() if v == max_value]
    else:
        gns_mostFreq = np.nan
        
    if isinstance(gns_mostFreq, list):
        return gns_mostFreq[0]
        
    else:
        return gn_list[0]
    
    
def get_uniprot_unified(index, row):
    """
    Count number of UniProt id occurences, return candidate unified UniProt ids
    """
    uni_list0 = row["Uniprot"].replace(" ","").split(";")
    uni_list = [x.split("-")[0] for x in uni_list0 if x is not None] 
    
    counts_occ = dict()
    for elem in uni_list:
        counts_occ[elem] = uniprot_study_dict.get(elem,0.0)
        
    # get keys with max value:
    max_value = max(counts_occ.values())
    
    if max_value >0:
        uni_mostFreq = [k for k,v in counts_occ.items() if v == max_value]
    else:
        uni_mostFreq = np.nan
        
    if isinstance(uni_mostFreq, list):
        return uni_mostFreq[0]
        
    else:
        return uni_list[0]
    
def check_GN_consistency_within_studies(df_all):
    """
    Check for GN consistency within studies. Get gene names repeated in different entries in combination with other gn. Example: UCHL3;UCHL4 + UCHL3
    """
    multiids = {}
    for i,study in enumerate(df_all['Study'].unique()):
        multiids[study] = []
        for gn in df_all.loc[df_all['Study']==study,'gene_names'].unique():
            if type(gn) == str:
                if ';' in gn:
                    gns = gn.split(';')
                    dupl = False
                    for ids in gns:
                        if ids in df_all.loc[df_all['Study']==study,'gene_names'].unique():
                            dupl = True
                    if dupl == True:
                        multiids[study].append(gn)
    return multiids



#check if at least one item in list exists in another list
any_in = lambda a, b: any(i in b for i in a)


def get_uniprot_raw_data(index,row,uniprots_list,uniprot_rev_values_for_uniprots_list,uniprot_unrev_values_for_uniprots_list,uniprot_ids_mrh_dict_values_for_uniprots_list):
    """
    Count occurences of ids to find the most frequent, consider id type and reliability (reviewed/unreviewed)
    """
    gene_name_unified = row['gene_name_unified']
    gene_names = row['gene_names']
    uniprot = row['Uniprot']
    
    #case 4
    
    #  before run gn alignment, pre-clean data from cases with no ids at all
    
    # 4.1: uniprot_rev_dict
    #uniprot_rev_values_for_uniprots_list =  list(set([item for sublist in [uniprot_rev_dict.get(i) for i in uniprots_list] for item in sublist])) # type: gene names
    # 4.2: uniprot_unrev_dict
    #uniprot_unrev_values_for_uniprots_list =  list(set([item for sublist in [uniprot_unrev_dict.get(i) for i in uniprots_list] for item in sublist]))  
    # 4.3: uniprot_ids_mrh_dict (based on df uniprot_gn)
    #uniprot_ids_mrh_dict_values_for_uniprots_list =  list(set([item for sublist in [uniprot_ids_mrh_dict.get(i) for i in uniprots_list] for item in sublist]))  
    
    
    
    
    # 4.1: uniprot_rev_dict
    if len(uniprot_rev_values_for_uniprots_list) >0:
        counts_occ = dict()
        for elem in uniprot_rev_values_for_uniprots_list:
            counts_occ[elem] = gn_study_count_dict.get(elem,0.0)  # count of studies with this gene name in df_fgn
                
        # get keys with max value:
        max_value = max(counts_occ.values())
        if max_value >0:
            gns_mostFreq41 = [k for k,v in counts_occ.items() if v == max_value]
        else:
            gns_mostFreq41 = np.nan
        if isinstance(gns_mostFreq41, list):
            if len(gns_mostFreq41)==1:
                return gns_mostFreq41[0]
            elif len(gns_mostFreq41) >1:
                return ";".join(gns_mostFreq41)
            else:
                return "#".join(['checkEntryOccurences1 case 4.1 ',index])
        else:
            return "@".join(uniprot_rev_values_for_uniprots_list)
        
    
    # 4.2: uniprot_unrev_dict
    elif len(uniprot_unrev_values_for_uniprots_list) >0:
        counts_occ = dict()
        for elem in uniprot_unrev_values_for_uniprots_list:
            counts_occ[elem] = gn_study_count_dict.get(elem,0.0)  # count of studies with this gene name in df_fgn
                
        # get keys with max value:
        max_value = max(counts_occ.values())
        if max_value >0:
            gns_mostFreq42 = [k for k,v in counts_occ.items() if v == max_value]
        else:
            gns_mostFreq42 = np.nan
        if isinstance(gns_mostFreq42, list):
            if len(gns_mostFreq42)==1:
                return gns_mostFreq42[0]
            elif len(gns_mostFreq42) >1:
                return ";".join(gns_mostFreq42)
            else:
                return "#".join(['checkEntryOccurences1 case 4.2 ',index])
        
        else:
            return "@".join(uniprot_unrev_values_for_uniprots_list)
            
            
    # 4.3: uniprot_ids_mrh_dict (based on df uniprot_gn)
    elif len(uniprot_ids_mrh_dict_values_for_uniprots_list) >0:
        counts_occ = dict()
        for elem in uniprot_ids_mrh_dict_values_for_uniprots_list:
            counts_occ[elem] = gn_study_count_dict.get(elem,0.0)  # count of studies with this gene name in df_fgn
                
        # get keys with max value:
        max_value = max(counts_occ.values())
        if max_value >0:
            gns_mostFreq43 = [k for k,v in counts_occ.items() if v == max_value]
        else:
            gns_mostFreq43 = np.nan
        if isinstance(gns_mostFreq43, list):
            if len(gns_mostFreq43)==1:
                return gns_mostFreq43[0]
            elif len(gns_mostFreq43) >1:
                return ";".join(gns_mostFreq43)
            else:
                return "#".join(['checkEntryOccurences1 case 4.3 ',index])
        
        else:
            return "@".join(uniprot_ids_mrh_dict_values_for_uniprots_list)
    
    
    
def get_gene_id_final(index,row):
    """
    Count occurences of ids to find the most frequent, consider id type and reliability (reviewed/unreviewed). Perform the main part of the ids alignment
    """
    
    uniprot = row['Uniprot']
    gene_name_unified = row['gene_name_unified']
    gene_names = row['gene_names']

    
    #case 4
    if ((isinstance(gene_names,float)) or (gene_name_unified == "NoMapping"))&(isinstance(uniprot,str)):
        if np.isnan(gene_names):
            uniprots_list0 = row["Uniprot"].replace(" ","").split(";")
            uniprots_list = [x.split("-")[0] for x in uniprots_list0 if x is not None] # check if its needed now!!!!

            if any(isinstance(i, list) for i in uniprots_list):
                uniprot_rev_values_for_uniprots_list0 = list(chain.from_iterable([uniprot_rev_dict.get(i) for i in uniprots_list]))
                uniprot_unrev_values_for_uniprots_list0 = list(chain.from_iterable([uniprot_unrev_dict.get(i) for i in uniprots_list]))
                uniprot_ids_mrh_dict_values_for_uniprots_list0 = list(chain.from_iterable([uniprot_ids_mrh_dict.get(i) for i in uniprots_list]))


            elif isinstance(uniprots_list, list):
                uniprot_rev_values_for_uniprots_list0 = [uniprot_rev_dict.get(i) for i in uniprots_list]
                uniprot_unrev_values_for_uniprots_list0 = [uniprot_unrev_dict.get(i) for i in uniprots_list]
                uniprot_ids_mrh_dict_values_for_uniprots_list0 = [uniprot_ids_mrh_dict.get(i) for i in uniprots_list]


            else:
                print("no uniprots")

            uniprot_rev_values_for_uniprots_list1 = [x for x in uniprot_rev_values_for_uniprots_list0 if x is not None]
            uniprot_rev_values_for_uniprots_list = list(chain.from_iterable(uniprot_rev_values_for_uniprots_list1))

            uniprot_unrev_values_for_uniprots_list1 = [x for x in uniprot_unrev_values_for_uniprots_list0 if x is not None]
            uniprot_unrev_values_for_uniprots_list = list(chain.from_iterable(uniprot_unrev_values_for_uniprots_list1))

            uniprot_ids_mrh_dict_values_for_uniprots_list1 = [x for x in uniprot_ids_mrh_dict_values_for_uniprots_list0 if x is not None]
            uniprot_ids_mrh_dict_values_for_uniprots_list = list(chain.from_iterable(uniprot_ids_mrh_dict_values_for_uniprots_list1))

            return get_uniprot_raw_data(index,row,uniprots_list,uniprot_rev_values_for_uniprots_list,uniprot_unrev_values_for_uniprots_list,uniprot_ids_mrh_dict_values_for_uniprots_list)
        
        
        else:
            print('CHECK NAN None', index)
    
    else:
        
        gene_names_list = row['gene_names'].replace(" ","").split(";")

        
        #case 3
        # 3.1: uniprot_rev_dict
        #uniprot_rev_values_for_uniprots_list =  #list(set([item for sublist in [uniprot_rev_dict.get(i) for i in uniprots_list] for item in sublist])) # type: gene names


        # 3.2: uniprot_unrev_dict
        #uniprot_unrev_values_for_uniprots_list =  #list(set([item for sublist in [uniprot_unrev_dict.get(i) for i in uniprots_list] for item in sublist]))  

        # 3.3: uniprot_ids_mrh_dict (based on df uniprot_gn)
        #uniprot_ids_mrh_dict_values_for_uniprots_list =  #list(set([item for sublist in [uniprot_ids_mrh_dict.get(i) for i in uniprots_list] for item in sublist]))  


        # cases 1 & 2
        if gene_name_unified in gene_names_list:
            return gene_name_unified
        

        
        elif isinstance(uniprot,str):
            uniprots_list0 = row["Uniprot"].replace(" ","").split(";")
            uniprots_list = [x.split("-")[0] for x in uniprots_list0 if x is not None]

            if any(isinstance(i, list) for i in uniprots_list):
                uniprot_rev_values_for_uniprots_list0 = list(chain.from_iterable([uniprot_rev_dict.get(i) for i in uniprots_list]))
                uniprot_unrev_values_for_uniprots_list0 = list(chain.from_iterable([uniprot_unrev_dict.get(i) for i in uniprots_list]))
                uniprot_ids_mrh_dict_values_for_uniprots_list0 = list(chain.from_iterable([uniprot_ids_mrh_dict.get(i) for i in uniprots_list]))


            elif isinstance(uniprots_list, list):
                uniprot_rev_values_for_uniprots_list0 = [uniprot_rev_dict.get(i) for i in uniprots_list]
                uniprot_unrev_values_for_uniprots_list0 = [uniprot_unrev_dict.get(i) for i in uniprots_list]
                uniprot_ids_mrh_dict_values_for_uniprots_list0 = [uniprot_ids_mrh_dict.get(i) for i in uniprots_list]


            else:
                print("no uniprots")

            uniprot_rev_values_for_uniprots_list1 = [x for x in uniprot_rev_values_for_uniprots_list0 if x is not None]
            uniprot_rev_values_for_uniprots_list = list(chain.from_iterable(uniprot_rev_values_for_uniprots_list1))

            uniprot_unrev_values_for_uniprots_list1 = [x for x in uniprot_unrev_values_for_uniprots_list0 if x is not None]
            uniprot_unrev_values_for_uniprots_list = list(chain.from_iterable(uniprot_unrev_values_for_uniprots_list1))

            uniprot_ids_mrh_dict_values_for_uniprots_list1 = [x for x in uniprot_ids_mrh_dict_values_for_uniprots_list0 if x is not None]
            uniprot_ids_mrh_dict_values_for_uniprots_list = list(chain.from_iterable(uniprot_ids_mrh_dict_values_for_uniprots_list1))



            # 3.1: uniprot_rev_dict
            if ((gene_name_unified in uniprot_rev_values_for_uniprots_list) or ( any_in(gene_names_list, uniprot_rev_values_for_uniprots_list) )  ):

                occurences01 = [x for x in uniprot_rev_values_for_uniprots_list if x in gene_names_list] # type: gene names
                occurences1 = [x for x in occurences01 if x is not None]

                occurences02 = [x for x in uniprot_rev_values_for_uniprots_list if x == gene_name_unified] # type: gene names
                occurences2 = [x for x in occurences02 if x is not None]

                occurences = occurences1+occurences2  # type: gene names 
                #occurences = [x for x in occurences1+occurences2 if x is not None]

                if len(occurences) == 1:
                    return occurences[0]

                elif len(occurences) > 1:
                    if len(occurences1) == 1:
                        return occurences1[0]
                    elif len(occurences1) > 1:
                        counts_occ = dict()
                        for elem in occurences1:
                            counts_occ[elem] = gn_study_count_dict.get(elem,0.0)  # count of studies with this gene name in df_fgn

                        # get keys with max value:
                        max_value = max(counts_occ.values())
                        if max_value >0:
                            gns_mostFreq = [k for k,v in counts_occ.items() if v == max_value]
                        else:
                            gns_mostFreq = np.nan
                        if isinstance(gns_mostFreq, list):
                            if len(gns_mostFreq)==1:
                                return gns_mostFreq[0]
                            elif len(gns_mostFreq) >1:
                                return ";".join(gns_mostFreq)
                            else:
                                return "#".join(['checkEntryOccurences1 case 3.1 ',index])


                    elif len(occurences2) > 1:
                        counts_occ = dict()
                        for elem in occurences2:
                            counts_occ[elem] = gn_study_count_dict.get(elem,0.0)  # count of studies with this gene name in df_fgn

                        # get keys with max value:
                        max_value = max(counts_occ.values())
                        if max_value >0:
                            gns_mostFreq = [k for k,v in counts_occ.items() if v == max_value]
                        else:
                            gns_mostFreq = np.nan
                        if isinstance(gns_mostFreq, list):
                            if len(gns_mostFreq)==1:
                                return gns_mostFreq[0]
                            elif len(gns_mostFreq) >1:
                                return ";".join(gns_mostFreq)
                            else:
                                return "#".join(['check entry occurences 2  case 3.1 ',index])


                    else:
                        return "#".join(["check for errors: occurences in case 3.1",index])

                else:
                    return "#".join(["check for errors in case 3.1",index])



            # 3.2: uniprot_unrev_dict    
            elif ((gene_name_unified in uniprot_unrev_values_for_uniprots_list) or ( any_in(gene_names_list, uniprot_unrev_values_for_uniprots_list) )  ):

                occurences012 = [x for x in uniprot_unrev_values_for_uniprots_list if x in gene_names_list] # type: gene names
                occurences12 = [x for x in occurences012 if x is not None]

                occurences022 = [x for x in uniprot_unrev_values_for_uniprots_list if x == gene_name_unified] # type: gene names
                occurences22 = [x for x in occurences022 if x is not None]

                occurences32 = occurences12+occurences22  # type: gene names 
                #occurences = [x for x in occurences1+occurences2 if x is not None]

                if len(occurences32) == 1:
                    return occurences32[0]

                elif len(occurences32) > 1:
                    if len(occurences12) == 1:
                        return occurences12[0]
                    elif len(occurences12) > 1:
                        counts_occ = dict()
                        for elem in occurences12:
                            counts_occ[elem] = gn_study_count_dict.get(elem,0.0)  # count of studies with this gene name in df_fgn

                        # get keys with max value:
                        max_value = max(counts_occ.values())
                        if max_value >0:
                            gns_mostFreq32 = [k for k,v in counts_occ.items() if v == max_value]
                        else:
                            gns_mostFreq32 = np.nan
                        if isinstance(gns_mostFreq32, list):
                            if len(gns_mostFreq32)==1:
                                return gns_mostFreq32[0]
                            elif len(gns_mostFreq32) >1:
                                return ";".join(gns_mostFreq32)
                            else:
                                return "#".join(['check entry occurences 1 case 3.2 ',index])


                    elif len(occurences22) > 1:
                        counts_occ = dict()
                        for elem in occurences22:
                            counts_occ[elem] = gn_study_count_dict.get(elem,0.0)  # count of studies with this gene name in df_fgn

                        # get keys with max value:
                        max_value = max(counts_occ.values())
                        if max_value >0:
                            gns_mostFreq32 = [k for k,v in counts_occ.items() if v == max_value]
                        else:
                            gns_mostFreq32 = np.nan
                        if isinstance(gns_mostFreq32, list):
                            if len(gns_mostFreq32)==1:
                                return gns_mostFreq32[0]
                            elif len(gns_mostFreq32) >1:
                                return ";".join(gns_mostFreq32)
                            else:
                                return "#".join(['check entry occurences 2 case 3.2 ',index])


                    else:
                        return "#".join(["check for errors: occurences in case 3.2",index])

                else:
                    return "#".join(["check for errors in case 3.2",index])



            # 3.3: uniprot_ids_mrh_dict (based on df uniprot_gn)
            elif ((gene_name_unified in uniprot_ids_mrh_dict_values_for_uniprots_list) or ( any_in(gene_names_list, uniprot_ids_mrh_dict_values_for_uniprots_list) )  ):

                occurences013 = [x for x in uniprot_ids_mrh_dict_values_for_uniprots_list if x in gene_names_list] # type: gene names
                occurences13 = [x for x in occurences013 if x is not None]

                occurences023 = [x for x in uniprot_rev_values_for_uniprots_list if x == gene_name_unified] # type: gene names
                occurences23 = [x for x in occurences023 if x is not None]

                occurences33 = occurences13+occurences23  # type: gene names 
                #occurences = [x for x in occurences1+occurences2 if x is not None]

                if len(occurences33) == 1:
                    return occurences33[0]

                elif len(occurences33) > 1:
                    if len(occurences13) == 1:
                        return occurences13[0]
                    elif len(occurences13) > 1:
                        counts_occ = dict()
                        for elem in occurences13:
                            counts_occ[elem] = gn_study_count_dict.get(elem,0.0)  # count of studies with this gene name in df_fgn

                        # get keys with max value:
                        max_value = max(counts_occ.values())
                        if max_value >0:
                            gns_mostFreq33 = [k for k,v in counts_occ.items() if v == max_value]
                        else:
                            gns_mostFreq33 = np.nan
                        if isinstance(gns_mostFreq33, list):
                            if len(gns_mostFreq33)==1:
                                return gns_mostFreq33[0]
                            elif len(gns_mostFreq33) >1:
                                return ";".join(gns_mostFreq33)
                            else:
                                return "#".join(['check entry occurences 1 case 3.3 ',index])


                    elif len(occurences2) > 1:
                        counts_occ = dict()
                        for elem in occurences2:
                            counts_occ[elem] = gn_study_count_dict.get(elem,0.0)  # count of studies with this gene name in df_fgn

                        # get keys with max value:
                        max_value = max(counts_occ.values())
                        if max_value >0:
                            gns_mostFreq33 = [k for k,v in counts_occ.items() if v == max_value]
                        else:
                            gns_mostFreq33 = np.nan
                        if isinstance(gns_mostFreq33, list):
                            if len(gns_mostFreq33)==1:
                                return gns_mostFreq33[0]
                            elif len(gns_mostFreq33) >1:
                                return ";".join(gns_mostFreq33)
                            else:
                                return "#".join(['check entry occurences 2 case 3.3 ',index])

                    else:
                        return "#".join(["check for errors: occurences in case 3.3",index])

                else:
                    return "#".join(["check for errors in case 3.3",index])
            
            else:
                print("attention ",index)
                return "&".join([gene_names, gene_name_unified])

