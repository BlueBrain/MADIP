#!/usr/bin/env python
# coding: utf-8
"""
This file is part of MADIP: Molecular Atlas Data Integration Pipeline

This module loads the data for step_1_collect_protein_data.ipynb jupyter notebook.
Age in days is only approximate and not involved in the downstream analysis. Qualitative age category is defined based on the original data sources.

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

def get_hamezah_2019_dataframe():
    """
    Return pandas dataframe for Hamezah 2019
    :return:
    pandas.core.frame.DataFrame: dataframe containing Hamezah 2019 data.
    """
    print("Importing Hamezah 2019 pandas dataframe.")
    # NEWDATA 1. Hamezah 2019 (Mice Hippocampus, Medial Prefrontal Cortex, and Striatum)

    hamezah_2019_f = pd.ExcelFile('../data/source_data/Hameezah_2019_formatted.xlsx')

    hamezah_2019_hippocampus = hamezah_2019_f.parse("Hippocampus", index_col=None)
    hamezah_2019_pfc = hamezah_2019_f.parse("MedialPrefrontalCortex", index_col=None)
    hamezah_2019_striatum = hamezah_2019_f.parse("Striatum", index_col=None)

    hamezah_2019_hippocampus['location'] = 'hippocampus'
    hamezah_2019_pfc['location'] = 'cortex'  # prefrontal cortex
    hamezah_2019_striatum['location'] = 'striatum'

    hamezah_2019 = pd.concat([hamezah_2019_hippocampus, hamezah_2019_pfc, hamezah_2019_striatum])
    hamezah_2019['Gene names'] = hamezah_2019['Gene names'].str.upper()

    hamezah_2019['Calc: LFQ intensity WT'] = 2 ** hamezah_2019['Average LFQ intensity (Log2) WT-Ctrl']
    hamezah_2019['Calc: LFQ intensity Alzheimer Tg'] = 2 ** hamezah_2019['Average LFQ intensity (Log2) Tg-ctrl']  #

    hamezah_2019 = hamezah_2019.reset_index(drop=True)

    hamezah_2019['Gene names'] = hamezah_2019['Gene names'].str.upper()

    hamezah_2019 = hamezah_2019.drop(columns=['Protein names', 'Average LFQ intensity (Log2) WT-Ctrl',
                                              'Average LFQ intensity (Log2) Tg-ctrl',
                                              'Number of\npeptides',
                                              'Maxquant\nscore', 'MS/MS\ncount', 'p-value', 'q-value'
                                              ])
    hamezah_2019 = hamezah_2019.rename(columns={'Protein\naccession': 'Uniprot',
                                                'Gene names': 'gene_names',
                                                'Molecular\nweight (kDa)': 'molecular_weight_kDa',
                                                'Calc: LFQ intensity WT': 'LFQintensity_WT',
                                                'Calc: LFQ intensity Alzheimer Tg': 'LFQintensity_Alzheimer'
                                                })

    hamezah_2019_df = pd.wide_to_long(hamezah_2019, stubnames='LFQintensity',
                                      i=['Uniprot', 'gene_names', 'molecular_weight_kDa', 'location'],
                                      j='condition', sep='_', suffix=r'\w+')

    hamezah_2019_df = hamezah_2019_df.reset_index()

    hamezah_2019_df['Study'] = 'Hamezah 2019'
    hamezah_2019_df['Organism'] = 'mouse'
    hamezah_2019_df['Age_days'] = 365 + 3 * 30 + 21  
    # Five-month-old mice ... for a duration of 10 months -> 15 months

    hamezah_2019_df['raw_data_units'] = 'LFQintensity'

    hamezah_2019_df = hamezah_2019_df.rename(columns={'LFQintensity': 'raw_data'})
    return hamezah_2019_df


def get_hamezah_2018_dataframe():
    """
    Return pandas dataframe for Hamezah 2018
    :return:
     pandas.core.frame.DataFrame: dataframe containing Hamezah 2018 data
    """
    # ### Hamezah_2018
    print("Importing pandas Hamezah 2018 dataframe.")
    hamezah_2018f = pd.ExcelFile('../data/source_data/1-s2.0-S0531556518303097-mmc2.xlsx')

    hamezah_2018_hippocampus = hamezah_2018f.parse('Sheet1')
    hamezah_2018_pfc = hamezah_2018f.parse('Sheet2')  # medial prefrontal cortex
    hamezah_2018_striatum = hamezah_2018f.parse('Sheet3')

    hamezah_2018_hippocampus['location'] = 'hippocampus'
    hamezah_2018_pfc['location'] = 'cortex'  # prefrontal cortex
    hamezah_2018_striatum['location'] = 'striatum'

    hamezah_2018 = pd.concat([hamezah_2018_hippocampus, hamezah_2018_pfc, hamezah_2018_striatum])
    hamezah_2018['Gene names'] = hamezah_2018['Gene names'].str.upper()

    hamezah_2018['Calc: LFQ intensity 14 months'] = 2 ** hamezah_2018['Average LFQ intensity (Log2) 14 months']
    hamezah_2018['Calc: LFQ intensity 18 months'] = 2 ** hamezah_2018['Average LFQ intensity (Log2) 18 months']
    hamezah_2018['Calc: LFQ intensity 23 months'] = 2 ** hamezah_2018['Average LFQ intensity (Log2) 23 months']
    hamezah_2018['Calc: LFQ intensity 27 months'] = 2 ** hamezah_2018['Average LFQ intensity (Log2) 27 months']

    hamezah_2018 = hamezah_2018.reset_index(drop=True)

    hamezah_2018['Gene names'] = hamezah_2018['Gene names'].str.upper()

    hamezah_2018 = hamezah_2018.drop(columns=['Protein names', 'Average LFQ intensity (Log2) 14 months',
                                              'Average LFQ intensity (Log2) 18 months',
                                              'Average LFQ intensity (Log2) 23 months',
                                              'Average LFQ intensity (Log2) 27 months',
                                              'Number of\npeptides',
                                              'Maxquant\nscore', 'MS/MS\ncount', 'ANOVA significant'
                                              ])

    hamezah_2018 = hamezah_2018.rename(columns={'Protein\naccession': 'Uniprot',
                                                'Gene names': 'gene_names',
                                                'Molecular\nweight (kDa)': 'molecular_weight_kDa',

                                                'Calc: LFQ intensity 14 months': 'LFQintensity_14months',
                                                'Calc: LFQ intensity 18 months': 'LFQintensity_18months',
                                                'Calc: LFQ intensity 23 months': 'LFQintensity_23months',
                                                'Calc: LFQ intensity 27 months': 'LFQintensity_27months'
                                                })

    hamezah_2018_df = pd.wide_to_long(hamezah_2018, stubnames='LFQintensity',
                                      i=['Uniprot', 'gene_names', 'molecular_weight_kDa', 'location'],
                                      j='sample_id', sep='_', suffix=r'\w+')

    hamezah_2018_df = hamezah_2018_df.reset_index()

    hamezah_2018_df['Study'] = 'Hamezah 2018'
    hamezah_2018_df['Organism'] = 'rat'

    hamezah_2018_df.loc[hamezah_2018_df['sample_id'] == '14months', 'Age_days'] = 365 + 2 * 30 + 21  # 446 # 365 + 2*30 + 21 # '14 months'
    hamezah_2018_df.loc[hamezah_2018_df['sample_id'] == '18months', 'Age_days'] = 365 + 6 * 30 + 21  # 365 + 6*30 +21 # '18 months'
    hamezah_2018_df.loc[hamezah_2018_df['sample_id'] == '23months', 'Age_days'] = 721  # 365*2 -30 +21 # '23 months'
    hamezah_2018_df.loc[hamezah_2018_df['sample_id'] == '27months', 'Age_days'] = 841  # 365*2 + 30*3 +21 # '27 months'

    hamezah_2018_df['raw_data_units'] = 'LFQintensity'

    hamezah_2018_df = hamezah_2018_df.rename(columns={'LFQintensity': 'raw_data'})
    return hamezah_2018_df


def get_chuang_2018_dataframe():
    """
    Return pandas dataframe for Chuang 2018
    :return:
     pandas.core.frame.DataFrame: dataframe containing Chuang 2018 data.
    """
    print("Importing Chuang 2018 pandas dataframe.")
    chuang2018f = pd.ExcelFile('../data/source_data/Supporting File S-3_The lists of the proteins identified  in the axon and whole-cell samples.xlsx')

    chuang2018_axon = chuang2018f.parse('Axon samples, 2548', skiprows=3, index_col=None)
    chuang2018_wholecell = chuang2018f.parse('Whole-cell samples, 2752', skiprows=3, index_col=None)

    chuang2018_wholecell = chuang2018_wholecell.drop(
        ['Fasta headers', 'Number of proteins', 'Peptides axon', 'Peptides whole-cell',
         'Razor + unique peptides axon', 'Razor + unique peptides whole-cell',
         'Score', 'Sequence coverage axon [%]',
         'Sequence coverage whole-cell [%]',
         'Fasta headers.1', 'Number of proteins.1',
         'Peptides axon.1', 'Peptides whole-cell.1',
         'Razor + unique peptides axon.1',
         'Razor + unique peptides whole-cell.1', 'Score.1',
         'Sequence coverage axon [%].1', 'Sequence coverage whole-cell [%].1'], axis=1)
    chuang2018_axon = chuang2018_axon.drop(['Fasta headers', 'Number of proteins', 'Peptides axon',
                                            'Peptides whole-cell',
                                            'Razor + unique peptides axon', 'Razor + unique peptides whole-cell',
                                            'Score', 'Sequence coverage axon [%]',
                                            'Sequence coverage whole-cell [%]',
                                            'Fasta headers.1', 'Number of proteins.1',
                                            'Peptides axon.1', 'Peptides whole-cell.1',
                                            'Razor + unique peptides axon.1',
                                            'Razor + unique peptides whole-cell.1', 'Score.1',
                                            'Sequence coverage axon [%].1', 'Sequence coverage whole-cell [%].1'],
                                           axis=1)

    chuang2018_axon = chuang2018_axon.rename(columns={'GN': 'gene_names',
                                                      'Accession': 'Uniprot',
                                                      'Protein IDs': 'Experiment1:Protein IDs',
                                                      'Protein IDs.1': 'Experiment2:Protein IDs.1',
                                                      'iBAQ axon': 'iBAQ_Experiment1',
                                                      'iBAQ axon.1': 'iBAQ_Experiment2'})

    chuang2018_wholecell = chuang2018_wholecell.rename(columns={'GN': 'gene_names',
                                                                'Accession': 'Uniprot',
                                                                'Protein IDs': 'Experiment1:Protein IDs',
                                                                'Protein IDs.1': 'Experiment2:Protein IDs.1',
                                                                'iBAQ whole-cell': 'iBAQ_Experiment1',
                                                                'iBAQ whole-cell.1': 'iBAQ_Experiment2'})

    chuang2018_axon['gene_names'] = chuang2018_axon['gene_names'].str.upper()
    chuang2018_wholecell['gene_names'] = chuang2018_wholecell['gene_names'].str.upper()

    chuang2018_axon['location'] = 'axon'
    chuang2018_wholecell['location'] = 'neurons'  # 'neuron_whole_cell'

    chuang2018 = pd.concat([chuang2018_axon, chuang2018_wholecell], sort=False)
    chuang2018 = chuang2018.reset_index(drop=True)

    chuang2018 = chuang2018.drop(['Description', 'Experiment1:Protein IDs', 'Experiment2:Protein IDs.1'], axis=1)

    chuang2018 = chuang2018.rename(columns={'Mol. weight [kDa]': 'molecular_weight_kDa'})

    chuang2018_df = pd.wide_to_long(chuang2018, stubnames='iBAQ',
                                    i=['Uniprot', 'gene_names', 'molecular_weight_kDa', 'location'],
                                    j='sample_id', sep='_', suffix=r'\w+')

    chuang2018_df = chuang2018_df.reset_index()

    chuang2018_df['Study'] = 'Chuang 2018'
    chuang2018_df['Organism'] = 'rat'
    chuang2018_df['Age_days'] = 18  # E18

    chuang2018_df['Age_cat'] = 'embr'

    chuang2018_df['raw_data_units'] = 'iBAQ'

    chuang2018_df = chuang2018_df.rename(columns={'iBAQ': 'raw_data'})
    return chuang2018_df


def get_duda_2018_dataframe():
    """
    Return pandas dataframe for Duda 2018
    :return:
     pandas.core.frame.DataFrame: dataframe containing Duda 2018 data.
    """
    print("Importing Duda 2018 pandas dataframe.")
    dudaf = pd.ExcelFile('../data/source_data/dataProt.xlsx')

    duda_hippocampus = dudaf.parse('hippocampus')
    duda_cerebellum = dudaf.parse('cerebellum')
    duda_cortex = dudaf.parse('cortex')

    # fill merged cells with the same values
    # merged cells processed manually to avoid artefacts
    # duda_hippocampus = duda_hippocampus.fillna(method='ffill')
    # duda_cerebellum = duda_cerebellum.fillna(method='ffill')
    # duda_cortex = duda_cortex.fillna(method='ffill')

    duda_hippocampus['location'] = 'hippocampus'
    duda_cerebellum['location'] = 'cerebellum'
    duda_cortex['location'] = 'cortex'

    duda = pd.concat([duda_hippocampus, duda_cerebellum, duda_cortex], sort=False)
    duda = duda.reset_index(drop=True)

    duda['gene_names'] = duda['gene_names'].str.upper()

    # young = 1 month; old = 12 months

    duda['duplicated'] = duda.duplicated(subset=['Young Mean concentration',
                                                 'Adult Mean concentration', 'location'], keep=False)

    duda = duda.drop(columns='duplicated')
    duda = duda.rename(columns={'Young Mean concentration': 'MeanConcentration_young',
                                'Adult Mean concentration': 'MeanConcentration_adult'})

    duda_2018_df = pd.wide_to_long(duda, stubnames='MeanConcentration',
                                   i=['gene_names', 'location'],
                                   j='condition', sep='_', suffix=r'\w+')

    duda_2018_df = duda_2018_df.reset_index()

    duda_2018_df['Study'] = 'Duda 2018'
    duda_2018_df['Organism'] = 'mouse'

    duda_2018_df.loc[duda_2018_df['condition'] == 'young', 'Age_days'] = 51  # P30 = 21 embryonic days + 30 postnatal
    duda_2018_df.loc[duda_2018_df['condition'] == 'adult', 'Age_days'] = 386  # 365 + 21 embryonic days #'12 months'

    duda_2018_df['raw_data_units'] = 'Mean concentration [mol/(g total protein)]'
    duda_2018_df = duda_2018_df.rename(columns={'MeanConcentration': 'raw_data'})
    return duda_2018_df


def get_krogager_2018_dataframe():
    """
    Return pandas dataframe for Krogager 2018
    :return:
     krogager_df :pandas.core.frame.DataFrame: dataframe containing Krogager 2018 data.
    """
    print("Importing Krogager 2018 pandas dataframe.")
    krogagerf = pd.ExcelFile('../data/source_data/MouseBrainProteomeKrogager2018_supp.xlsx')
    krogager = krogagerf.parse('Sheet1')

    krogager = krogager.drop(columns=['Significant (S0:1, FDR:0.05)', '-LOG(P-value)',
                                      'Log2(SORT Output / Control Output)', 'Protein names',
                                      'Intensity', 'MS/MS Count'])

    # in this case we combine samples due to many NaN in individual samples

    col1 = krogager.loc[:, ['Log2(LFQ) Control Output 1', 'Log2(LFQ) Control Output 2', 'Log2(LFQ) Control Output 3',
                            'Log2(LFQ) Control Output 4', 'Log2(LFQ) Control Output 5', 'Log2(LFQ) Control Output 6']]
    krogager['Log2(LFQ) Control median'] = col1.median(axis=1)

    col2 = krogager.loc[:, ['Log2(LFQ) SORT Output 1', 'Log2(LFQ) SORT Output 2', 'Log2(LFQ) SORT Output 3']]
    krogager['Log2(LFQ) SORT median'] = col2.median(axis=1)

    krogager['LFQintensity_control'] = 2 ** krogager['Log2(LFQ) Control median']
    krogager['LFQintensity_SORT'] = 2 ** krogager['Log2(LFQ) SORT median']

    krogager['Gene names'] = krogager['Gene names'].str.upper()

    krogager = krogager.rename(columns={'Gene names': 'gene_names',
                                        'Majority protein IDs': 'Uniprot'
                                        })

    krogager_drop = krogager.drop(['Log2(LFQ) Control Output 1',
                                   'Log2(LFQ) Control Output 2', 'Log2(LFQ) Control Output 3',
                                   'Log2(LFQ) Control Output 4', 'Log2(LFQ) Control Output 5',
                                   'Log2(LFQ) Control Output 6', 'Log2(LFQ) SORT Output 1',
                                   'Log2(LFQ) SORT Output 2', 'Log2(LFQ) SORT Output 3',
                                   'Log2(LFQ) Control median', 'Log2(LFQ) SORT median'], axis=1)

    krogager_df = pd.wide_to_long(krogager_drop, stubnames='LFQintensity',
                                  i=['Uniprot', 'gene_names'],
                                  j='condition', sep='_', suffix=r'\w+')

    krogager_df = krogager_df.reset_index()

    krogager_df['Study'] = 'Krogager 2018'
    krogager_df['Organism'] = 'mouse'

    krogager_df['Age_days'] = 13 * 7 + 21  # 13*7 +21 # 10 weeks + 2weeks after surgery + 1 week treatment
    krogager_df.loc[krogager_df['condition'] == 'SORT', 'location'] = 'neurons'  # striatum neurons
    krogager_df.loc[krogager_df['condition'] == 'control', 'location'] = 'striatum'  # striatum neurons

    krogager_df['raw_data_units'] = 'LFQintensity'

    krogager_df = krogager_df.rename(columns={'LFQintensity': 'raw_data'})
    return krogager_df


def get_hosp_2017_dataframe():
    """
    Return pandas dataframe for Hosp 2017
    :return:
     pandas.core.frame.DataFrame: dataframe containing Hosp 2017 data.
    """
    print("Importing Hosp 2017 pandas dataframe. This can last a while.")
    hosp_solf = pd.ExcelFile('../data/source_data/1-s2.0-S2211124717315772-mmc2.xlsx')
    hosp_sol = hosp_solf.parse('S1A_soluble proteome') 

    hosp_sol2f = pd.ExcelFile('../data/source_data/1-s2.0-S2211124717315772-mmc3.xlsx')
    hosp_sol2 = hosp_sol2f.parse('S2A_CSF_proteome')

    hosp_insolf = pd.ExcelFile('../data/source_data/1-s2.0-S2211124717315772-mmc4.xlsx')
    hosp_insol = hosp_insolf.parse('S3A_insolube_proteome_data')

    hosp_sol = hosp_sol.drop(
        ['GOBP name', 'GOMF name', 'GOCC name', 'KEGG name', 'Pfam', 'GSEA', 'Keywords', 'Corum', 'Peptides',
         'Razor + unique peptides', 'Razor + unique peptides', 'Sequence coverage [%]',
         'Unique + razor sequence coverage [%]', 'Unique sequence coverage [%]', 'Q-value'], axis=1)
    hosp_sol = hosp_sol[hosp_sol.columns.drop(list(hosp_sol.filter(regex=r'LFQ')))]
    hosp_sol = hosp_sol[hosp_sol.columns.drop(list(hosp_sol.filter(regex=r'_R6\/2_')))]

    hosp_sol2 = hosp_sol2.drop(
        ['GOBP name', 'GOMF name', 'GOCC name', 'KEGG name', 'Pfam', 'GSEA', 'Fasta headers', 'Corum', 'Peptides',
         'Razor + unique peptides', 'Razor + unique peptides', 'Sequence coverage [%]',
         'Unique + razor sequence coverage [%]', 'Unique sequence coverage [%]', 'Q-value'], axis=1)
    hosp_sol2 = hosp_sol2[hosp_sol2.columns.drop(list(hosp_sol2.filter(regex=r'LFQ')))]
    hosp_sol2 = hosp_sol2[hosp_sol2.columns.drop(list(hosp_sol2.filter(regex=r'_R6\/2_')))]

    hosp_insol = hosp_insol.drop(
        ['GOBP name', 'GOMF name', 'GOCC name', 'KEGG name', 'Pfam', 'GSEA', 'Fasta headers', 'Corum',
         'Coiled-coil domain',
         'LCR motif', 'polyQ domain', 'coiled-coil length', 'LCR length', 'polyQ length', 'Peptides',
         'Razor + unique peptides', 'Razor + unique peptides', 'Sequence coverage [%]',
         'Unique + razor sequence coverage [%]', 'Unique sequence coverage [%]', 'Q-value'], axis=1)
    hosp_insol = hosp_insol[hosp_insol.columns.drop(list(hosp_insol.filter(regex=r'LFQ')))]
    hosp_insol = hosp_insol[hosp_insol.columns.drop(list(hosp_insol.filter(regex=r'_R6\/2_')))]

    hosp_sol['Gene names'] = hosp_sol['Gene names'].str.upper()
    hosp_sol2['Gene names'] = hosp_sol2['Gene names'].str.upper()
    hosp_insol['Gene names'] = hosp_insol['Gene names'].str.upper()

    ###

    hosp_sol = hosp_sol.rename(columns={'Gene names': 'gene_names',
                                        'Majority protein IDs': 'Uniprot',
                                        'Mol. weight [kDa]': 'molecular_weight_kDa'})

    hosp_sol = hosp_sol.drop([
        'Protein IDs', 'Protein names', 'Unique peptides', 'Intensity', 'MS/MS Count',
        'iBAQ', 'iBAQ library'], axis=1)

    hosp_sol2 = hosp_sol2.rename(columns={'Gene names': 'gene_names',
                                          'Majority protein IDs': 'Uniprot',
                                          'Mol. weight [kDa]': 'molecular_weight_kDa'})

    hosp_sol2 = hosp_sol2.drop([
        'Protein IDs', 'Protein names', 'Unique peptides', 'Score', 'Intensity', 'MS/MS Count',
        'iBAQ_total'], axis=1)

    hosp_insol = hosp_insol.rename(columns={'Gene names': 'gene_names',
                                            'Majority protein IDs': 'Uniprot',
                                            'Mol. weight [kDa]': 'molecular_weight_kDa'})

    hosp_insol = hosp_insol.drop([
        'Protein IDs', 'Protein names', 'Unique peptides', 'Score', 'Intensity', 'MS/MS Count',
        'iBAQ'], axis=1)

    hosp_sol.columns = ['Uniprot', 'gene_names', 'molecular_weight_kDa', 'iBAQ_5wWTce1', 'iBAQ_5wWTce2', 'iBAQ_5wWTce3',
                        'iBAQ_5wWTce4',
                        'iBAQ_5wWTco1', 'iBAQ_5wWTco2', 'iBAQ_5wWTco3',
                        'iBAQ_5wWTco4', 'iBAQ_5wWThc1', 'iBAQ_5wWThc2',
                        'iBAQ_5wWThc3', 'iBAQ_5wWThc4', 'iBAQ_5wWTst1',
                        'iBAQ_5wWTst2', 'iBAQ_5wWTst3', 'iBAQ_5wWTst4',
                        'iBAQ_8wWTce1', 'iBAQ_8wWTce2', 'iBAQ_8wWTce3',
                        'iBAQ_8wWTco1', 'iBAQ_8wWTco2', 'iBAQ_8wWTco3',
                        'iBAQ_8wWThc1', 'iBAQ_8wWThc2', 'iBAQ_8wWThc3',
                        'iBAQ_8wWTst1', 'iBAQ_8wWTst2', 'iBAQ_8wWTst3',
                        'iBAQ_12wWTce1', 'iBAQ_12wWTce2', 'iBAQ_12wWTce3',
                        'iBAQ_12wWTco1', 'iBAQ_12wWTco2', 'iBAQ_12wWTco3',
                        'iBAQ_12wWThc1', 'iBAQ_12wWThc2', 'iBAQ_12wWThc3',
                        'iBAQ_12wWTst1', 'iBAQ_12wWTst2', 'iBAQ_12wWTst3']

    hosp_sol2.columns = ['Uniprot', 'gene_names', 'molecular_weight_kDa', 'iBAQ_5wWT1', 'iBAQ_5wWT2', 'iBAQ_5wWT3',
                         'iBAQ_8wWT1', 'iBAQ_8wWT2',
                         'iBAQ_8wWT3', 'iBAQ_12wWT1', 'iBAQ_12wWT2', 'iBAQ_12wWT3'
                         ]

    hosp_insol.columns = ['Uniprot', 'gene_names', 'molecular_weight_kDa',
                          'iBAQ_5wWTce1', 'iBAQ_5wWTce2', 'iBAQ_5wWTce3', 'iBAQ_5wWTce4',
                          'iBAQ_5wWTco1', 'iBAQ_5wWTco2', 'iBAQ_5wWTco3',
                          'iBAQ_5wWTco4', 'iBAQ_5wWThc1', 'iBAQ_5wWThc2',
                          'iBAQ_5wWThc3', 'iBAQ_5wWThc4', 'iBAQ_5wWTst1',
                          'iBAQ_5wWTst2', 'iBAQ_5wWTst3', 'iBAQ_5wWTst4',
                          'iBAQ_8wWTce1', 'iBAQ_8wWTce2', 'iBAQ_8wWTce3',
                          'iBAQ_8wWTco1', 'iBAQ_8wWTco2', 'iBAQ_8wWTco3',
                          'iBAQ_8wWThc1', 'iBAQ_8wWThc2', 'iBAQ_8wWThc3',
                          'iBAQ_8wWTst1', 'iBAQ_8wWTst2', 'iBAQ_8wWTst3',
                          'iBAQ_12wWTce1', 'iBAQ_12wWTce2', 'iBAQ_12wWTce3',
                          'iBAQ_12wWTco1', 'iBAQ_12wWTco2', 'iBAQ_12wWTco3',
                          'iBAQ_12wWThc1', 'iBAQ_12wWThc2', 'iBAQ_12wWThc3',
                          'iBAQ_12wWTst1', 'iBAQ_12wWTst2', 'iBAQ_12wWTst3',
                          'iBAQ_5wWTcex/yiBAQ_inc5wWTce',
                          'iBAQ_5wWTcox/yiBAQ_inc5wWTco',
                          'iBAQ_5wWThcx/yiBAQ_inc5wWThc',
                          'iBAQ_5wWTstx/yiBAQ_inc5wWTst',
                          'iBAQ_8wWTcex/yiBAQ_inc8wWTce',
                          'iBAQ_8wWTcox/yiBAQ_inc8wWTco',
                          'iBAQ_8wWThcx/yiBAQ_inc8wWThc',
                          'iBAQ_8wWTstx/yiBAQ_inc8wWTst',
                          'iBAQ_12wWTcex/yiBAQ_inc12wWTce',
                          'iBAQ_12wWTcox/yiBAQ_inc12wWTco',
                          'iBAQ_12wWThcx/yiBAQ_inc12wWThc',
                          'iBAQ_12wWTstx/yiBAQ_inc12wWTst']

    hosp_insol = hosp_insol.drop(['iBAQ_5wWTcex/yiBAQ_inc5wWTce',
                                  'iBAQ_5wWTcox/yiBAQ_inc5wWTco',
                                  'iBAQ_5wWThcx/yiBAQ_inc5wWThc',
                                  'iBAQ_5wWTstx/yiBAQ_inc5wWTst',
                                  'iBAQ_8wWTcex/yiBAQ_inc8wWTce',
                                  'iBAQ_8wWTcox/yiBAQ_inc8wWTco',
                                  'iBAQ_8wWThcx/yiBAQ_inc8wWThc',
                                  'iBAQ_8wWTstx/yiBAQ_inc8wWTst',
                                  'iBAQ_12wWTcex/yiBAQ_inc12wWTce',
                                  'iBAQ_12wWTcox/yiBAQ_inc12wWTco',
                                  'iBAQ_12wWThcx/yiBAQ_inc12wWThc',
                                  'iBAQ_12wWTstx/yiBAQ_inc12wWTst'], axis=1)

    hosp_sol_df = pd.wide_to_long(hosp_sol, stubnames='iBAQ',
                                  i=['Uniprot', 'gene_names', 'molecular_weight_kDa'],
                                  j='sample_id', sep='_', suffix=r'\w+')

    hosp_sol_df = hosp_sol_df.reset_index()

    hosp_sol_df['Study'] = 'Hosp 2017, soluble'

    hosp_sol_df['Organism'] = 'mouse'

    hosp_sol_df['raw_data_units'] = 'iBAQ'

    hosp_sol_df = hosp_sol_df.rename(columns={'iBAQ': 'raw_data'})

    hosp_sol2_df = pd.wide_to_long(hosp_sol2, stubnames='iBAQ',
                                   i=['Uniprot', 'gene_names', 'molecular_weight_kDa'],
                                   j='sample_id', sep='_', suffix=r'\w+')

    hosp_sol2_df = hosp_sol2_df.reset_index()

    hosp_sol2_df['Study'] = 'Hosp 2017, CSF'

    hosp_sol2_df['Organism'] = 'mouse'

    hosp_sol2_df['raw_data_units'] = 'iBAQ'

    hosp_sol2_df = hosp_sol2_df.rename(columns={'iBAQ': 'raw_data'})

    hosp_insol_df = pd.wide_to_long(hosp_insol, stubnames='iBAQ',
                                    i=['Uniprot', 'gene_names', 'molecular_weight_kDa'],
                                    j='sample_id', sep='_', suffix=r'\w+')

    hosp_insol_df = hosp_insol_df.reset_index()

    hosp_insol_df['Study'] = 'Hosp 2017, insoluble'

    hosp_insol_df['Organism'] = 'mouse'

    hosp_insol_df['raw_data_units'] = 'iBAQ'

    hosp_insol_df = hosp_insol_df.rename(columns={'iBAQ': 'raw_data'})

    hosp_3 = pd.concat([hosp_sol_df, hosp_sol2_df, hosp_insol_df
                        ], ignore_index=True, sort=False)

    hosp_3.loc[hosp_3['sample_id'].isin(['5wWTce1', '5wWTce2', '5wWTce3', '5wWTce4', '5wWTco1', '5wWTco2',
                                         '5wWTco3', '5wWTco4', '5wWThc1', '5wWThc2', '5wWThc3', '5wWThc4',
                                         '5wWTst1', '5wWTst2', '5wWTst3', '5wWTst4', '5wWT1', '5wWT2',
                                         '5wWT3']), 'Age_days'] = 35 + 21  # 35 +21 #'5 weeks'

    hosp_3.loc[hosp_3['sample_id'].isin(['8wWTce1', '8wWTce2',
                                         '8wWTce3', '8wWTco1', '8wWTco2', '8wWTco3', '8wWThc1', '8wWThc2',
                                         '8wWThc3', '8wWTst1', '8wWTst2', '8wWTst3', '8wWT1', '8wWT2',
                                         '8wWT3']), 'Age_days'] = 56 + 21  # 56 +21 #'8 weeks'

    hosp_3.loc[
        hosp_3['sample_id'].isin(['12wWTce1', '12wWTce2', '12wWTce3', '12wWTco1', '12wWTco2', '12wWTco3', '12wWThc1',
                                  '12wWThc2', '12wWThc3', '12wWTst1', '12wWTst2', '12wWTst3', '12wWT1', '12wWT2',
                                  '12wWT3']), 'Age_days'] = 12 * 7 + 21  # 12*7 +21 #'12 weeks'

    hosp_3.loc[hosp_3['sample_id'].isin(['5wWTce1', '5wWTce2', '5wWTce3', '5wWTce4', '8wWTce1', '8wWTce2',
                                         '8wWTce3', '12wWTce1', '12wWTce2',
                                         '12wWTce3']), 'location'] = 'cerebellum'

    hosp_3.loc[hosp_3['sample_id'].isin(['5wWTco1', '5wWTco2', '5wWTco3', '5wWTco4', '8wWTco1', '8wWTco2', '8wWTco3',
                                         '12wWTco1', '12wWTco2', '12wWTco3']), 'location'] = 'cortex'

    hosp_3.loc[hosp_3['sample_id'].isin(['5wWThc1', '5wWThc2', '5wWThc3', '5wWThc4',
                                         '8wWThc1', '8wWThc2', '8wWThc3', '12wWThc1', '12wWThc2',
                                         '12wWThc3']), 'location'] = 'hippocampus'

    hosp_3.loc[hosp_3['sample_id'].isin(
        ['5wWTst1', '5wWTst2', '5wWTst3', '5wWTst4', '8wWTst1', '8wWTst2', '8wWTst3', '12wWTst1', '12wWTst2',
         '12wWTst3']), 'location'] = 'striatum'

    hosp_3.loc[hosp_3['sample_id'].isin(
        ['5wWT1', '5wWT2', '5wWT3', '8wWT1', '8wWT2', '8wWT3', '12wWT1', '12wWT2', '12wWT3']), 'location'] = 'csf'
    return hosp_3


def get_itzhak_2017_dataframe():
    """
    Return pandas dataframe for Itzhak 2017
    :return:
     pandas.core.frame.DataFrame: dataframe containing Itzhak 2017 data.
    """
    print("Importing itzhak 2017 pandas dataframe. This can last a while.")

    itzhak_concf = pd.ExcelFile('../data/source_data/1-s2.0-S2211124717311889-mmc4.xlsx')
    itzhak_conc = itzhak_concf.parse('Mouse Neuron Spatial Proteome')

    itzhak_conc = itzhak_conc.rename(columns={'Lead gene name': 'gene_names',
                                              'Majority protein IDs': 'Uniprot',
                                              'Prediction': 'location',
                                              'MW (kD)': 'molecular_weight_kDa',
                                              'Median cellular concentration [nM]': 'raw_data'})

    itzhak_conc = itzhak_conc.replace([np.inf, -np.inf], np.nan)

    itzhak_conc = itzhak_conc.drop(['Canonical lead protein ID', 'Protein names', 'Marker protein?',
                                    'SVM score', 'Subcellular distribution',
                                    'Prediction confidence class', 'Median copy number', 'Max copy number',
                                    'Min copy number',
                                    'Abundance percentile (copy numbers)', 'Max cellular concentration [nM]',
                                    'Min cellular concentration [nM]', 'ppm of total cell mass',
                                    'Average Nuclear Fraction', 'Average Membrane',
                                    'Average Cytosol Fraction', 'MSMS count'], axis=1)

    itzhak_conc['raw_data_units'] = 'Median cellular concentration [nM]'

    itzhak_conc['Study'] = 'Itzhak 2017'
    itzhak_conc['Organism'] = 'mouse'
    itzhak_conc['Age_days'] = 15
    itzhak_conc['Age_cat'] = 'embr'
    itzhak_conc.loc[itzhak_conc['location'].isna(),'location'] = 'subcellular not specified'
    return itzhak_conc


def get_beltran_2016_dataframe():
    """
    Return pandas dataframe for Beltran 2016
    :return:
     pandas.core.frame.DataFrame: dataframe containing Beltran 2016 data.
    """
    print("Importing Beltran 2016 pandas dataframe")

    beltranf = pd.ExcelFile('../data/source_data/1-s2.0-S2405471216302897-mmc2.xlsx')
    beltran = beltranf.parse('TableS1A', skiprows=2, index_col=None)

    beltran = beltran.rename(columns={'Gene': 'gene_names',
                                      'Uniprot Accession': 'Uniprot',
                                      '24hpi': 'infected 24hpi',
                                      '48hpi': 'infected 48hpi',
                                      '72hpi': 'infected 72hpi', '96hpi': 'infected 96hpi',
                                      '120hpi': 'infected 120hpi', '24hpi.1': 'Uninfected (Mock) 24hpi'})

    beltran = beltran.drop(['infected 24hpi', 'infected 48hpi', 'infected 72hpi', 'infected 96hpi', 'infected 120hpi'],
                           axis=1)

    beltran = beltran.rename(columns={'Uninfected (Mock) 24hpi': 'raw_data'})

    beltran = beltran[beltran['Organism'] == 'Homo sapiens (Human)']

    beltran['raw_data_units'] = 'iBAQ'

    beltran['Study'] = 'Beltran 2016'
    beltran['Age_days'] = 0
    beltran['Age_cat'] = 'embr'

    # Protein abundance (iBAQ values)

    beltran['Organism'] = 'human'
    beltran['gene_names'] = beltran['gene_names'].str.upper()

    beltranLocf = pd.ExcelFile('../data/source_data/1-s2.0-S2405471216302897-mmc5.xlsx')
    beltranLoc = beltranLocf.parse('Label-Free Localization', skiprows=2, index_col=None)

    beltranLoc = beltranLoc[beltranLoc['Organism'] == 'Homo sapiens (Human)']
    beltranLoc = beltranLoc[['Uniprot Accession', 'Gene', '24hpi.1']]

    beltranLoc = beltranLoc.rename(
        columns={'Uniprot Accession': 'Uniprot', 'Gene': 'gene_names', '24hpi.1': 'location'})

    beltranLoc['gene_names'] = beltranLoc['gene_names'].str.upper()
    beltranfin = pd.merge(beltran, beltranLoc, on=['Uniprot', 'gene_names'], how='inner')

    beltranfin.loc[beltranfin['location'].isna(), 'location'] = 'subcellular not specified'
    return beltranfin


def get_sharma_2015_dataframe():
    """
    Return pandas dataframe for Sharma 2015
    :return:
     pandas.core.frame.DataFrame: dataframe containing Sharma 2015 data.
    """
    print("Importing Sharma 2015 pandas dataframe. This can last a while.")
    sharma4fo = pd.ExcelFile('../data/source_data/nn.4160-S4.xlsx')
    sharma4o = sharma4fo.parse('For Suppl table 2 proteinGroups', skiprows=2)

    sharmaf = pd.ExcelFile('../data/source_data/nn.4160-S7.xlsx')
    sharma = sharmaf.parse('Sheet1')

    sharma4o.columns = ['gene_names', 'Protein names',
                        'Log2LFQintensity_IsolatedAstrocytes',
                        'Log2LFQintensity_IsolatedMicroglia',
                        'Log2LFQintensity_IsolatedNeurons',
                        'Log2LFQintensity_IsolatedOligodendrocytes',
                        'Log2LFQintensity_Brain',
                        'Log2LFQintensity_Brainstem',
                        'Log2LFQintensity_Cerebellum',
                        'Log2LFQintensity_CorpusCallosum',
                        'Log2LFQintensity_MotorCortex',
                        'Log2LFQintensity_OlfactoryBulb',
                        'Log2LFQintensity_OpticNerve',
                        'Log2LFQintensity_PrefrontalCortex',
                        'Log2LFQintensity_Striatum',
                        'Log2LFQintensity_Thalamus',
                        'Log2LFQintensity_VentralHippocampus',
                        'Log2LFQintensity_CerebellumP05',
                        'Log2LFQintensity_CerebellumP14',
                        'Log2LFQintensity_CerebellumP24',
                        'Standard Deviation Isolated Astrocytes', 'Standard Deviation Isolated Microglia',
                        'Standard Deviation Isolated Neurons', 'Standard Deviation Isolated Oligodendrocytes',
                        'Standard Deviation Brain',
                        'Standard Deviation Brainstem', 'Standard Deviation Cerebellum',
                        'Standard Deviation Corpus Callosum', 'Standard Deviation MotorCortex',
                        'Standard Deviation Olfactory Bulb', 'Standard Deviation Optic Nerve',
                        'Standard Deviation Prefrontal Cortex', 'Standard Deviation Striatum',
                        'Standard Deviation Thalamus', 'Standard Deviation Ventral Hippocampus',
                        'Standard Deviation Cerebellum P05',
                        'Standard Deviation Cerebellum P14', 'Standard Deviation Cerebellum P24', 'Peptides',
                        'Sequence coverage [%]', 'molecular_weight_kDa', 'Score',
                        'Uniprot']

    sharma4o = sharma4o.drop(columns=['Protein names',
                                      'Standard Deviation Isolated Astrocytes',
                                      'Standard Deviation Isolated Microglia',
                                      'Standard Deviation Isolated Neurons',
                                      'Standard Deviation Isolated Oligodendrocytes',
                                      'Standard Deviation Brain', 'Standard Deviation Brainstem',
                                      'Standard Deviation Cerebellum', 'Standard Deviation Corpus Callosum',
                                      'Standard Deviation MotorCortex', 'Standard Deviation Olfactory Bulb',
                                      'Standard Deviation Optic Nerve',
                                      'Standard Deviation Prefrontal Cortex', 'Standard Deviation Striatum',
                                      'Standard Deviation Thalamus', 'Standard Deviation Ventral Hippocampus',
                                      'Standard Deviation Cerebellum P05',
                                      'Standard Deviation Cerebellum P14',
                                      'Standard Deviation Cerebellum P24', 'Peptides',
                                      'Sequence coverage [%]', 'Score'])

    sharma.columns = ['gene_names', 'Protein names', 'LFQintensity_adultMicroglia1',
                      'LFQintensity_adultMicroglia2', 'LFQintensity_adultMicroglia3',
                      'LFQintensity_youngMicroglia1', 'LFQintensity_youngMicroglia2',
                      'LFQintensity_youngMicroglia3', 'LFQintensity_Astrocytes1',
                      'LFQintensity_Astrocytes2', 'LFQintensity_Astrocytes3',
                      'LFQintensity_Neuronsdiv051',
                      'LFQintensity_Neuronsdiv052',
                      'LFQintensity_Neuronsdiv053',
                      'LFQintensity_Neuronsdiv101',
                      'LFQintensity_Neuronsdiv102',
                      'LFQintensity_Neuronsdiv103',
                      'LFQintensity_Neuronsdiv151',
                      'LFQintensity_Neuronsdiv152',
                      'LFQintensity_Neuronsdiv153',
                      'LFQintensity_Oligodendrocytesdiv11',
                      'LFQintensity_Oligodendrocytesdiv12',
                      'LFQintensity_Oligodendrocytesdiv13',
                      'LFQintensity_Oligodendrocytesdiv251',
                      'LFQintensity_Oligodendrocytesdiv252',
                      'LFQintensity_Oligodendrocytesdiv253',
                      'LFQintensity_Oligodendrocytesdiv41',
                      'LFQintensity_Oligodendrocytesdiv42',
                      'LFQintensity_Oligodendrocytesdiv43', 'PEP',
                      'molecular_weight_kDa', 'Sequence coverage [%]', 'Protein IDs',
                      'Uniprot']

    sharma = sharma.drop(['Protein names', 'PEP', 'Sequence coverage [%]', 'Protein IDs'], axis=1)

    sharma4o_df = pd.wide_to_long(sharma4o, stubnames='Log2LFQintensity',
                                  i=['Uniprot', 'gene_names', 'molecular_weight_kDa'],
                                  j='sample_id', sep='_', suffix=r'\w+')

    sharma4o_df = sharma4o_df.reset_index()

    sharma4o_df['raw_data'] = 2 ** sharma4o_df['Log2LFQintensity']

    sharma4o_df['Study'] = 'Sharma 2015, isolated'

    sharma4o_df['Organism'] = 'mouse'

    sharma4o_df['raw_data_units'] = 'LFQintensity'

    sharma4o_df = sharma4o_df.drop('Log2LFQintensity', axis=1)

    sharma_df = pd.wide_to_long(sharma, stubnames='LFQintensity',
                                i=['Uniprot', 'gene_names', 'molecular_weight_kDa'],
                                j='sample_id', sep='_', suffix=r'\w+')

    sharma_df = sharma_df.reset_index()

    sharma_df['Study'] = 'Sharma 2015, cultured'

    sharma_df['Organism'] = 'mouse'

    sharma_df['raw_data_units'] = 'LFQintensity'

    sharma_df['Age_days'] = 0  # 'cultured cells'

    sharma_df = sharma_df.rename(columns={'LFQintensity': 'raw_data'})

    sharma_df.loc[sharma_df['sample_id'].isin(['Neuronsdiv051',
                                               'Neuronsdiv052', 'Neuronsdiv053', 'Neuronsdiv101', 'Neuronsdiv102',
                                               'Neuronsdiv103', 'Neuronsdiv151', 'Neuronsdiv152',
                                               'Neuronsdiv153']), 'location'] = 'neurons'

    sharma_df.loc[sharma_df['sample_id'].isin(['Astrocytes1', 'Astrocytes2', 'Astrocytes3']), 'location'] = 'astrocytes'

    sharma_df.loc[sharma_df['sample_id'].isin(['adultMicroglia1', 'adultMicroglia2', 'adultMicroglia3',
                                               'youngMicroglia1', 'youngMicroglia2',
                                               'youngMicroglia3']), 'location'] = 'microglia'

    sharma_df.loc[sharma_df['sample_id'].isin(['Oligodendrocytesdiv11', 'Oligodendrocytesdiv12',
                                               'Oligodendrocytesdiv13', 'Oligodendrocytesdiv251',
                                               'Oligodendrocytesdiv252', 'Oligodendrocytesdiv253',
                                               'Oligodendrocytesdiv41', 'Oligodendrocytesdiv42',
                                               'Oligodendrocytesdiv43']), 'location'] = 'oligodendrocytes'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['IsolatedAstrocytes', 'IsolatedMicroglia', 'IsolatedNeurons',
                                                   'IsolatedOligodendrocytes']), 'Age_days'] = 29  # 8 + 21 # 'P8'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['Brain', 'Brainstem', 'Cerebellum',
                                                   'CorpusCallosum', 'MotorCortex', 'OlfactoryBulb', 'OpticNerve',
                                                   'PrefrontalCortex', 'Striatum', 'Thalamus',
                                                   'VentralHippocampus', ]), 'Age_days'] = 81  # 60 + 21 'P60'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['CerebellumP05']), 'Age_days'] = 26  # 5 + 21 # 'P5'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['CerebellumP14']), 'Age_days'] = 35  # 14 + 21 # 'P14'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['CerebellumP24']), 'Age_days'] = 45  # 24 + 21 # 'P24'

    ###

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['IsolatedAstrocytes']), 'location'] = 'astrocytes'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['IsolatedMicroglia']), 'location'] = 'microglia'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['IsolatedNeurons']), 'location'] = 'neurons'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['IsolatedOligodendrocytes']), 'location'] = 'oligodendrocytes'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['Brain']), 'location'] = 'brain'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['Brainstem']), 'location'] = 'brainstem'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(
        ['Cerebellum', 'CerebellumP05', 'CerebellumP14', 'CerebellumP24']), 'location'] = 'cerebellum'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['CorpusCallosum']), 'location'] = 'corpus callosum'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['MotorCortex', 'PrefrontalCortex']), 'location'] = 'cortex'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['VentralHippocampus']), 'location'] = 'hippocampus'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['OlfactoryBulb']), 'location'] = 'olfactory bulb'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['Striatum']), 'location'] = 'striatum'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['OpticNerve']), 'location'] = 'optic nerve'

    sharma4o_df.loc[sharma4o_df['sample_id'].isin(['Thalamus']), 'location'] = 'thalamus'
    return sharma4o_df, sharma_df


def get_wisniewski_2015_dataframe():
    """
    Return pandas dataframe for Wisńiewski 2015
    :return:
     pandas.core.frame.DataFrame: dataframe containing Wisniewski 2015 data.
    """
    print("Importing Wisńiewski 2015 pandas dataframe")

    wisniewskif = pd.ExcelFile('../data/source_data/pr5b00276_si_001.xlsx')
    wisniewski = wisniewskif.parse('Sheet1')

    wisniewski = wisniewski[~((wisniewski['Protein concentration (mol/g protein) Brain1'] == 0) &
                              (wisniewski['Protein concentration (mol/g protein) Brain2'] == 0) &
                              (wisniewski['Protein concentration (mol/g protein) Brain3'] == 0))]

    wisniewski = wisniewski.drop(
        ['Protein IDs', 'Protein names', 'Total protein Brain1', 'Total protein Brain2', 'Total protein Brain3'],
        axis=1)
    wisniewski = wisniewski.rename(columns={'Majority protein IDs': 'Uniprot',
                                            'Gene names': 'gene_names',
                                            'Protein concentration (mol/g protein) Brain1': 'Conc_1',
                                            'Protein concentration (mol/g protein) Brain2': 'Conc_2',
                                            'Protein concentration (mol/g protein) Brain3': 'Conc_3'})

    wisniewski_df = pd.wide_to_long(wisniewski, stubnames='Conc',
                                    i=['Uniprot', 'gene_names'],
                                    j='sample_id', sep='_', suffix=r'\w+')

    wisniewski_df = wisniewski_df.reset_index()

    wisniewski_df['Study'] = 'Wisniewski 2015'
    wisniewski_df['Organism'] = 'mouse'
    wisniewski_df['location'] = 'brain'

    wisniewski_df['Age_days'] = 10 * 7 + 21  # adult

    wisniewski_df['raw_data_units'] = 'Protein concentration (mol/g protein)'

    wisniewski_df = wisniewski_df.rename(columns={'Conc': 'raw_data'})

    return wisniewski_df


def get_han_2014_dataframe():
    """
    Return pandas dataframe for Han 2014
    :return:
     pandas.core.frame.DataFrame: dataframe containing  Han 2014 data.
    """
    print("Importing Han 2014 pandas dataframe")

    hanf = pd.ExcelFile('../data/source_data/pmic7746-sup-0001-tables1.xlsx')
    han = hanf.parse('Sheet1')

    # six technical replicates of two samples (conditioned media (CM) and whole-cell lysates (WCL))

    han = han.rename(columns={'Gene symbol': 'gene_names',
                              'LFQ intensity WCL_set1_tech1': 'LFQintensity_WCLset1tech1',
                              'LFQ intensity WCL_set1_tech2': 'LFQintensity_WCLset1tech2',
                              'LFQ intensity WCL_set1_tech3': 'LFQintensity_WCLset1tech3',
                              'LFQ intensity WCL_set2_tech1': 'LFQintensity_WCLset2tech1',
                              'LFQ intensity WCL_set2_tech2': 'LFQintensity_WCLset2tech2',
                              'LFQ intensity WCL_set2_tech3': 'LFQintensity_WCLset2tech3',
                              'LFQ intensity WCL_set3_tech1': 'LFQintensity_WCLset3tech1',
                              'LFQ intensity WCL_set3_tech2': 'LFQintensity_WCLset3tech2',
                              'LFQ intensity WCL_set3_tech3': 'LFQintensity_WCLset3tech3'
                              })

    han = han.drop(['Protein IDs', 'Majority protein IDs', 'Leading protein', 'Intensity', 'Intensity CM_set1_tech1',
                    'Intensity CM_set1_tech2', 'Intensity CM_set1_tech3',
                    'Intensity CM_set2_tech1', 'Intensity CM_set2_tech2',
                    'Intensity CM_set2_tech3', 'Intensity CM_set3_tech1',
                    'Intensity CM_set3_tech2', 'Intensity CM_set3_tech3', 'LFQ intensity CM_set1_tech1',
                    'LFQ intensity CM_set1_tech2', 'LFQ intensity CM_set1_tech3',
                    'LFQ intensity CM_set2_tech1', 'LFQ intensity CM_set2_tech2',
                    'LFQ intensity CM_set2_tech3', 'LFQ intensity CM_set3_tech1',
                    'LFQ intensity CM_set3_tech2', 'LFQ intensity CM_set3_tech3', 'MS/MS Count CM_set1_tech1',
                    'MS/MS Count CM_set1_tech2', 'MS/MS Count CM_set1_tech3',
                    'MS/MS Count CM_set2_tech1', 'MS/MS Count CM_set2_tech2',
                    'MS/MS Count CM_set2_tech3', 'MS/MS Count CM_set3_tech1',
                    'MS/MS Count CM_set3_tech2', 'MS/MS Count CM_set3_tech3',
                    'MS/MS Count WCL_set1_tech1', 'MS/MS Count WCL_set1_tech2',
                    'MS/MS Count WCL_set1_tech3', 'MS/MS Count WCL_set2_tech1',
                    'MS/MS Count WCL_set2_tech2', 'MS/MS Count WCL_set2_tech3',
                    'MS/MS Count WCL_set3_tech1', 'MS/MS Count WCL_set3_tech2',
                    'MS/MS Count WCL_set3_tech3', 'Intensity WCL_set1_tech1',
                    'Intensity WCL_set1_tech2', 'Intensity WCL_set1_tech3', 'Intensity WCL_set2_tech1',
                    'Intensity WCL_set2_tech2', 'Intensity WCL_set2_tech3', 'Intensity WCL_set3_tech1',
                    'Intensity WCL_set3_tech2', 'Intensity WCL_set3_tech3'], axis=1)

    han = han[han['gene_names'] != '-']

    han_df = pd.wide_to_long(han, stubnames='LFQintensity',
                             i=['Uniprot', 'gene_names'],
                             j='sample_id', sep='_', suffix=r'\w+')

    han_df = han_df.reset_index()

    han_df['Study'] = 'Han 2014'

    han_df['Organism'] = 'mouse'

    han_df['raw_data_units'] = 'LFQintensity'

    han_df['Age_days'] = 0  # 'cultured cells'

    han_df['location'] = 'astrocytes'

    han_df = han_df.rename(columns={'LFQintensity': 'raw_data'})
    return han_df


def get_geiger_2013_dataframe():
    """
    Return pandas dataframe for Geiger 2013
    :return:
     pandas.core.frame.DataFrame: dataframe containing Geiger 2013 data.
    """
    print("Importing Geiger 2013 pandas dataframe. This operation can last a while.")

    geigerf = pd.ExcelFile('../data/source_data/mcp.M112.024919-2.xlsx')
    geiger = geigerf.parse('Suppl Table S1', skiprows=1, index_col=None)

    geiger = geiger.drop(['Protein IDs', 'Protein names', 'Peptides', 'Razor + unique peptides', 'Unique peptides',
                          'Sequence coverage [%]', 'PEP',
                          'Ratio H/L normalized', 'Ratio H/L normalized Adrenal gland',
                          'Ratio H/L normalized Brain cortex',
                          'Ratio H/L normalized Brain medulla', 'Ratio H/L normalized Brown fat',
                          'Ratio H/L normalized Cerebellum', 'Ratio H/L normalized Colon',
                          'Ratio H/L normalized Diaphragm', 'Ratio H/L normalized Duodenum',
                          'Ratio H/L normalized Embryonic tissue', 'Ratio H/L normalized Eye',
                          'Ratio H/L normalized Heart', 'Ratio H/L normalized Ileum',
                          'Ratio H/L normalized Jejunum', 'Ratio H/L normalized Kidney cortex',
                          'Ratio H/L normalized Kidney medulla', 'Ratio H/L normalized Liver',
                          'Ratio H/L normalized Lung', 'Ratio H/L normalized Midbrain',
                          'Ratio H/L normalized Muscle', 'Ratio H/L normalized Olfactory bulb',
                          'Ratio H/L normalized Ovary', 'Ratio H/L normalized Pancreas',
                          'Ratio H/L normalized Salivary gland', 'Ratio H/L normalized Spleeen',
                          'Ratio H/L normalized Stomach', 'Ratio H/L normalized Thymus',
                          'Ratio H/L normalized Uterus', 'Ratio H/L normalized White fat',
                          'Intensity L Adrenal gland', 'Intensity L Brown fat', 'Intensity L Colon',
                          'Intensity L Diaphragm',
                          'Intensity L Duodenum', 'Intensity L Embryonic tissue',
                          'Intensity L Eye', 'Intensity L Heart', 'Intensity L Ileum',
                          'Intensity L Jejunum', 'Intensity L Kidney cortex',
                          'Intensity L Kidney medulla', 'Intensity L Liver', 'Intensity L Lung', 'Intensity L Muscle',
                          'Intensity L Ovary',
                          'Intensity L Pancreas', 'Intensity L Salivary gland',
                          'Intensity L Spleeen', 'Intensity L Stomach', 'Intensity L Thymus',
                          'Intensity L Uterus', 'Intensity L White fat'], axis=1)

    geiger = geiger.rename(columns={'Majority protein IDs': 'Uniprot',
                                    'Gene names': 'gene_names',
                                    'Mol. weight [kDa]': 'molecular_weight_kDa',
                                    'Intensity L Brain cortex': 'IntensityL_cortex',
                                    'Intensity L Brain medulla': 'IntensityL_medulla',
                                    'Intensity L Cerebellum': 'IntensityL_cerebellum',
                                    'Intensity L Midbrain': 'IntensityL_midbrain',
                                    'Intensity L Olfactory bulb': 'IntensityL_olfactorybulb'
                                    })

    # Lys-C

    geiger = geiger[~geiger['Uniprot'].isna()]

    geiger_df = pd.wide_to_long(geiger, stubnames='IntensityL',
                                i=['Uniprot', 'gene_names', 'molecular_weight_kDa'],
                                j='location', sep='_', suffix=r'\w+')

    geiger_df = geiger_df.reset_index()

    geiger_df['Study'] = 'Geiger 2013'

    geiger_df['Organism'] = 'mouse'

    geiger_df['raw_data_units'] = 'IntensityL'

    geiger_df['Age_days'] = 10 * 7 + 21  # adult

    geiger_df.loc[geiger_df['location'] == 'olfactorybulb', 'location'] = 'olfactory bulb'

    geiger_df = geiger_df.rename(columns={'IntensityL': 'raw_data'})
    return geiger_df


def get_bai_2020_dataframe():
    """
    Return pandas dataframe for Bai 2020
    :return:
     pandas.core.frame.DataFrame: dataframe containing Bai 2020 data.
    """
    print("Importing Bai 2020 pandas dataframe.")

    # mouse samples:
    bai2020_f = pd.ExcelFile('../data/source_data/1-s2.0-S089662731931058X-mmc7.xlsx')
    bai2020 = bai2020_f.parse('Sheet1', skiprows=4)

    bai2020 = bai2020.drop(['LPC1', 'LPC2', 'HPC1', 'HPC2', 'MCI1', 'MCI2', 'AD1', 'AD2', 'PSP1', 'PSP2'],
                           axis=1)  # these are human samples

    bai2020 = bai2020.rename(columns={'Human gene name': 'gene_names',
                                      'Human protein accession': 'Uniprot_human',
                                      'Mouse protein accession': 'Uniprot_mouse'
                                      })

    # bai2020[bai2020['Uniprot_mouse'].str.contains('CON_ENSBTAP00000024146')] #Q61838
    bai2020['Uniprot_human'] = bai2020['Uniprot_human'].str.split("|").str[1]
    bai2020['Uniprot_mouse'] = bai2020['Uniprot_mouse'].str.split("|").str[1]
    # [i for i in bai2020['Uniprot_mouse'].unique() if len(i)>6 ]
    bai2020.loc[bai2020['Uniprot_mouse']=='CON_ENSBTAP00000024146','Uniprot_mouse'] = 'Q61838' #CON_ENSBTAP00000024146 fix by human uniprot
    # [i for i in bai2020['Uniprot_human'].unique() if len(i)>6 ]
    bai2020['Uniprot'] = bai2020['Uniprot_mouse'] + ";" + bai2020['Uniprot_human'] 
    bai2020 = bai2020.drop(['Uniprot_mouse','Uniprot_human'],axis=1)

    bai2020 = bai2020.drop_duplicates(keep='first')
    bai2020[bai2020[['Uniprot','gene_names']].duplicated(keep=False)]

    bai2020_df = pd.wide_to_long(bai2020, stubnames='tmt', 
                                         i=['Uniprot','gene_names'],
                                         j='sample_id',sep='_',suffix='\w+')

    bai2020_df = bai2020_df.reset_index()

    bai2020_df['Organism'] = 'mouse' # there were human and mouse data, and here I used only mouse data, human data is imported by function get_human_samples_bai_2020_dataframe()
    bai2020_df['location'] = 'cortex'
    bai2020_df['raw_data_units'] = 'Protein Abundance (Summerized TMT Reporter Ion Intensities)'
    bai2020_df['Study'] = 'Bai 2020'
    bai2020_df = bai2020_df.rename(columns = {'tmt':'raw_data'})

    bai2020_df.loc[bai2020_df['sample_id'].isin(['WT3M1', 'WT3M2','AD3M1', 'AD3M2']),'Age_days'] = 111 # 3*30+21 # 3 months

    bai2020_df.loc[bai2020_df['sample_id'].isin(['WT6M1', 'WT6M2', 'WT6M3', 'WT6M4', 'AD6M1', 'AD6M2', 'AD6M3', 'AD6M4']),'Age_days'] = 201 # 6*30+21 # 6 months

    bai2020_df.loc[bai2020_df['sample_id'].isin(['WT12M1', 'WT12M2', 'AD12M1', 'AD12M2']),'Age_days'] = 386 # 365+21 # 12 months

    bai2020_df.loc[bai2020_df['sample_id'].isin(['WT3M1', 'WT3M2', 'WT6M1', 'WT6M2', 'WT6M3', 'WT6M4', 'WT12M1', 'WT12M2']),'condition'] = 'control'

    bai2020_df.loc[bai2020_df['sample_id'].isin(['AD3M1', 'AD3M2', 'AD6M1', 'AD6M2', 'AD6M3', 'AD6M4','AD12M1', 'AD12M2']),'condition'] = 'AD' # 5xFAD mice Alzheimer model

    return bai2020_df


def get_human_samples_bai_2020_dataframe():
    """
    Return pandas dataframe for human samples Bai 2020
    :return:
     pandas.core.frame.DataFrame: dataframe containing human samples 2020 data.
    """
    print("Importing human samples Bai 2020 pandas dataframe.")

    bai2020human_f = pd.ExcelFile('../data/source_data/1-s2.0-S089662731931058X-mmc7.xlsx')
    bai2020human = bai2020human_f.parse('Sheet1', skiprows=4)

    bai2020human = bai2020human.drop(['tmt_WT3M1', 'tmt_WT3M2', 'tmt_WT6M1', 'tmt_WT6M2', 'tmt_WT6M3',
                                      'tmt_WT6M4', 'tmt_WT12M1', 'tmt_WT12M2', 'tmt_AD3M1', 'tmt_AD3M2',
                                      'tmt_AD6M1', 'tmt_AD6M2', 'tmt_AD6M3', 'tmt_AD6M4', 'tmt_AD12M1',
                                      'tmt_AD12M2'], axis=1)  # these are mouse samples

    bai2020human = bai2020human.rename(columns={'Human gene name': 'gene_names',
                                                'Human protein accession': 'Uniprot_human',
                                                'Mouse protein accession': 'Uniprot_mouse'
                                                })

    bai2020human['Uniprot_human'] = bai2020human['Uniprot_human'].str.split("|").str[1]
    bai2020human['Uniprot_mouse'] = bai2020human['Uniprot_mouse'].str.split("|").str[1]

    # [i for i in bai2020human['Uniprot_mouse'].unique() if len(i)>6 ]
    # bai2020human[bai2020human['Uniprot_mouse'] == 'CON_ENSBTAP00000024146']

    bai2020human.loc[bai2020human[
                         'Uniprot_mouse'] == 'CON_ENSBTAP00000024146', 'Uniprot_mouse'] = 'Q61838'
    # CON_ENSBTAP00000024146 fix by human uniprot

    # [i for i in bai2020human['Uniprot_human'].unique() if len(i)>6 ]

    bai2020human['Uniprot'] = bai2020human['Uniprot_mouse'] + ";" + bai2020human['Uniprot_human']
    bai2020human = bai2020human.drop(['Uniprot_mouse', 'Uniprot_human'], axis=1)
    bai2020human = bai2020human.drop_duplicates(keep='first')
    bai2020human = bai2020human.rename(columns = { 'LPC1':'humandat_LPC1',
                                                 'LPC2':'humandat_LPC2',
                                                  'HPC1':'humandat_HPC1',
                                                  'HPC2':'humandat_HPC2',
                                                  'MCI1':'humandat_MCI1',
                                                  'MCI2':'humandat_MCI2',
                                                  'AD1':'humandat_AD1',
                                                  'AD2':'humandat_AD2',
                                                  'PSP1':'humandat_PSP1',
                                                  'PSP2':'humandat_PSP2'
                                                 })
    bai2020human_df = pd.wide_to_long(bai2020human, stubnames='humandat', 
                                     i=['Uniprot','gene_names'],
                                     j='sample_id',sep='_',suffix='\w+')

    bai2020human_df = bai2020human_df.reset_index()

    bai2020human_df['Organism'] = 'human' # there were human and mouse data, and I used only mouse data
    bai2020human_df['location'] = 'cortex' # frontal cortical samples of 100 human cases
    bai2020human_df['raw_data_units'] = 'Protein Abundance (Summerized TMT Reporter Ion Intensities)'

    bai2020human_df['Study'] = 'Bai 2020'

    bai2020human_df = bai2020human_df.rename(columns = {'humandat':'raw_data'})
    bai2020human_df.loc[bai2020human_df['sample_id'].isin(['LPC1', 'LPC2']),'condition'] = 'LPC: low pathology of plaques and tangles. AD'

    bai2020human_df.loc[bai2020human_df['sample_id'].isin(['HPC1', 'HPC2']),'condition'] = 'HPC: high Ab pathology but no detectable cognitive defects. AD'

    bai2020human_df.loc[bai2020human_df['sample_id'].isin(['MCI1', 'MCI2']),'condition'] = 'MCI: mild cognitive impairment with Ab pathology and a slight but measurable defect in cognition. AD'

    bai2020human_df.loc[bai2020human_df['sample_id'].isin(['AD1', 'AD2']),'condition'] = 'AD: late-stage AD with high pathology scores of plaques and tangles'

    bai2020human_df.loc[bai2020human_df['sample_id'].isin(['PSP1', 'PSP2']),'condition'] = 'PSP: progressive supranuclear palsy, another neurodegenerative disorder of tauopathy'

    bai2020human_df.loc[:,'Age_days'] =  'post-mortem'
    bai2020human_df = bai2020human_df.reset_index(drop=True)

    return bai2020human_df


def get_carlyle_2017_dataframe():
    """
    Return pandas dataframe for Carlyle 2017
    :return:
     pandas.core.frame.DataFrame: dataframe containing Carlyle 2017 data.
    """
    print("Importing Carlyle 2017 pandas dataframe.")
    # from 41593_2017_11_MOESM12_ESM data in 41593_2017_11_MOESM6_ESM is log10LFQ

    carlyle2017 = pd.read_csv('../data/source_data/41593_2017_11_MOESM6_ESM.txt', sep='\t')

    carlyle2017 = carlyle2017.drop(['EnsemblID'], axis=1)
    carlyle2017 = carlyle2017.rename(columns={'GeneSymbol': 'gene_names',
                                              'HSB118_AMY': 'Log10LFQ_118AMY',
                                              'HSB118_CBC': 'Log10LFQ_118CBC',
                                              'HSB118_DFC': 'Log10LFQ_118DFC',
                                              'HSB118_HIP': 'Log10LFQ_118HIP',
                                              'HSB118_MD': 'Log10LFQ_118MD',
                                              'HSB118_STR': 'Log10LFQ_118STR',
                                              'HSB118_V1C': 'Log10LFQ_118V1C',
                                              'HSB119_AMY': 'Log10LFQ_119AMY',
                                              'HSB121_AMY': 'Log10LFQ_121AMY',
                                              'HSB122_AMY': 'Log10LFQ_122AMY',
                                              'HSB122_CBC': 'Log10LFQ_122CBC',
                                              'HSB122_DFC': 'Log10LFQ_122DFC',
                                              'HSB122_HIP': 'Log10LFQ_122HIP',
                                              'HSB122_MD': 'Log10LFQ_122MD',
                                              'HSB122_STR': 'Log10LFQ_122STR',
                                              'HSB122_V1C': 'Log10LFQ_122V1C',
                                              'HSB123_AMY': 'Log10LFQ_123AMY',
                                              'HSB123_CBC': 'Log10LFQ_123CBC',
                                              'HSB123_DFC': 'Log10LFQ_123DFC',
                                              'HSB123_HIP': 'Log10LFQ_123HIP',
                                              'HSB123_STR': 'Log10LFQ_123STR',
                                              'HSB123_V1C': 'Log10LFQ_123V1C',
                                              'HSB126_AMY': 'Log10LFQ_126AMY',
                                              'HSB126_CBC': 'Log10LFQ_126CBC',
                                              'HSB126_DFC': 'Log10LFQ_126DFC',
                                              'HSB126_HIP': 'Log10LFQ_126HIP',
                                              'HSB126_MD': 'Log10LFQ_126MD',
                                              'HSB126_STR': 'Log10LFQ_126STR',
                                              'HSB126_V1C': 'Log10LFQ_126V1C',
                                              'HSB127_CBC': 'Log10LFQ_127CBC',
                                              'HSB127_DFC': 'Log10LFQ_127DFC',
                                              'HSB127_HIP': 'Log10LFQ_127HIP',
                                              'HSB127_MD': 'Log10LFQ_127MD',
                                              'HSB127_STR': 'Log10LFQ_127STR',
                                              'HSB127_V1C': 'Log10LFQ_127V1C',
                                              'HSB135_AMY': 'Log10LFQ_135AMY',
                                              'HSB135_CBC': 'Log10LFQ_135CBC',
                                              'HSB135_DFC': 'Log10LFQ_135DFC',
                                              'HSB135_HIP': 'Log10LFQ_135HIP',
                                              'HSB135_MD': 'Log10LFQ_135MD',
                                              'HSB135_STR': 'Log10LFQ_135STR',
                                              'HSB135_V1C': 'Log10LFQ_135V1C',
                                              'HSB136_AMY': 'Log10LFQ_136AMY',
                                              'HSB136_CBC': 'Log10LFQ_136CBC',
                                              'HSB136_DFC': 'Log10LFQ_136DFC',
                                              'HSB136_HIP': 'Log10LFQ_136HIP',
                                              'HSB136_MD': 'Log10LFQ_136MD',
                                              'HSB136_STR': 'Log10LFQ_136STR',
                                              'HSB136_V1C': 'Log10LFQ_136V1C',
                                              'HSB139_CBC': 'Log10LFQ_139CBC',
                                              'HSB139_DFC': 'Log10LFQ_139DFC',
                                              'HSB139_HIP': 'Log10LFQ_139HIP',
                                              'HSB139_MD': 'Log10LFQ_139MD',
                                              'HSB139_STR': 'Log10LFQ_139STR',
                                              'HSB139_V1C': 'Log10LFQ_139V1C',
                                              'HSB141_CBC': 'Log10LFQ_141CBC',
                                              'HSB141_DFC': 'Log10LFQ_141DFC',
                                              'HSB141_HIP': 'Log10LFQ_141HIP',
                                              'HSB141_MD': 'Log10LFQ_141MD',
                                              'HSB141_STR': 'Log10LFQ_141STR',
                                              'HSB141_V1C': 'Log10LFQ_141V1C',
                                              'HSB143_CBC': 'Log10LFQ_143CBC',
                                              'HSB143_DFC': 'Log10LFQ_143DFC',
                                              'HSB143_HIP': 'Log10LFQ_143HIP',
                                              'HSB143_MD': 'Log10LFQ_143MD',
                                              'HSB143_STR': 'Log10LFQ_143STR',
                                              'HSB143_V1C': 'Log10LFQ_143V1C',
                                              'HSB145_AMY': 'Log10LFQ_145AMY',
                                              'HSB145_CBC': 'Log10LFQ_145CBC',
                                              'HSB145_DFC': 'Log10LFQ_145DFC',
                                              'HSB145_HIP': 'Log10LFQ_145HIP',
                                              'HSB145_MD': 'Log10LFQ_145MD',
                                              'HSB145_STR': 'Log10LFQ_145STR',
                                              'HSB145_V1C': 'Log10LFQ_145V1C',
                                              'HSB173_AMY': 'Log10LFQ_173AMY',
                                              'HSB174_AMY': 'Log10LFQ_174AMY',
                                              'HSB175_AMY': 'Log10LFQ_175AMY'})

    carlyle2017_df = pd.wide_to_long(carlyle2017, stubnames='Log10LFQ',
                                     i=['gene_names'],
                                     j='sample_id', sep='_', suffix=r'\w+')

    carlyle2017_df = carlyle2017_df.reset_index()

    carlyle2017_df['Organism'] = 'human'
    carlyle2017_df['Study'] = 'Carlyle 2017'
    carlyle2017_df['raw_data'] = 10 ** carlyle2017_df['Log10LFQ']
    carlyle2017_df['raw_data_units'] = 'LFQintensity'

    carlyle2017_df = carlyle2017_df.drop(['Log10LFQ'], axis=1)

    carlyle_samples = pd.read_csv('../data/source_data/41593_2017_11_MOESM3_ESM.txt',sep='\t')
    
    # cerebellar cortex (CBC), striatum (STR), mediodorsal thalamic
    # nucleus (MD), amygdala (AMY), hippocampus (HIP), 
    # primary visual cortex (V1C), dorsolateral prefrontal cortex (DFC)

    #from early infancy (1 year after conception) to adulthood (42 years)
    #LFQs from 7 brain regions of 16 individuals

    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('AMY'),'location'] = 'amygdala'
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('CBC'),'location'] = 'cerebellum'
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('DFC'),'location'] = 'cortex'
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('HIP'),'location'] = 'hippocampus'
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('MD'),'location'] = 'thalamus'
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('STR'),'location'] = 'striatum'
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('V1C'),'location'] = 'cortex'

    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('118'),'Age_days'] = 1726
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('119'),'Age_days'] = 5741
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('121'),'Age_days'] = 386
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('122'),'Age_days'] = 631
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('123'),'Age_days'] = 13771
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('126'),'Age_days'] = 11216
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('127'),'Age_days'] = 7201
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('135'),'Age_days'] = 14866
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('136'),'Age_days'] = 8661
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('139'),'Age_days'] = 386
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('141'),'Age_days'] = 3186
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('143'),'Age_days'] = 996
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('145'),'Age_days'] = 13406
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('173'),'Age_days'] = 1361
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('174'),'Age_days'] = 3186
    carlyle2017_df.loc[carlyle2017_df['sample_id'].str.contains('175'),'Age_days'] = 4281

    # due to mismatches in some pairs of gene names (KRT86;GFAP as an example) and absence of Uniprot ids in the data, 8.2 % of gn (containing more than one gene id per entry) are removed:

    carlyle2017_df_filt = carlyle2017_df.loc[~carlyle2017_df['gene_names'].str.contains(';')]
    
    return carlyle2017_df_filt


def get_davis_2019_dataframe():
    """
    Return pandas dataframe for Davis 2019
    :return:
     pandas.core.frame.DataFrame: dataframe containing Davis 2019 data.
    """
    print("Importing Davis 2019 pandas dataframe.")

    # Davis 2019 (human brain, LFQ and iBAQ)
    # primary motor cortex (Betz cells) and cerebellar cortex (Purkinje cells) - individual neurons

    davis2019f = pd.ExcelFile('../data/source_data/pr8b00981_si_003.xlsx')
    davis2019 = davis2019f.parse('proteinGroups', skiprows=2)

    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'Peptide')))]
    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'Oxidation')))]
    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'Unique peptides')))]

    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'Razor')))]
    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'Identification type')))]

    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'Sequence coverage')))]
    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'Intensity')))]

    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'LFQ')))]

    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'Library')))]

    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'MS/MS')))]

    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'Molecular Buffer')))]

    davis2019 = davis2019[davis2019.columns.drop(list(davis2019.filter(regex=r'Molecular Cap')))]

    davis2019 = davis2019.drop(
        ['Protein IDs', 'Protein names', 'Fasta headers', 'Only identified by site', 'Reverse', 'Potential contaminant',
         'id', 'Mod. peptide IDs', 'Evidence IDs',
         'Number of proteins', 'Unique + razor sequence coverage [%]', 'Unique sequence coverage [%]',
         'Sequence length', 'Sequence lengths',
         'Fraction average', 'Fraction 1', 'Fraction 2', 'Fraction 3',
         'Fraction 10', 'Q-value', 'Score', 'iBAQ'], axis=1)
    
    davis2019 = davis2019.rename(columns = {
        'Majority protein IDs':'Uniprot', 'Gene names':'gene_names', 'Mol. weight [kDa]':'molecular_weight_kDa',
        'iBAQ Betz Cap SP3 1':'iBAQ_BetzCapSP31',
        'iBAQ Betz Cap SP3 2':'iBAQ_BetzCapSP32',
        'iBAQ Betz Cap SP3 3':'iBAQ_BetzCapSP33',
        'iBAQ Purkinje Cap SP3 1':'iBAQ_PurkinjeCapSP31',
        'iBAQ Purkinje Cap SP3 2':'iBAQ_PurkinjeCapSP32',
        'iBAQ Purkinje Cap SP3 3':'iBAQ_PurkinjeCapSP33',
        'iBAQ Purkinje InCap 100':'iBAQ_PurkinjeInCap100',
        'iBAQ Purkinje InCap 200':'iBAQ_PurkinjeInCap200',
        'iBAQ Purkinje InCap 400':'iBAQ_PurkinjeInCap400',
        'iBAQ Purkinje InCap 800':'iBAQ_PurkinjeInCap800',
        'iBAQ Purkinje Spin 10':'iBAQ_PurkinjeSpin10',
        'iBAQ Purkinje Spin 100':'iBAQ_PurkinjeSpin100',
        'iBAQ Purkinje Spin 200':'iBAQ_PurkinjeSpin200',
        'iBAQ Purkinje Spin 400':'iBAQ_PurkinjeSpin400',
        'iBAQ Purkinje Spin 50':'iBAQ_PurkinjeSpin50',
        'iBAQ Purkinje Spin 800':'iBAQ_PurkinjeSpin800'
    })

    davis2019_df = pd.wide_to_long(davis2019, stubnames='iBAQ', 
                                         i=['Uniprot', 'gene_names', 'molecular_weight_kDa',],
                                         j='sample_id',sep='_',suffix='\w+')

    davis2019_df = davis2019_df.reset_index()

    davis2019_df['Organism'] = 'human'

    davis2019_df['Study'] = 'Davis 2019'

    davis2019_df['raw_data_units'] = 'iBAQ'

    davis2019_df = davis2019_df.rename(columns = {'iBAQ':'raw_data'})

    davis2019_df['location'] = 'neurons'

    davis2019_df['Age_days'] = 'post-mortem'
    davis2019_df['Age_cat'] = 'post-mortem'
    return davis2019_df


def get_fecher_2019_dataframe():
    """
    Return pandas dataframe for Fecher 2019
    :return:
     pandas.core.frame.DataFrame: dataframe containing Fecher 2019 data.
    """
    print("Importing Fecher 2019 pandas dataframe.")

    fecher2019f = pd.ExcelFile('../data/source_data/41593_2019_479_MOESM3_ESM.xlsx')
    fecher2019 = fecher2019f.parse('Proteomics from PC, GC & A', skiprows=3)

    fecher2019 = fecher2019.drop(
        ['ENSG', 'Protein name', 'MitoCarta', 'Candidate', 'N   (IC GFP / IC Tom) Purkinje cell mito', 'Unnamed: 18',
         'Peptides Purkinje cell mito',
         'Sequence coverage [%] Purkinje cell mito', 'N  (IC GFP / IC Tom) Granule cell mito', 'Unnamed: 33',
         'Peptides Granule cell mito', 'Sequence coverage [%] Granule cell mito', 'N  (IC GFP / IC Tom)Astrocytic mito',
         'Unnamed: 50',
         'Peptides Astrocytic mito', 'Sequence coverage [%] Astrocytic mito',
         'Shared', 'Localization_LocTree3',
         'Core mitochondrial function', 'mitochondria related', 'MIM morbid #',
         'MolweightkDa_Purkinjecellmito',
         'MolweightkDa_Granulecellmito',
         'MolweightkDa_Astrocyticmito'], axis=1)

    fecher2019 = fecher2019.rename(columns={
        'Gene name': 'gene_names'

    })
    fecher2019_df = pd.wide_to_long(fecher2019, stubnames='log2LFQ', 
                                     i=['gene_names'],
                                     j='sample_id',sep='_',suffix='\w+')
    fecher2019_df = fecher2019_df.reset_index()
    fecher2019_df['Organism'] = 'mouse' #
    fecher2019_df['Study'] = 'Fecher 2019'

    fecher2019_df['raw_data_units'] = 'LFQintensity'
    fecher2019_df['raw_data'] = 2 ** fecher2019_df['log2LFQ']

    fecher2019_df['location'] = 'mitochondria'

    fecher2019_df['Age_days'] = 8 * 7 + 21 # 8- to 9-week-old male mice
    fecher2019_df = fecher2019_df.drop(['log2LFQ'],axis=1)
    fecher2019_df = fecher2019_df[~fecher2019_df['raw_data'].isna()]
    return fecher2019_df


def get_fornasiero_2018_dataframe():
    """
    Return pandas dataframe for Fornasiero 2018
    :return:
     pandas.core.frame.DataFrame: dataframe containing Fornasiero 2018 data.
    """
    print("Importing Fornasiero 2018 pandas dataframe.")

    fornasierof = pd.ExcelFile('../data/source_data/41467_2018_6519_MOESM3_ESM.xlsx')
    fornasiero2018 = fornasierof.parse('Data')

    fornasiero2018 = fornasiero2018[
        ['Uniprot identifier', 'gene_names', 'Protein abundance in brain cortex (average iBAQ expressed as log10+10)']]

    fornasiero2018['raw_data'] = 10 ** fornasiero2018[
        'Protein abundance in brain cortex (average iBAQ expressed as log10+10)'] - 10

    fornasiero2018['raw_data_units'] = 'iBAQ'

    fornasiero2018 = fornasiero2018.drop(['Protein abundance in brain cortex (average iBAQ expressed as log10+10)'],
                                         axis=1)

    fornasiero2018 = fornasiero2018.rename(columns={'Uniprot identifier': 'Uniprot'})

    fornasiero2018['Study'] = 'Fornasiero 2018'

    fornasiero2018['Organism'] = 'mouse'
    fornasiero2018['location'] = 'cortex'
    fornasiero2018['Age_cat'] = 'adult'  # adult by source
    fornasiero2018[
        'Age_days'] = 3.5 * 30 + 21 + 6*7 # 3.5 months # approximate, age in days is not used anywhere in the analysis, only age category appears in the analysis later  
    return fornasiero2018


def get_guergues_2019_dataframe():
    """
    Return pandas dataframe for Guergues 2019
    :return:
     pandas.core.frame.DataFrame: dataframe containing Guergues 2019 data.
    """
    print("Importing Guergues 2019 pandas dataframe.")

    guergues2019f = pd.ExcelFile('../data/source_data/pmic13102-sup-0003-tables1.xlsx')

    guergues2019 = guergues2019f.parse('proteinGroups')

    guergues2019 = guergues2019.drop(
        ['Protein IDs', 'Protein IDs 1', 'Peptide counts (all)', 'Peptide counts (razor+unique)',
         'Peptide counts (unique)', 'Protein names', 'Fasta headers', 'Number of proteins', 'Peptides',
         'Razor + unique peptides', 'Unique peptides', 'Peptides STrap_300K_1',
         'Peptides STrap_300K_2', 'Peptides STrap_300K_3',
         'Razor + unique peptides STrap_300K_1',
         'Razor + unique peptides STrap_300K_2',
         'Razor + unique peptides STrap_300K_3', 'Unique peptides STrap_300K_1',
         'Unique peptides STrap_300K_2', 'Unique peptides STrap_300K_3',
         'Sequence coverage [%]', 'Unique + razor sequence coverage [%]',
         'Unique sequence coverage [%]', 'Sequence length',
         'Sequence lengths', 'Q-value', 'Score',
         'Identification type STrap_300K_1', 'Identification type STrap_300K_2',
         'Identification type STrap_300K_3',
         'Sequence coverage STrap_300K_1 [%]',
         'Sequence coverage STrap_300K_2 [%]',
         'Sequence coverage STrap_300K_3 [%]', 'Intensity',
         'Intensity STrap_300K_1', 'Intensity STrap_300K_2',
         'Intensity STrap_300K_3',
         'no imp LFQ intensity STrap_300K_1',
         'no imp LFQ intensity STrap_300K_2',
         'no imp LFQ intensity STrap_300K_3', 'MS/MS count STrap_300K_1',
         'MS/MS count STrap_300K_2', 'MS/MS count STrap_300K_3', 'MS/MS count',
         'id', 'Peptide IDs', 'Peptide is razor', 'Mod. peptide IDs',
         'Evidence IDs', 'MS/MS IDs', 'Best MS/MS', 'Oxidation (M) site IDs',
         'Oxidation (M) site positions', 'Gene names sort check'], axis=1)

    guergues2019 = guergues2019.rename(columns={'Majority protein IDs': 'Uniprot',
                                                'Gene names': 'gene_names',
                                                'Mol. weight [kDa]': 'molecular_weight_kDa',

                                                'LFQ intensity STrap_300K_1': 'LFQintensity_STrap300K1',
                                                'LFQ intensity STrap_300K_2': 'LFQintensity_STrap300K2',
                                                'LFQ intensity STrap_300K_3': 'LFQintensity_STrap300K3'

                                                })

    guergues2019_df = pd.wide_to_long(guergues2019, stubnames='LFQintensity',
                                      i=['Uniprot', 'gene_names', 'molecular_weight_kDa'],
                                      j='sample_id', sep='_', suffix=r'\w+')

    guergues2019_df = guergues2019_df.reset_index()

    guergues2019_df['Organism'] = 'mouse'

    guergues2019_df['Study'] = 'Guergues 2019'

    guergues2019_df['raw_data_units'] = 'LFQintensity'

    guergues2019_df = guergues2019_df.rename(columns={'LFQintensity': 'raw_data'})

    guergues2019_df['location'] = 'microglia'

    guergues2019_df['Age_cat'] = 'adult'
    guergues2019_df[
        'Age_days'] = 8 * 7 + 21  # cultured cells from 8-week-old C57BL/6J mice #
    # by ref [4] R. C. McCarthy et al. J. Neuroinflamm. 2016
    return guergues2019_df


def get_mcketney_2019_dataframe():
    """
    Return pandas dataframe for McKetney 2019
    :return:
     pandas.core.frame.DataFrame: dataframe containing McKetney 2019 data.
    """
    print("Importing McKetney 2019 pandas dataframe.")

    mc_ketney2019f = pd.ExcelFile('../data/source_data/pr9b00004_si_002.xlsx')

    mc_ketney2019 = mc_ketney2019f.parse('LFQ intensities', skiprows=1)

    mc_ketney2019 = mc_ketney2019.drop(['Protein.IDs', 'Majority.protein.IDs', 'entrez_ids', 'geneNames',
                                        'Peptides',
                                        'Sequence coverage [%]', 'Peptide counts (all)',
                                        'Peptide counts (razor+unique)', 'Peptide counts (unique)'], axis=1)

    mc_ketney2019 = mc_ketney2019.rename(columns={'Reference Uniprot ID': 'Uniprot',
                                                  'geneSym': 'gene_names',
                                                  'LFQ.intensity.146_ AMY': 'LFQintensity_s146AMY',
                                                  'LFQ.intensity.146_ CNC': 'LFQintensity_s146CNC',
                                                  'LFQ.intensity.146_CBM': 'LFQintensity_s146CBM',
                                                  'LFQ.intensity.146_ECX': 'LFQintensity_s146ECX',
                                                  'LFQ.intensity.146_MFG': 'LFQintensity_s146MFG',
                                                  'LFQ.intensity.146_IPL': 'LFQintensity_s146IPL',
                                                  'LFQ.intensity.146_STG': 'LFQintensity_s146STG',
                                                  'LFQ.intensity.146_THA': 'LFQintensity_s146THA',
                                                  'LFQ.intensity.146_VCX': 'LFQintensity_s146VCX',
                                                  'LFQ.intensity.383_ AMY': 'LFQintensity_s383AMY',
                                                  'LFQ.intensity.383_ CNC': 'LFQintensity_s383CNC',
                                                  'LFQ.intensity.383_CBM': 'LFQintensity_s383CBM',
                                                  'LFQ.intensity.383_ECX': 'LFQintensity_s383ECX',
                                                  'LFQ.intensity.383_MFG': 'LFQintensity_s383MFG',
                                                  'LFQ.intensity.383_IPL': 'LFQintensity_s383IPL',
                                                  'LFQ.intensity.383_STG': 'LFQintensity_s383STG',
                                                  'LFQ.intensity.383_THA': 'LFQintensity_s383THA',
                                                  'LFQ.intensity.383_VCX': 'LFQintensity_s383VCX',
                                                  'LFQ.intensity.405_ AMY': 'LFQintensity_s405AMY',
                                                  'LFQ.intensity.405_ CNC': 'LFQintensity_s405CNC',
                                                  'LFQ.intensity.405_ECX': 'LFQintensity_s405ECX',
                                                  'LFQ.intensity.405_MFG': 'LFQintensity_s405MFG',
                                                  'LFQ.intensity.405_IPL': 'LFQintensity_s405IPL',
                                                  'LFQ.intensity.405_STG': 'LFQintensity_s405STG',
                                                  'LFQ.intensity.405_THA': 'LFQintensity_s405THA',
                                                  'LFQ.intensity.405_VCX': 'LFQintensity_s405VCX'
                                                  })

    mc_ketney2019_df = pd.wide_to_long(mc_ketney2019, stubnames='LFQintensity',
                                       i=['Uniprot', 'gene_names'],
                                       j='sample_id', sep='_', suffix=r'\w+')

    mc_ketney2019_df = mc_ketney2019_df.reset_index()

    mc_ketney2019_df['Organism'] = 'human'

    mc_ketney2019_df['Study'] = 'McKetney 2019'

    mc_ketney2019_df['raw_data_units'] = 'LFQintensity'

    mc_ketney2019_df = mc_ketney2019_df.rename(columns={'LFQintensity': 'raw_data'})

    mc_ketney2019_df['Age_cat'] = 'adult'  # “adult” samples (23−40 years old)
    mc_ketney2019_df['Age_days'] = 31 * 365  # “adult” samples (23−40 years old)

    mc_ketney2019_df.loc[mc_ketney2019_df['sample_id'].str.contains('AMY'), 'location'] = 'amygdala'
    mc_ketney2019_df.loc[mc_ketney2019_df['sample_id'].str.contains('CNC'), 'location'] = 'striatum'  # Caudate Nucleus
    mc_ketney2019_df.loc[mc_ketney2019_df['sample_id'].str.contains('CBM'), 'location'] = 'cerebellum'
    mc_ketney2019_df.loc[mc_ketney2019_df['sample_id'].str.contains('ECX'), 'location'] = 'cortex'
    mc_ketney2019_df.loc[
        mc_ketney2019_df['sample_id'].str.contains('IPL'), 'location'] = 'cortex'  # Inferior Parietal Lobule;
    mc_ketney2019_df.loc[
        mc_ketney2019_df['sample_id'].str.contains('MFG'), 'location'] = 'cortex'  # Middle Frontal Gyrus
    mc_ketney2019_df.loc[
        mc_ketney2019_df['sample_id'].str.contains('STG'), 'location'] = 'cortex'  # Superior Temporal Gyrus
    mc_ketney2019_df.loc[mc_ketney2019_df['sample_id'].str.contains('THA'), 'location'] = 'thalamus'
    mc_ketney2019_df.loc[mc_ketney2019_df['sample_id'].str.contains('VCX'), 'location'] = 'cortex'

    mc_ketney2019_df.loc[mc_ketney2019_df['sample_id'].str.contains('146'), 'condition'] = 'control'
    mc_ketney2019_df.loc[mc_ketney2019_df['sample_id'].str.contains('383'), 'condition'] = 'AD_severe'
    mc_ketney2019_df.loc[mc_ketney2019_df['sample_id'].str.contains('405'), 'condition'] = 'AD_intermediate'

    # LFQ Intensities for all proteins identified in at least one sample.

    # Case: 146 = No Tau Tangle; 383 = Severe Tangles; 405 = Intermediate Tangles.
    #    Sections: AMY=Amygdala; CNC=Caudate Nucleus  CBM=Cerebellum; ECX=Entorhinal Cortex;
    #    IPL=Inferior Parietal Lobule; MFG=Middle Frontal Gyrus;  STG=Superior Temporal Gyrus;  
    #    THA=Thalamus;  VCX=Visual Cortex;
    
    return mc_ketney2019_df


def get_hasan_2020_dataframe():
    """
    Return pandas dataframe for  Hasan 2020
    :return:
     pandas.core.frame.DataFrame: dataframe containing  Hasan 2020 data.
    """
    print("Importing  Hasan 2020 pandas dataframe.")

    # Hasan 2020 (TMT 6 brain regions mouse)

    hasan2020f = pd.ExcelFile('../data/source_data/pmic13060-sup-0002-tables1.xlsx')

    hasan2020_BS = hasan2020f.parse('BS', skiprows=1)
    hasan2020_CB = hasan2020f.parse('CB', skiprows=1)
    hasan2020_CN = hasan2020f.parse('CN', skiprows=1)
    hasan2020_FC = hasan2020f.parse('FC', skiprows=1)
    hasan2020_HP = hasan2020f.parse('HP', skiprows=1)
    hasan2020_SpC = hasan2020f.parse('SpC', skiprows=1)

    hasan2020_BS['location'] = 'brainstem'
    hasan2020_CB['location'] = 'cerebellum'
    hasan2020_CN['location'] = 'striatum'  # caudate nucleus is part of corpus striatum
    hasan2020_FC['location'] = 'cortex'
    hasan2020_HP['location'] = 'hippocampus'
    hasan2020_SpC['location'] = 'spinal cord'

    hasan2020 = pd.concat([hasan2020_BS, hasan2020_CB, hasan2020_CN, hasan2020_FC, hasan2020_HP, hasan2020_SpC
                           ], ignore_index=True, sort=False)

    hasan2020 = hasan2020.drop(['Description'], axis=1)
    hasan2020 = hasan2020.rename(columns={
        'Abundance: F1: 127, Sample, EAE1, n/a, n/a': 'Abundance_EAE1',
        'Abundance: F1: 129, Sample, n/a, EAE2, n/a': 'Abundance_EAE2',
        'Abundance: F1: 131, Sample, n/a, n/a, EAE3': 'Abundance_EAE3',
        'Abundance: F1: 126, Sample, CON1, n/a, n/a': 'Abundance_CON1',
        'Abundance: F1: 128, Sample, n/a, CON2, n/a': 'Abundance_CON2',
        'Abundance: F1: 130, Sample, n/a, n/a, CON3': 'Abundance_CON3',
        'Accession': 'Uniprot'
    })

    hasan2020_df = pd.wide_to_long(hasan2020, stubnames='Abundance',
                                   i=['Uniprot', 'location'],
                                   j='sample_id', sep='_', suffix=r'\w+')

    hasan2020_df = hasan2020_df.reset_index()

    hasan2020_df['Organism'] = 'mouse'

    hasan2020_df['Study'] = 'Hasan 2020'

    hasan2020_df['raw_data_units'] = 'tmt abundance'

    hasan2020_df = hasan2020_df.rename(columns={'Abundance': 'raw_data'})

    hasan2020_df['Age_days'] = 14 * 7 + 37 + 21  # 4 weeks + 10 weeks + 37 days #age is based on ref [12] M. Hasan, .. O.-S. Kwon, Neuroscience 2017

    hasan2020_df.loc[hasan2020_df['sample_id'].str.contains('EAE'), 'condition'] = 'EAE'
    hasan2020_df.loc[hasan2020_df['sample_id'].str.contains('CON'), 'condition'] = 'control'
    return hasan2020_df


# experimental autoimmune encephalomyelitis (EAE) mouse model - Multiple sclerosis model


def get_zhu_2018_dataframe():
    """
    Return pandas dataframe for Zhu 2018
    :return:
     pandas.core.frame.DataFrame: dataframe containing Zhu 2018 data.
    """
    print("Importing Zhu 2018 pandas dataframe.")

    zhu2018f = pd.ExcelFile('../data/source_data/136011_1_supp_135012_p8hbjd.xlsx')
    zhu2018 = zhu2018f.parse('ALl proteins with raw LFQ', skiprows=1)

    zhu2018 = zhu2018.drop(['Protein IDs', 'Protein name'], axis=1)

    zhu2018 = zhu2018.rename(columns={'Majority protein IDs': 'Uniprot',
                                      'Gene name': 'gene_names',
                                      'CTX_1': 'LFQ_CTX1',
                                      'CTX_2': 'LFQ_CTX2',
                                      'CTX_3': 'LFQ_CTX3',
                                      'CTX_4': 'LFQ_CTX4',
                                      'CC_1': 'LFQ_CC1',
                                      'CC_2': 'LFQ_CC2',
                                      'CC_3': 'LFQ_CC3',
                                      'CC_4': 'LFQ_CC4',
                                      'CP_1': 'LFQ_CP1',
                                      'CP_2': 'LFQ_CP2',
                                      'CP_3': 'LFQ_CP3',
                                      'CP_4': 'LFQ_CP4'

                                      })

    zhu2018_df = pd.wide_to_long(zhu2018, stubnames='LFQ',
                                 i=['Uniprot', 'gene_names'],
                                 j='sample_id', sep='_', suffix=r'\w+')

    zhu2018_df = zhu2018_df.reset_index()

    zhu2018_df['Organism'] = 'rat'

    zhu2018_df['Study'] = 'Zhu 2018'

    zhu2018_df['raw_data'] = 2 ** zhu2018_df['LFQ']

    zhu2018_df = zhu2018_df.drop(['LFQ'], axis=1)

    zhu2018_df['raw_data_units'] = 'LFQ'

    zhu2018_df['Age_days'] = 17 + 21  # rat P17 # Birth occurs (22nd day in rat, 19th day in mouse)
    # https://embryology.med.unsw.edu.au/embryology/index.php/Rat_Development_Stages

    zhu2018_df.loc[zhu2018_df['sample_id'].str.contains('CTX'), 'location'] = 'cortex'
    zhu2018_df.loc[zhu2018_df['sample_id'].str.contains('CC'), 'location'] = 'corpus callosum'
    zhu2018_df.loc[zhu2018_df['sample_id'].str.contains('CP'), 'location'] = 'striatum'  # caudoputamen
    return zhu2018_df


def get_kjell_2020_dataframe():
    """
    Return pandas dataframe for Kjell 2020
    :return:
     pandas.core.frame.DataFrame: dataframe containing Kjell 2020 data.
    """
    print("Importing Kjell 2020 pandas dataframe.")
    # CTX, cerebral cortex;
    # CC, corpus callosum;
    # CP,caudoputamen.

    kjell_2020LMSSf = pd.ExcelFile('../data/source_data/1-s2.0-S1934590920300023-mmc2.xlsx')
    kjell_2020LMSS = kjell_2020LMSSf.parse('LMSS data')

    kjell_2020QDSPf = pd.ExcelFile('../data/source_data/1-s2.0-S1934590920300023-mmc4.xlsx')
    kjell_2020QDSP = kjell_2020QDSPf.parse('QDSP data for table.txt')

    kjell_2020LMSS = kjell_2020LMSS[['Gene names', 'Cx - Mean abundance (log2) ',
                                     'OB - Mean abundance (log2) ', 'SEZ - Mean abundance (log2) ',
                                     'MEZ - Mean abundance (log2) ']]

    kjell_2020LMSS = kjell_2020LMSS.rename(columns={'Gene names': 'gene_names',
                                                    'Cx - Mean abundance (log2) ': 'Abund_cortexLMSS',
                                                    # Cerebral cortex
                                                    'OB - Mean abundance (log2) ': 'Abund_olfbulbLMSS',
                                                    # olfactory bulb
                                                    'SEZ - Mean abundance (log2) ': 'Abund_sezLMSS',
                                                    # largest NSC niche, the subependymal zone (SEZ)
                                                    'MEZ - Mean abundance (log2) ': 'Abund_mezLMSS'
                                                    # medial sub-ependymal zone
                                                    })

    #####
    kjell_2020LMSS = kjell_2020LMSS[~kjell_2020LMSS['gene_names'].isna()]
    kjell_2020LMSS = kjell_2020LMSS.groupby(
        'gene_names').median().reset_index()  # because duplicated genes with no info about proteins and diff values

    kjell_2020_dfLMSS = pd.wide_to_long(kjell_2020LMSS, stubnames='Abund',
                                        i=['gene_names'], j='sample_id', sep='_', suffix=r'\w+')

    kjell_2020_dfLMSS = kjell_2020_dfLMSS.reset_index()

    kjell_2020_dfLMSS.loc[kjell_2020_dfLMSS['sample_id'] == 'cortexLMSS', 'location'] = 'cortex'
    kjell_2020_dfLMSS.loc[kjell_2020_dfLMSS['sample_id'] == 'olfbulbLMSS', 'location'] = 'olfactory bulb'
    kjell_2020_dfLMSS.loc[kjell_2020_dfLMSS['sample_id'] == 'sezLMSS', 'location'] = 'subependymal zone'
    kjell_2020_dfLMSS.loc[kjell_2020_dfLMSS['sample_id'] == 'mezLMSS', 'location'] = 'medial subependymal zone'

    ###############################################################################################################

    kjell_2020QDSP = kjell_2020QDSP[['Gene names', 'Cx_FR1 (mean abundance (log2))',
                                     'Cx_FR2 (mean abundance (log2))', 'Cx_FR3 (mean abundance (log2))',
                                     'Cx_FR4 (mean abundance (log2))', 'OB_FR1 (mean abundance (log2))',
                                     'OB_FR2 (mean abundance (log2))', 'OB_FR3 (mean abundance (log2))',
                                     'OB_FR4 (mean abundance (log2))', 'SEZ_FR1 (mean abundance (log2))',
                                     'SEZ_FR2 (mean abundance (log2))', 'SEZ_FR3 (mean abundance (log2))',
                                     'SEZ_FR4 (mean abundance (log2))', 'Majority protein IDs'
                                     ]]
    kjell_2020QDSP = kjell_2020QDSP.rename(columns={'Gene names': 'gene_names',
                                                    'Majority protein IDs': 'Uniprot',
                                                    'Cx_FR1 (mean abundance (log2))': 'Abund_cortex1',
                                                    'Cx_FR2 (mean abundance (log2))': 'Abund_cortex2',
                                                    'Cx_FR3 (mean abundance (log2))': 'Abund_cortex3',
                                                    'Cx_FR4 (mean abundance (log2))': 'Abund_cortex4',
                                                    'OB_FR1 (mean abundance (log2))': 'Abund_olfbulb1',
                                                    'OB_FR2 (mean abundance (log2))': 'Abund_olfbulb2',
                                                    'OB_FR3 (mean abundance (log2))': 'Abund_olfbulb3',
                                                    'OB_FR4 (mean abundance (log2))': 'Abund_olfbulb4',
                                                    'SEZ_FR1 (mean abundance (log2))': 'Abund_sez1',
                                                    'SEZ_FR2 (mean abundance (log2))': 'Abund_sez2',
                                                    'SEZ_FR3 (mean abundance (log2))': 'Abund_sez3',
                                                    'SEZ_FR4 (mean abundance (log2))': 'Abund_sez4'
                                                    })

    kjell_2020QDSP = kjell_2020QDSP[~kjell_2020QDSP['gene_names'].isna()]
    # kjell_2020QDSP = kjell_2020QDSP.groupby(['gene_names','Uniprot']).median().reset_index()
    # because duplicated genes with no info about proteins and diff values

    kjell_2020QDSP = pd.wide_to_long(kjell_2020QDSP, stubnames='Abund',
                                     i=['gene_names', 'Uniprot'], j='sample_id', sep='_', suffix=r'\w+')

    kjell_2020QDSP = kjell_2020QDSP.reset_index()

    kjell_2020QDSP.loc[kjell_2020QDSP['sample_id'].str.contains('cortex'), 'location'] = 'cortex'
    kjell_2020QDSP.loc[kjell_2020QDSP['sample_id'].str.contains('olfbulb'), 'location'] = 'olfactory bulb'
    kjell_2020QDSP.loc[kjell_2020QDSP['sample_id'].str.contains('sez'), 'location'] = 'subependymal zone'

    kjell_2020_df = pd.concat([kjell_2020_dfLMSS, kjell_2020QDSP], ignore_index=True, sort=False)
    kjell_2020_df['Organism'] = 'mouse'
    kjell_2020_df['Study'] = 'Kjell 2020'
    kjell_2020_df['Age_days'] = 9 * 7 + 21  # 8-10 weeks
    kjell_2020_df['raw_data_units'] = 'LFQ'

    kjell_2020_df['raw_data'] = 2 ** kjell_2020_df['Abund']

    kjell_2020_df = kjell_2020_df.drop(['Abund'], axis=1)
    return kjell_2020_df
