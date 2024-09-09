import pandas as pd
import numpy as np
import networkx as nx
import itertools

original_database=pd.read_excel('E:/课题/drug_combination_database_web/upload/drug_comdb_data/drug_combination_database_check_noinfo.xlsx')
drug_target=pd.read_csv('E:/课题/drug_combination_database_web/upload/drug_comdb_data/dgidb_score_cutoff_15.txt',sep='\t')

def load_network_from_file(network_file):
    fh = open(network_file, "rb")
    G = nx.read_edgelist(fh)
    fh.close()
    return G
background_network=load_network_from_file('E:/课题/drug_combination_database_web/upload/drug_comdb_data/human_ppi_v03.tsv')

drugdb_eg=original_database.iloc[162:163,:]

def drug_target_plot(drugdb_eg):
    target_drug_name=drugdb_eg['Target_Drug'].str.split('; ',expand=True).T
    target_drug_id=drugdb_eg['Target_Drug_Drugbank_id'].str.split('; ',expand=True).T


    other_drug_df=drugdb_eg[['Other_Type_Drug','Other_Type_Drug_Drugbank_id']].dropna()
    if other_drug_df.shape[0]==0:
        other_drug_name=pd.DataFrame([np.nan])
        other_drug_id=pd.DataFrame([np.nan])
    else:
        other_drug_name=drugdb_eg['Other_Type_Drug'].str.split('; ',expand=True).T
        other_drug_id=drugdb_eg['Other_Type_Drug_Drugbank_id'].str.split('; ',expand=True).T

    other_drug_name.columns=target_drug_name.columns=['all_drug_name']
    other_drug_id.columns=target_drug_id.columns=['all_drug_id']
    other_drug_all=pd.concat([other_drug_name,other_drug_id],axis=1)
    target_drug_all=pd.concat([target_drug_name,target_drug_id],axis=1)
    all_drug_df=pd.concat([other_drug_all,target_drug_all]).drop_duplicates().dropna()

    all_drug_target_df=pd.merge(drug_target,all_drug_df,left_on='DrugBank',right_on='all_drug_id')
    drug_name_all=all_drug_target_df['all_drug_name'].drop_duplicates().tolist()
    target_all=all_drug_target_df['target'].drop_duplicates().tolist()
    
    drug_name_all_combination=list(itertools.combinations(drug_name_all, 2))
    
    drug_name_all_combination_str_all=[]
    for drug_name_all_combination_eg in drug_name_all_combination:
        drug_name_all_combination_str=['_overlap_'.join(drug_name_all_combination_eg)]
        drug_name_all_combination_str_all=drug_name_all_combination_str_all+drug_name_all_combination_str

    category_all=drug_name_all+drug_name_all_combination_str_all

    nodes_dict_all=pd.DataFrame()
    id=0
    for nodes_eg in target_all:
        nodes_dict_eg=pd.DataFrame()
        nodes_dict_eg['id']=str(id)
        id=id+1
        nodes_dict_eg['name']=nodes_eg
        nodes_dict_eg['symbolSize']


