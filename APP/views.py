from flask import Blueprint,request,send_file,Response
from flask import render_template
from flask_bootstrap import Bootstrap
from APP.models import DrugComdb,Match_VCF,db
from sqlalchemy import Table
import pandas as pd
import csv
import io
from datetime import datetime
import math
from itertools import combinations
import json
import numpy as np
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import subprocess
import functools
from multiprocessing import Pool
from flask import redirect, url_for
import warnings
warnings.filterwarnings("ignore")

blue=Blueprint('web_index',__name__)

original_database=pd.read_excel('/h/tianyi/oncodrug/upload/drug_comdb_data/drug_combination_database.xlsx')
drug_target=pd.read_csv('/h/tianyi/oncodrug/upload/drug_comdb_data/dgidb_score_cutoff_15.txt',sep='\t')
drug_info_detail_all=pd.read_excel('/h/tianyi/oncodrug/upload/drug_comdb_data/each_drug_info_indb.xlsx')
drugbank_id_raw_data=pd.read_csv('/h/tianyi/oncodrug/upload/drug_comdb_data/drugbank_vovabulary_extract_all.txt',sep='\t')


def transvar_translation(mutation_id,mut_type):
    if mut_type=='hgvsp':
        cmd_line=f'transvar panno -i \'{mutation_id}\' --ensembl '
    elif mut_type=='hgvsc':
        cmd_line=f'transvar canno --ccds -i \'{mutation_id}\' '
    elif mut_type=='hgvsg':
        cmd_line=f'transvar ganno --ccds -i \'{mutation_id}\' '
    process_completed = subprocess.run(
    cmd_line, shell =True, encoding='utf-8', stdout = subprocess.PIPE,
    stderr = subprocess.PIPE
    )
    process_completed_pro=str(process_completed).split('stdout=')[1].split('\', ')[0].replace('\'','')
    process_completed_pro1=process_completed_pro.split('\\n')
    process_completed_all=pd.DataFrame()
    for process_completed_eg in process_completed_pro1:
        process_completed_process_eg=pd.DataFrame(process_completed_eg.split('\\t')).T
        process_completed_all=pd.concat([process_completed_process_eg,process_completed_all])
    process_completed_all=process_completed_all.dropna()
    process_completed_all=process_completed_all[process_completed_all[0]!='input']
    vcf_input=pd.DataFrame()
    vcf_input['Action_mut_gene']=process_completed_all[2].tolist()
    vcf_input['Action_mut_event']=process_completed_all[4].str.split('/').str[2].str.replace('p.','').tolist()
    return(vcf_input)
            


def save_to_file(data,filename):
    # 将数据保存到文件
    with open(f'upload/{filename}', 'w') as file:
        file.write(data)

def create_csv(data):
    # 创建CSV文件并写入数据
    current_time = datetime.now()
    formatted_time = current_time.strftime("%Y%m%d_%H%M%S")
    filename = f'upload/{formatted_time}_users.txt'

    with open(filename, 'w') as file:
        file.write('Drug_combination_ID\tEvidence_level\tPrioritization_score\tTargeted_drug\tNon_targeted_drug\tCancer_type_oncotree_Level2\tBiomarkers\tResponse\tAdverse_effect\tDrug_dosage\n')

        for item in data:
            file.write(f"{item['Drug_combination_ID']}\t{item['Evidence_level']}\t{item['Prioritization_score']}\t{item['Targeted_drug']}\t{item['Non_targeted_drug']}\t{item['Cancer_type_oncotree_Level2']}\t{item['Biomarkers']}\t{item['Response']}\t{item['Adverse_effect']}\t{item['Drug_dosage']}\n")  # 写入数据行
    return filename

#####数据超链接
def x_y_drug_pos(drug_name_all):
    if drug_name_all.shape[0]==2:
        x_all=[300,800]
        y_all=[300,300]
    elif drug_name_all.shape[0]==3:
        x_all=[300,800,550]
        y_all=[300,300,500]
    elif drug_name_all.shape[0]==4:
        x_all=[300,800,550,550]
        y_all=[300,300,100,500]
    elif drug_name_all.shape[0]==5:
        x_all=[300,400,450,350,250]
        y_all=[100,300,400,500,250]
    elif drug_name_all.shape[0]==6:
        x_all=[300,400,250,300,400,450]
        y_all=[300,300,400,500,500,400]
    elif drug_name_all.shape[0]==6:
        x_all=[300,400,450,400,300,250,350]
        y_all=[300,300,400,500,500,400,400]
    return(x_all,y_all)


def single_drug_pos(node_count, radius=300):
    if node_count <= 0:
        return [0,0]

    # 计算每个节点之间的角度间隔
    angle_interval = 2 * math.pi / node_count

    # 生成节点坐标
    coordinates = []
    for i in range(node_count):
        angle = i * angle_interval
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        coordinates.append((x, y))
    coordinates=pd.DataFrame(coordinates)
    
    return coordinates[0].tolist(),coordinates[1].tolist()

@blue.route('/oncodrug')
def home():
    ###
    return render_template('navbar.html')

@blue.route('/oncodrug/statistics')
def statistics():
    ###
    return render_template('statistics.html')

@blue.route('/oncodrug/biomarker')
def biomarker():
    ###
    return render_template('biomarker.html')


@blue.route('/oncodrug/navbar_error_gene')
def navbar_error_gene():
    ###
    return render_template('navbar_error_gene.html')

@blue.route('/oncodrug/discovery_drugcomb',methods=['POST'])
def discovery_drugcomb():
    try:
        searchtype = request.form['inlineRadioOptions']
        if searchtype=='Drug: e.g. Afatinib':
            searchinput = request.form['query']
            searchcancer_type = request.form['searchcancer_type']
            all_combination_pairs_part1=original_database[original_database['Targeted drug'].str.contains(searchinput,case=False).fillna(False)]
            all_combination_pairs_part2=original_database[original_database['Non-targeted drug'].str.contains(searchinput,case=False).fillna(False)]
            all_combination_pairs=pd.concat([all_combination_pairs_part1,all_combination_pairs_part2])
            if all_combination_pairs.shape[0]!=0:
                return redirect(url_for('web_index.each_drug_cancer_search_info', drugname=searchinput,cancer_type=searchcancer_type))
            else:
                return render_template('navbar_error.html' )
        elif searchtype=="Gene: e.g. BRAF or \'BRAF; KRAS\'":
            searchinput = request.form['query']
            searchcancer_type = request.form['searchcancer_type']
            searchinput=searchinput.replace("\'",'')
            if ';' in searchinput:
                searchinput_split=searchinput.replace(" ",'').split(';')
            else:
                searchinput_split=searchinput.split(' ')
            action_mutation_data=original_database[original_database['Action mutations'].astype(str).apply(lambda x: all(gene in x for gene in searchinput_split))]
            if action_mutation_data.shape[0]!=0:
                return redirect(url_for('web_index.each_action_mutation_cancer_search_info', actmut=searchinput,cancer_type=searchcancer_type))
            else:
                return render_template('navbar_error_gene.html' )
        elif searchtype=="Biomarker: e.g. BRAF V600E or \'TP53 V197G; KRAS G12D\'":
            searchinput = request.form['query']
            searchcancer_type = request.form['searchcancer_type']
            if '; ' in searchinput:
                searchinput_split=searchinput.split('; ')
            else:
                searchinput_split=searchinput.split(';')
            action_mutation_data=original_database[original_database['Action mutations'].astype(str).apply(lambda x: all(gene in x for gene in searchinput_split))]
            if action_mutation_data.shape[0]!=0:
                return redirect(url_for('web_index.each_biomarker_cancer_search_info', actmut=searchinput,cancer_type=searchcancer_type))
            else:
                return render_template('navbar_error_biomarker.html' )
    except:
      return render_template('navbar_error.html' )

       


@blue.route('/oncodrug/each_drug_cancer_search_info/<drugname>_<cancer_type>')
def each_drug_cancer_search_info(drugname,cancer_type):
    try:
        ####对应models中person()
            ####对应models中person()
        all_combination_pairs_part1=original_database[original_database['Targeted drug'].str.contains(drugname,case=False).fillna(False)]
        all_combination_pairs_part2=original_database[original_database['Non-targeted drug'].str.contains(drugname,case=False).fillna(False)]
        all_combination_pairs=pd.concat([all_combination_pairs_part1,all_combination_pairs_part2])
        if cancer_type!='All':
            all_combination_pairs=all_combination_pairs[all_combination_pairs['Cancer type (oncotree Level2)']==cancer_type]
        else:
            all_combination_pairs=all_combination_pairs
        all_combination_pairs.columns=pd.Series(all_combination_pairs.columns.tolist()).str.replace(' ','_').str.replace('-','_').str.replace('(','').str.replace(')','')

        drug_name_all=np.union1d(all_combination_pairs['Targeted_drug'].astype(str).tolist(),all_combination_pairs['Non_targeted_drug'].astype(str).tolist())
        drug_name_all=pd.Series(drug_name_all).dropna()

        drug_name_all_str=[]
        for drug_name_eg in drug_name_all:
            drug_name_eg_str=drug_name_eg.split('; ')
            drug_name_all_str=drug_name_eg_str+drug_name_all_str

        drug_name_all=pd.DataFrame(drug_name_all_str).drop_duplicates()
        drug_name_all.columns=['drug_name']
        drug_name_all['names']='name'
        drug_name_all=drug_name_all[(drug_name_all['drug_name']!='nan') & (drug_name_all['drug_name']!='No info')]
        ##drugbank id###
        all_drugbank_info_all=pd.merge(drug_name_all,drugbank_id_raw_data,left_on=drug_name_all['drug_name'].str.upper(),right_on=drugbank_id_raw_data['All_Drug'].str.upper(),how='left')[['drug_name','DrugBank']]
        all_drugbank_info_all.columns=['drugbank_name','drugbank_id']
        all_drugbank_info_all=all_drugbank_info_all.dropna().drop_duplicates(['drugbank_id'])

        link_combi=all_drugbank_info_all[['drugbank_name']]
        link_combi['names']=drugname
        link_combi.columns=['node1','node2']
        link_combi=link_combi[link_combi['node1'] != link_combi['node2']]

    ####核心节点
        drugbank_core_data=all_drugbank_info_all[all_drugbank_info_all['drugbank_name']==drugname]
        if drugbank_core_data.shape[0]!=0:
            core_drug_name_eg=''.join(drugbank_core_data['drugbank_name'])
            core_drugbank_eg=''.join(drugbank_core_data['drugbank_id'])
        else:
            core_drug_name_eg=drugname
            core_drugbank_eg=''
        core_node_data={
                'name': f'{core_drug_name_eg}',
                'x': 0,
                'y': 0,'symbolSize': [25, 25],'itemStyle': { 'color': '#67C1C4' } ,
                'tooltip': f'{core_drug_name_eg}: {core_drugbank_eg}'
                }
    ###其他节点
        #####nodes position####
        if len(all_drugbank_info_all['drugbank_name'].drop_duplicates())!=1:
            all_comb_drug_name_drop=np.setdiff1d(all_drugbank_info_all['drugbank_name'],drugname)
        else:
            all_comb_drug_name_drop=drugname
        x_all,y_all=single_drug_pos(len(all_comb_drug_name_drop))
        node_data_nocore=[]
        for pos in range(len(all_comb_drug_name_drop)):
            drugname_eg=all_comb_drug_name_drop[pos]
            all_drugbank_info_eg=all_drugbank_info_all[all_drugbank_info_all['drugbank_name']==drugname_eg]
            if all_drugbank_info_eg.shape[0]==0:
                drug_name_eg=drugname_eg
                drugbank_eg='No info'
            else:
                drug_name_eg=''.join(all_drugbank_info_eg['drugbank_name'])
                drugbank_eg=''.join(all_drugbank_info_eg['drugbank_id'])
            x_eg=x_all[pos]
            y_eg=y_all[pos]
            node_data_eg={
                'name': f'{drug_name_eg}',
                'x': x_eg,
                'y': y_eg,
                'tooltip': f'{drug_name_eg}: {drugbank_eg}'
                }
            node_data_nocore=node_data_nocore+[node_data_eg]

        node_data=[core_node_data]+node_data_nocore
        link_data=[]
        for link_pos in range(link_combi.shape[0]):
            source_eg=''.join(link_combi.iloc[link_pos:(link_pos+1),:]['node1'])
            target_eg=''.join(link_combi.iloc[link_pos:(link_pos+1),:]['node2'])
            link_data_eg={
            'source': f'{source_eg}',
            'target': f'{target_eg}',
            'label': {
            'show': 'true','formatter': ''
            },
            'lineStyle': {
            'curveness': 0,
            }
            }
            link_data=link_data+[link_data_eg]
        
        option = {
                'title': {
                'text': 'Drug combination therapy regimens','left': 'center','top': 25, 
                },
                    'tooltip': {},
                    'animationDurationUpdate': 1500,
                    'animationEasingUpdate': 'quinticInOut',
                    'series': [
                    {
                        'type': 'graph',
                        'layout': 'none',
                        'symbolSize': [25/2, 25/2],
                        'roam': 'true',
                        'label': {
                        'show': 'true','formatter': '{b}'
                        },
                        'edgeLabel': {
                        'fontSize': 20
                        },
                        'itemStyle': {
                        'borderColor': '#000',
                        'borderWidth': 1,
                        'borderType': 'solid',
                        'shadowBlur': 10,
                        'shadowColor': 'rgba(0, 0, 0, 0.3)'},
                        
                        'data': node_data,
                        'links': link_data,
                        'lineStyle': {
                        'opacity': 0.9,
                        'width': 2,'smooth': 'false'
                        }
                    }
                    ]
                    }
        # 将 option 对象传递给模板渲染
        option_json = json.dumps(option)
        if core_drugbank_eg!='':
            drugbank_link=f'https://go.drugbank.com/drugs/{core_drugbank_eg}'
        else:
            drugbank_link=''
        #####all drug info detail#####
        drug_info_detail_eg=drug_info_detail_all[drug_info_detail_all['drugbank_id']==core_drugbank_eg]
        drug_info_detail_description1=drug_info_detail_eg['description'].dropna()
        if len(drug_info_detail_description1)==0:
            drug_info_detail_description2=''
        else:
            drug_info_detail_description2=''.join(drug_info_detail_eg['description'])

        drug_info_detail_target1=drug_info_detail_eg['target'].dropna()
        if len(drug_info_detail_target1)==0:
            drug_info_detail_target2=''
        else:
            drug_info_detail_target2=''.join(drug_info_detail_eg['target'])
        
        drug_info_detail_moa1=drug_info_detail_eg['moa'].dropna()
        if len(drug_info_detail_moa1)==0:
            drug_info_detail_moa2=''
        else:
            drug_info_detail_moa2=''.join(drug_info_detail_eg['moa'])

        drug_info_detail_pharmacodynamics1=drug_info_detail_eg['pharmacodynamics'].dropna()
        if len(drug_info_detail_pharmacodynamics1)==0:
            drug_info_detail_pharmacodynamics2=''
        else:
            drug_info_detail_pharmacodynamics2=''.join(drug_info_detail_eg['pharmacodynamics'])

        drug_info_detail_indication1=drug_info_detail_eg['indication'].dropna()
        if len(drug_info_detail_indication1)==0:
            drug_info_detail_indication2=''
        else:
            drug_info_detail_indication2=''.join(drug_info_detail_eg['indication'])
            
        all_cancer_type_id_all=pd.DataFrame(all_combination_pairs['Cancer_type_oncotree_Level2'].value_counts()).iloc[0:10,:]
        cancer_type_link_data=[]
        for pos in range(all_cancer_type_id_all.shape[0]):
            all_cancer_type_id_all_eg=all_cancer_type_id_all.iloc[pos:(pos+1),:]
            cancer_type_eg=str(all_cancer_type_id_all_eg.index.tolist()[0])
            count_eg=int(all_cancer_type_id_all_eg['count'])
            cancer_type_link_data_eg={
                'name':cancer_type_eg,
                'value':count_eg
            }
            cancer_type_link_data=cancer_type_link_data+[cancer_type_link_data_eg]

            cancer_type_option_json={
                'title': {
                        'text': 'Cancer types','left': 'center','top': 25, 
                        }
                ,'tooltip': {},
        'series': [ {
        'type': 'wordCloud',
        'gridSize': 2,
        'sizeRange': [10, 30],
        'rotationRange': [-90, 90],
        'shape': 'apple',
        'width': 600,
        'height': 400,
        'drawOutOfBound': 'true',
        'textStyle': {
        'fontFamily': 'sans-serif',
        'fontWeight': 'bold',
            'emphasis': {
            'shadowBlur': 10,
            'shadowColor': '#333'
        }
        },'data':cancer_type_link_data}]
                                };

        cancer_type_option_json = json.dumps(cancer_type_option_json)
            
        all_evidence_id_all=pd.DataFrame(all_combination_pairs['Evidence_level'].value_counts())
        all_evidence_id_all.columns=['count']
        all_evidence_id_all['Evidence_level']=all_evidence_id_all.index.tolist()
        all_evidence_id_all.reset_index(drop=True,inplace=True)
        evidence_link_data=[]
        for evidence_level in ['Level A','Level B','Level C','Level D']:
            all_evidence_id_all_eg=all_evidence_id_all[all_evidence_id_all['Evidence_level']==evidence_level]
            if all_evidence_id_all_eg.shape[0]!=0:
                evidence_eg=evidence_level
                count_eg=int(all_evidence_id_all_eg['count'])
                evidence_link_data_eg={
                    'name':evidence_eg,
                    'value':count_eg
                }
            else:
                evidence_link_data_eg={
                    'name':evidence_level,
                    'value':0
                }
            evidence_link_data=evidence_link_data+[evidence_link_data_eg]
        evidence_option_json={
        'title': {
        'text': 'Statistical evidence levels for data entries','left': 'center',
        'top': 25,    
        'left': 'center'
        },
            'tooltip': {
            'trigger': 'item'
            },
            'legend': {
            'bottom': '12.5%',
        'left': 'center'
            },
            'series': [
            {
                'name': 'Access From',
                'type': 'pie',
                'radius': ['30%', '60%'],
                'avoidLabelOverlap': 'false',
                'itemStyle': {
                'borderRadius': 10,
                'borderColor': '#fff',
                'borderWidth': 2
                },
                'label': {
                'show': 'false',
                'position': 'center'
                },
                'emphasis': {
                'label': {
                    'show': 'true',
                    'fontSize': 20,
                    'fontWeight': 'bold'
                }
                },
                'labelLine': {
                'show': 'false'
                },
                'data': evidence_link_data
            }
            ]
        };

        evidence_option_json = json.dumps(evidence_option_json)
        
        all_mutation_id_all=pd.DataFrame(all_combination_pairs['Action_mutations'].value_counts())
        all_mutation_id_all=all_mutation_id_all.iloc[0:10,:]
        all_mutation_id_all=all_mutation_id_all.loc[list(~(pd.Series(all_mutation_id_all.index.tolist()).str.contains('Unavailable|No info',case=False))),:]
        mutation_link_data=[]
        for pos in range(all_mutation_id_all.shape[0]):
            all_mutation_id_all_eg=all_mutation_id_all.iloc[pos:(pos+1),:]
            mutation_eg=str(all_mutation_id_all_eg.index.tolist()[0])
            count_eg=int(all_mutation_id_all_eg['count'])
            mutation_link_data_eg={
                'name':mutation_eg,
                'value':count_eg
            }
            mutation_link_data=mutation_link_data+[mutation_link_data_eg]

            mutation_option_json={
                'title': {
                        'text': 'Biomarkers predicting drug responses','left': 'center','top': 25,
                        }
                ,'tooltip': {},
        'series': [ {
        'type': 'wordCloud',
        'gridSize': 2,
        'sizeRange': [10, 30],
        'rotationRange': [-90, 90],
        'shape': 'apple',
        'width': 600,
        'height': 400,
        'drawOutOfBound': 'true',
        'textStyle': {
        'fontFamily': 'sans-serif',
        'fontWeight': 'bold',
            'emphasis': {
            'shadowBlur': 10,
            'shadowColor': '#333'
        }
        },'data':mutation_link_data}]
                                };

        mutation_option_json = json.dumps(mutation_option_json)
        return render_template('each_drug_cancer_search_info.html',core_drugbank_eg=core_drugbank_eg,drugbank_link=drugbank_link,drugname=drugname,drug_info_detail_description=drug_info_detail_description2,drug_info_detail_target=drug_info_detail_target2,drug_info_detail_moa=drug_info_detail_moa2,drug_info_detail_pharmacodynamics=drug_info_detail_pharmacodynamics2,drug_info_detail_indication=drug_info_detail_indication2,all_combination_pairs=all_combination_pairs.to_dict(orient='records'),option_json=option_json,cancer_type_option_json=cancer_type_option_json,evidence_option_json=evidence_option_json,mutation_option_json=mutation_option_json)
    except:
        return render_template('navbar_error.html')

@blue.route('/oncodrug/each_action_mutation_cancer_search_info/<actmut>_<cancer_type>')
def each_action_mutation_cancer_search_info(actmut,cancer_type):
    ####对应models中person()
        ####对应models中person()
    try:
        actmut=actmut.replace("\'",'')
        if ';' in actmut:
            actmut_split=actmut.replace(" ",'').split(';')
        else:
            actmut_split=actmut.split(' ')
        action_mutation_data=original_database[original_database['Action mutations'].astype(str).apply(lambda x: all(gene in x for gene in actmut_split))]
        
        if cancer_type!='All':
            action_mutation_data=action_mutation_data[action_mutation_data['Cancer type (oncotree Level2)']==cancer_type]
        else:
            action_mutation_data=action_mutation_data
        
        action_mutation_data.columns=pd.Series(action_mutation_data.columns.tolist()).str.replace(' ','_').str.replace('-','_').str.replace('(','').str.replace(')','')

        all_action_mutation_id_all=[]
        for action_mut_eg in action_mutation_data['Action_mutations'].drop_duplicates():
            action_mut_eg_split=action_mut_eg.split('; ')
            all_action_mutation_id_all=np.union1d(action_mut_eg_split,all_action_mutation_id_all)
        all_action_mutation_id_all=pd.Series(all_action_mutation_id_all)
        all_action_mutation_id_all=all_action_mutation_id_all[~(all_action_mutation_id_all.str.contains(actmut,case=False)) ].tolist()
        ##drugbank id###
        link_act_combi=pd.DataFrame()
        link_act_combi['action_mut1']=all_action_mutation_id_all
        link_act_combi['action_mut2']=actmut
        link_act_combi=link_act_combi[link_act_combi['action_mut1'] != link_act_combi['action_mut2']]

    ####核心节点

        core_node_data={
                'name': f'{actmut}',
                'x': 0,
                'y': 0,'itemStyle': { 'color': '#F4B2B0' } ,'symbolSize': [25, 25]
                }
    ###其他节点
        #####nodes position####
        if len(all_action_mutation_id_all)==1:
            all_comb_mut_name_drop=all_action_mutation_id_all
        else:
            all_comb_mut_name_drop=np.setdiff1d(all_action_mutation_id_all,actmut)
        
        x_all,y_all=single_drug_pos(len(all_comb_mut_name_drop))
        node_data_nocore=[]
        for pos in range(len(all_comb_mut_name_drop)):
            action_mut_eg=all_comb_mut_name_drop[pos]
            x_eg=x_all[pos]
            y_eg=y_all[pos]
            node_data_eg={
                'name': f'{action_mut_eg}',
                'x': x_eg,
                'y': y_eg
                }
            node_data_nocore=node_data_nocore+[node_data_eg]
        node_data=[core_node_data]+node_data_nocore

        link_data=[]
        for link_pos in range(link_act_combi.shape[0]):
            source_eg=''.join(link_act_combi.iloc[link_pos:(link_pos+1),:]['action_mut1'])
            target_eg=''.join(link_act_combi.iloc[link_pos:(link_pos+1),:]['action_mut2'])
            link_data_eg={
            'source': f'{source_eg}',
            'target': f'{target_eg}',
            'label': {
            'show': 'true','formatter': ''
            },
            'lineStyle': {
            'curveness': 0,
            }
            }
            link_data=link_data+[link_data_eg]
        
        option = {            'title': {
                        'text': 'Biomarkers predicting drug responses','left': 'center','top':25
                        },
                    'animationDurationUpdate': 1500,
                    'animationEasingUpdate': 'quinticInOut',
                    'series': [
                    {
                        'type': 'graph',
                        'layout': 'none',
                        'symbolSize': [25/2, 25/2],
                        'roam': 'true',
                        'label': {
                        'show': 'true','formatter': '{b}'
                        },
                        'edgeLabel': {
                        'fontSize': 20
                        },
                        'itemStyle': {
                        'borderColor': '#000',
                        'borderWidth': 1,
                        'borderType': 'solid',
                        'shadowBlur': 10,
                        'shadowColor': 'rgba(0, 0, 0, 0.3)'},
                        
                        'data': node_data,
                        'links': link_data,
                        'lineStyle': {
                        'opacity': 0.9,
                        'width': 2,'smooth': 'false'
                        },'tooltip': {'show': 'false'},
                    }
                    ]
                    };
        # 将 option 对象传递给模板渲染
        option_json = json.dumps(option)
        
        
        ##drug name####
        comb_drug_all=pd.DataFrame()
        for pos in range(action_mutation_data.shape[0]):
            action_mutation_data_eg=action_mutation_data.iloc[pos:(pos+1),:]
            drug_name_eg=action_mutation_data_eg['Targeted_drug'].str.split('; ').iloc[0]+action_mutation_data_eg['Targeted_drug'].str.split('; ').iloc[0]
            drug_name_eg=np.setdiff1d(drug_name_eg,['No info'])
            drug_name_sort_eg=pd.Series(drug_name_eg).drop_duplicates().sort_values()
            drug_name_sort_str_eg='; '.join(drug_name_sort_eg)
            drug_name_sort_str_df_eg=pd.DataFrame([drug_name_sort_str_eg])
            comb_drug_all=pd.concat([drug_name_sort_str_df_eg,comb_drug_all])
        
        all_drug_id_all=pd.DataFrame(comb_drug_all[0].value_counts()).sort_values(0).iloc[0:10,:]

        drug_link_data=[]
        for pos in range(all_drug_id_all.shape[0]):
            all_drug_id_all_eg=all_drug_id_all.iloc[pos:(pos+1),:]
            all_drug_id_all_eg.columns=['count']
            drug_eg=str(all_drug_id_all_eg.index.tolist()[0])
            count_eg=int(all_drug_id_all_eg['count'])
            drug_link_data_eg={
                'name':drug_eg,
                'value':count_eg
            }
            drug_link_data=drug_link_data+[drug_link_data_eg]

            drug_option_json={
                'title': {
                        'text': 'Drugs associated with biomarkers','left': 'center','top':25
                        }
                ,'tooltip': {},
        'series': [ {
        'type': 'wordCloud',
        'gridSize': 2,
        'sizeRange': [10, 30],
        'rotationRange': [-90, 90],
        'shape': 'apple',
        'width': 600,
        'height': 400,
        'drawOutOfBound': 'true',
        'textStyle': {
        'fontFamily': 'sans-serif',
        'fontWeight': 'bold',
            'emphasis': {
            'shadowBlur': 10,
            'shadowColor': '#333'
        }
        },'data':drug_link_data}]
                                };

        drug_option_json = json.dumps(drug_option_json)
        
        all_evidence_id_all=pd.DataFrame(action_mutation_data['Evidence_level'].value_counts())
        all_evidence_id_all.columns=['count']
        all_evidence_id_all['Evidence_level']=all_evidence_id_all.index.tolist()
        all_evidence_id_all.reset_index(drop=True,inplace=True)
        evidence_link_data=[]
        for evidence_level in ['Level A','Level B','Level C','Level D']:
            all_evidence_id_all_eg=all_evidence_id_all[all_evidence_id_all['Evidence_level']==evidence_level]
            if all_evidence_id_all_eg.shape[0]!=0:
                evidence_eg=evidence_level
                count_eg=int(all_evidence_id_all_eg['count'])
                evidence_link_data_eg={
                    'name':evidence_eg,
                    'value':count_eg
                }
            else:
                evidence_link_data_eg={
                    'name':evidence_level,
                    'value':0
                }
            evidence_link_data=evidence_link_data+[evidence_link_data_eg]
        evidence_option_json={
        'title': {
        'text': 'Statistical evidence levels for data entries','left': 'center',
        'top': 25,    
        'left': 'center'
        },
            'tooltip': {
            'trigger': 'item'
            },
            'legend': {
            'bottom': '12.5%',
        'left': 'center'
            },
            'series': [
            {
                'name': 'Access From',
                'type': 'pie',
                'radius': ['30%', '60%'],
                'avoidLabelOverlap': 'false',
                'itemStyle': {
                'borderRadius': 10,
                'borderColor': '#fff',
                'borderWidth': 2
                },
                'label': {
                'show': 'false',
                'position': 'center'
                },
                'emphasis': {
                'label': {
                    'show': 'true',
                    'fontSize': 20,
                    'fontWeight': 'bold'
                }
                },
                'labelLine': {
                'show': 'false'
                },
                'data': evidence_link_data
            }
            ]
        };

        evidence_option_json = json.dumps(evidence_option_json)
        
        all_cancer_type_id_all=pd.DataFrame(action_mutation_data['Cancer_type_oncotree_Level2'].value_counts())
        all_cancer_type_id_all=all_cancer_type_id_all.iloc[0:10,:]
        cancer_type_link_data=[]
        for pos in range(all_cancer_type_id_all.shape[0]):
            all_cancer_type_id_all_eg=all_cancer_type_id_all.iloc[pos:(pos+1),:]
            cancer_type_eg=str(all_cancer_type_id_all_eg.index.tolist()[0])
            count_eg=int(all_cancer_type_id_all_eg['count'])
            cancer_type_link_data_eg={
                'name':cancer_type_eg,
                'value':count_eg
            }
            cancer_type_link_data=cancer_type_link_data+[cancer_type_link_data_eg]

            cancer_type_option_json={
                'title': {
                        'text': 'Cancer types associated with biomarkers','left': 'center','top':25
                        }
                ,'tooltip': {},
        'series': [ {
        'type': 'wordCloud',
        'gridSize': 2,
        'sizeRange': [10, 30],
        'rotationRange': [-90, 90],
        'shape': 'apple',
        'width': 600,
        'height': 400,
        'drawOutOfBound': 'true',
        'textStyle': {
        'fontFamily': 'sans-serif',
        'fontWeight': 'bold',
            'emphasis': {
            'shadowBlur': 10,
            'shadowColor': '#333'
        }
        },'data':cancer_type_link_data}]
                                };

        cancer_type_option_json = json.dumps(cancer_type_option_json)
        
        return render_template('each_action_mutation_cancer_search_info.html',actmut=actmut,action_mutation_data=action_mutation_data.to_dict(orient='records'),option_json=option_json,drug_option_json=drug_option_json,evidence_option_json=evidence_option_json,cancer_type_option_json=cancer_type_option_json)
    except:
        return render_template('navbar_error_gene.html')


@blue.route('/oncodrug/each_biomarker_cancer_search_info/<actmut>_<cancer_type>')
def each_biomarker_cancer_search_info(actmut,cancer_type):
    ####对应models中person()
        ####对应models中person()
    try:
        if '; ' in actmut:
            actmut_split=actmut.split('; ')
        else:
            actmut_split=actmut.split(';')
        action_mutation_data=original_database[original_database['Action mutations'].astype(str).apply(lambda x: all(gene in x for gene in actmut_split))]
        
        if cancer_type!='All':
            action_mutation_data=action_mutation_data[action_mutation_data['Cancer type (oncotree Level2)']==cancer_type]
        else:
            action_mutation_data=action_mutation_data
        
        action_mutation_data.columns=pd.Series(action_mutation_data.columns.tolist()).str.replace(' ','_').str.replace('-','_').str.replace('(','').str.replace(')','')

        all_action_mutation_id_all=[]
        for action_mut_eg in action_mutation_data['Action_mutations'].drop_duplicates():
            action_mut_eg_split=action_mut_eg.split('; ')
            all_action_mutation_id_all=np.union1d(action_mut_eg_split,all_action_mutation_id_all)
        all_action_mutation_id_all=pd.Series(all_action_mutation_id_all)
        all_action_mutation_id_all=all_action_mutation_id_all[~(all_action_mutation_id_all.str.contains(actmut,case=False)) ].tolist()
        ##drugbank id###
        link_act_combi=pd.DataFrame()
        link_act_combi['action_mut1']=all_action_mutation_id_all
        link_act_combi['action_mut2']=actmut
        link_act_combi=link_act_combi[link_act_combi['action_mut1'] != link_act_combi['action_mut2']]

    ####核心节点

        core_node_data={
                'name': f'{actmut}',
                'x': 0,
                'y': 0,'itemStyle': { 'color': '#F4B2B0' } ,'symbolSize': [25, 25]
                }
    ###其他节点
        #####nodes position####
        if len(all_action_mutation_id_all)==1:
            all_comb_mut_name_drop=all_action_mutation_id_all
        else:
            all_comb_mut_name_drop=np.setdiff1d(all_action_mutation_id_all,actmut)
        
        x_all,y_all=single_drug_pos(len(all_comb_mut_name_drop))
        node_data_nocore=[]
        for pos in range(len(all_comb_mut_name_drop)):
            action_mut_eg=all_comb_mut_name_drop[pos]
            x_eg=x_all[pos]
            y_eg=y_all[pos]
            node_data_eg={
                'name': f'{action_mut_eg}',
                'x': x_eg,
                'y': y_eg
                }
            node_data_nocore=node_data_nocore+[node_data_eg]
        node_data=[core_node_data]+node_data_nocore

        link_data=[]
        for link_pos in range(link_act_combi.shape[0]):
            source_eg=''.join(link_act_combi.iloc[link_pos:(link_pos+1),:]['action_mut1'])
            target_eg=''.join(link_act_combi.iloc[link_pos:(link_pos+1),:]['action_mut2'])
            link_data_eg={
            'source': f'{source_eg}',
            'target': f'{target_eg}',
            'label': {
            'show': 'true','formatter': ''
            },
            'lineStyle': {
            'curveness': 0,
            }
            }
            link_data=link_data+[link_data_eg]
        
        option = {            'title': {
                        'text': 'Biomarkers predicting drug responses','left': 'center','top':25
                        },
                    'animationDurationUpdate': 1500,
                    'animationEasingUpdate': 'quinticInOut',
                    'series': [
                    {
                        'type': 'graph',
                        'layout': 'none',
                        'symbolSize': [25/2, 25/2],
                        'roam': 'true',
                        'label': {
                        'show': 'true','formatter': '{b}'
                        },
                        'edgeLabel': {
                        'fontSize': 20
                        },
                        'itemStyle': {
                        'borderColor': '#000',
                        'borderWidth': 1,
                        'borderType': 'solid',
                        'shadowBlur': 10,
                        'shadowColor': 'rgba(0, 0, 0, 0.3)'},
                        
                        'data': node_data,
                        'links': link_data,
                        'lineStyle': {
                        'opacity': 0.9,
                        'width': 2,'smooth': 'false'
                        },'tooltip': {'show': 'false'},
                    }
                    ]
                    };
        # 将 option 对象传递给模板渲染
        option_json = json.dumps(option)
        
        
        ##drug name####
        comb_drug_all=pd.DataFrame()
        for pos in range(action_mutation_data.shape[0]):
            action_mutation_data_eg=action_mutation_data.iloc[pos:(pos+1),:]
            drug_name_eg=action_mutation_data_eg['Targeted_drug'].str.split('; ').iloc[0]+action_mutation_data_eg['Targeted_drug'].str.split('; ').iloc[0]
            drug_name_eg=np.setdiff1d(drug_name_eg,['No info'])
            drug_name_sort_eg=pd.Series(drug_name_eg).drop_duplicates().sort_values()
            drug_name_sort_str_eg='; '.join(drug_name_sort_eg)
            drug_name_sort_str_df_eg=pd.DataFrame([drug_name_sort_str_eg])
            comb_drug_all=pd.concat([drug_name_sort_str_df_eg,comb_drug_all])
        
        all_drug_id_all=pd.DataFrame(comb_drug_all[0].value_counts()).sort_values(0).iloc[0:10,:]

        drug_link_data=[]
        for pos in range(all_drug_id_all.shape[0]):
            all_drug_id_all_eg=all_drug_id_all.iloc[pos:(pos+1),:]
            all_drug_id_all_eg.columns=['count']
            drug_eg=str(all_drug_id_all_eg.index.tolist()[0])
            count_eg=int(all_drug_id_all_eg['count'])
            drug_link_data_eg={
                'name':drug_eg,
                'value':count_eg
            }
            drug_link_data=drug_link_data+[drug_link_data_eg]

            drug_option_json={
                'title': {
                        'text': 'Drugs associated with biomarkers','left': 'center','top':25
                        }
                ,'tooltip': {},
        'series': [ {
        'type': 'wordCloud',
        'gridSize': 2,
        'sizeRange': [10, 30],
        'rotationRange': [-90, 90],
        'shape': 'apple',
        'width': 600,
        'height': 400,
        'drawOutOfBound': 'true',
        'textStyle': {
        'fontFamily': 'sans-serif',
        'fontWeight': 'bold',
            'emphasis': {
            'shadowBlur': 10,
            'shadowColor': '#333'
        }
        },'data':drug_link_data}]
                                };

        drug_option_json = json.dumps(drug_option_json)
        
        all_evidence_id_all=pd.DataFrame(action_mutation_data['Evidence_level'].value_counts())
        all_evidence_id_all.columns=['count']
        all_evidence_id_all['Evidence_level']=all_evidence_id_all.index.tolist()
        all_evidence_id_all.reset_index(drop=True,inplace=True)
        evidence_link_data=[]
        for evidence_level in ['Level A','Level B','Level C','Level D']:
            all_evidence_id_all_eg=all_evidence_id_all[all_evidence_id_all['Evidence_level']==evidence_level]
            if all_evidence_id_all_eg.shape[0]!=0:
                evidence_eg=evidence_level
                count_eg=int(all_evidence_id_all_eg['count'])
                evidence_link_data_eg={
                    'name':evidence_eg,
                    'value':count_eg
                }
            else:
                evidence_link_data_eg={
                    'name':evidence_level,
                    'value':0
                }
            evidence_link_data=evidence_link_data+[evidence_link_data_eg]
        evidence_option_json={
        'title': {
        'text': 'Statistical evidence levels for data entries','left': 'center',
        'top': 25,    
        'left': 'center'
        },
            'tooltip': {
            'trigger': 'item'
            },
            'legend': {
            'bottom': '12.5%',
        'left': 'center'
            },
            'series': [
            {
                'name': 'Access From',
                'type': 'pie',
                'radius': ['30%', '60%'],
                'avoidLabelOverlap': 'false',
                'itemStyle': {
                'borderRadius': 10,
                'borderColor': '#fff',
                'borderWidth': 2
                },
                'label': {
                'show': 'false',
                'position': 'center'
                },
                'emphasis': {
                'label': {
                    'show': 'true',
                    'fontSize': 20,
                    'fontWeight': 'bold'
                }
                },
                'labelLine': {
                'show': 'false'
                },
                'data': evidence_link_data
            }
            ]
        };

        evidence_option_json = json.dumps(evidence_option_json)
        
        all_cancer_type_id_all=pd.DataFrame(action_mutation_data['Cancer_type_oncotree_Level2'].value_counts())
        all_cancer_type_id_all=all_cancer_type_id_all.iloc[0:10,:]
        cancer_type_link_data=[]
        for pos in range(all_cancer_type_id_all.shape[0]):
            all_cancer_type_id_all_eg=all_cancer_type_id_all.iloc[pos:(pos+1),:]
            cancer_type_eg=str(all_cancer_type_id_all_eg.index.tolist()[0])
            count_eg=int(all_cancer_type_id_all_eg['count'])
            cancer_type_link_data_eg={
                'name':cancer_type_eg,
                'value':count_eg
            }
            cancer_type_link_data=cancer_type_link_data+[cancer_type_link_data_eg]

            cancer_type_option_json={
                'title': {
                        'text': 'Cancer types associated with biomarkers','left': 'center','top':25
                        }
                ,'tooltip': {},
        'series': [ {
        'type': 'wordCloud',
        'gridSize': 2,
        'sizeRange': [10, 30],
        'rotationRange': [-90, 90],
        'shape': 'apple',
        'width': 600,
        'height': 400,
        'drawOutOfBound': 'true',
        'textStyle': {
        'fontFamily': 'sans-serif',
        'fontWeight': 'bold',
            'emphasis': {
            'shadowBlur': 10,
            'shadowColor': '#333'
        }
        },'data':cancer_type_link_data}]
                                };

        cancer_type_option_json = json.dumps(cancer_type_option_json)
        
        return render_template('each_biomarker_cancer_search_info.html',actmut=actmut,action_mutation_data=action_mutation_data.to_dict(orient='records'),option_json=option_json,drug_option_json=drug_option_json,evidence_option_json=evidence_option_json,cancer_type_option_json=cancer_type_option_json)
    except:
        return render_template('navbar_error_biomarker.html')
        
def vcf_match_drug(pos,original_database,transvar_all):
    original_database_eg=original_database.iloc[pos:(pos+1),:].fillna('No info')
    all_action_mut_gene_eg=original_database_eg['Action mutation genes'].str.split('; ',expand=True)
    all_action_mut_event_eg=original_database_eg['Action mutation events'].str.split('; ',expand=True)
    all_action_mut_gene_event=list(all_action_mut_gene_eg.iloc[0,:]+':p.'+all_action_mut_event_eg.iloc[0,:])
    if len(np.intersect1d(all_action_mut_gene_event,transvar_all['hgvsp'].tolist()))==len(all_action_mut_gene_event):
        return(original_database_eg)
    else:
        pass
@blue.route('/oncodrug/saveinput',methods=['POST'])
def saveinput():
    try:
        db.create_all()
    except:
        pass
    db.session.query(Match_VCF).delete()
    db.session.commit()
    try:
        if request.method=='POST':
            current_time = datetime.now()
            formatted_time = current_time.strftime("%Y%m%d_%H%M%S")
            filename = f'{formatted_time}_users.csv'
            input_text = request.form['inputText']
            uploaded_file =request.files['file']
            if input_text:
                save_to_file(input_text,filename)
            else:
                save_to_file(uploaded_file.read().decode('utf-8'), filename)
            ###上传文件后读取文件，然后进行后台数据分析###
            vcf_read=pd.read_csv(f'/h/tianyi/oncodrug/upload/{filename}',sep='\t',header=None)
            vcf_read[0]=vcf_read[0].str.replace(' ','')
            ###
            transvar_all=pd.DataFrame()
            for pos in range(vcf_read.shape[0]):
                vcf_read_eg=''.join(vcf_read.iloc[pos,:])
                if ':p.' in vcf_read_eg:
                    transvar_eg=pd.DataFrame(vcf_read_eg.split(':p.')).T
                    transvar_eg.columns=['Action_mut_gene','Action_mut_event']
                elif ':g.' in vcf_read_eg:
                    transvar_eg=transvar_translation(vcf_read_eg,'hgvsg')
                elif ':c.' in vcf_read_eg:
                    transvar_eg=transvar_translation(vcf_read_eg,'hgvsc')
                elif 'DEL' in vcf_read_eg:
                    transvar_eg=pd.DataFrame(vcf_read_eg.split(':')).T
                    transvar_eg.columns=['Action_mut_gene','Action_mut_event']
                elif 'AMP' in vcf_read_eg:
                    transvar_eg=pd.DataFrame(vcf_read_eg.split(':')).T
                    transvar_eg.columns=['Action_mut_gene','Action_mut_event']
                elif '__' in vcf_read_eg:
                    transvar_eg=pd.DataFrame()
                    transvar_eg['Action_mut_gene']=[vcf_read_eg]
                    transvar_eg['Action_mut_event']=[vcf_read_eg]
                transvar_all=pd.concat([transvar_eg,transvar_all])
            transvar_all=transvar_all.drop_duplicates()
            transvar_all['hgvsp']=transvar_all['Action_mut_gene']+':p.'+transvar_all['Action_mut_event']
            vcf_match_drug_funtools=functools.partial(vcf_match_drug,original_database=original_database,transvar_all=transvar_all)
            try:
                pool = Pool(5)
                drug_match =pd.concat(pool.map(vcf_match_drug_funtools,range(original_database.shape[0])) )
                pool.close()
                pool.join()
            except:
                drug_match=pd.DataFrame(columns=original_database.columns.tolist())
            drug_match_comdb=drug_match
            drug_match_comdb.to_csv(f'/h/tianyi/oncodrug/download/{filename}_download',sep='\t',index=False)
            
            for pos in range(drug_match_comdb.shape[0]):
                drugcomb_data_eg=drug_match_comdb.iloc[pos:(pos+1),]
                Drug_combination_ID= ''.join(drugcomb_data_eg['Drug combination ID'].astype(str))
                Evidence_level= ''.join(drugcomb_data_eg['Evidence level'].astype(str))
                Prioritization_score= int(drugcomb_data_eg['Prioritization score'].astype(int))
                Targeted_drug= ''.join(drugcomb_data_eg['Targeted drug'].astype(str))
                Non_targeted_drug= ''.join(drugcomb_data_eg['Non-targeted drug'].astype(str))
                Cancer_type_oncotree_Level2= ''.join(drugcomb_data_eg['Cancer type (oncotree Level2)'].astype(str))
                Action_mutations= ''.join(drugcomb_data_eg['Action mutations'].astype(str))
                Response= ''.join(drugcomb_data_eg['Response'].astype(str))
                Adverse_effect= ''.join(drugcomb_data_eg['Adverse effect'].astype(str))
                Drug_dosage= ''.join(drugcomb_data_eg['Drug dosage'].astype(str))
                Evidence_level_score= int(drugcomb_data_eg['Evidence level score'].astype(int))
                Actionable_mutation_precision_score= int(drugcomb_data_eg['Actionable mutation precision score'].astype(int))
                Action_mutation_level_score= int(drugcomb_data_eg['Action mutation level score'].astype(int))
                FDA_evidence_score= int(drugcomb_data_eg['FDA evidence score'].astype(int))
                Response_score=int(drugcomb_data_eg['Response score'].astype(int))


                drug_match_comdb_content=Match_VCF(Drug_combination_ID=Drug_combination_ID,
                    Evidence_level=Evidence_level,
                    Prioritization_score=Prioritization_score,
                    Targeted_drug=Targeted_drug,
                    Non_targeted_drug=Non_targeted_drug,
                    Cancer_type_oncotree_Level2=Cancer_type_oncotree_Level2,
                    Action_mutations=Action_mutations,
                    Response=Response,
                    Adverse_effect=Adverse_effect,
                    Drug_dosage=Drug_dosage,
                    Evidence_level_score=Evidence_level_score,
                    Action_mutation_level_score=Action_mutation_level_score,
                    FDA_evidence_score=FDA_evidence_score,
                    Actionable_mutation_precision_score=Actionable_mutation_precision_score,
                    Response_score=Response_score)
                db.session.add(drug_match_comdb_content)
            ####处理后的数据传入数据库并进行网页渲染
            db.session.commit()
            drug_match_comdb_content_all=Match_VCF.query.all()
            return render_template('biomarker.html',drug_match_comdb_content_all=drug_match_comdb_content_all )
        else:
            return render_template('biomarker_error.html')
    except:
        return render_template('biomarker_error.html' )


@blue.route('/oncodrug/download_file')
def download_file():
    download_databases_file = Match_VCF.query.all()
    download_databases_file_dic = []
    for download_databases_file_eg in download_databases_file:
        download_databases_file_dic_eg = {
            'Drug_combination_ID': download_databases_file_eg.Drug_combination_ID,
            'Evidence_level': download_databases_file_eg.Evidence_level,
            'Prioritization_score': download_databases_file_eg.Prioritization_score,
            'Targeted_drug': download_databases_file_eg.Targeted_drug,
            'Non_targeted_drug': download_databases_file_eg.Non_targeted_drug,
            'Cancer_type_oncotree_Level2': download_databases_file_eg.Cancer_type_oncotree_Level2,
            'Biomarkers': download_databases_file_eg.Action_mutations,
            'Response': download_databases_file_eg.Response,
            'Adverse_effect': download_databases_file_eg.Adverse_effect,
            'Drug_dosage': download_databases_file_eg.Drug_dosage
        }
        download_databases_file_dic.append(download_databases_file_dic_eg)
    filename=create_csv(download_databases_file_dic)
    return send_file(create_csv(download_databases_file_dic), as_attachment=True, download_name=f'{filename}')


#####数据库删除
@blue.route('/oncodrug/delete_matchvcf')
def delete_table():
    table = Table('match_vcf', db.metadata, autoload_with=db.engine)
    table.drop(db.engine)
    return '表已删除'

@blue.route('/oncodrug/delete_drugcom')
def delete_table_drugcom():
    table = Table('drug_comdb', db.metadata, autoload_with=db.engine)
    table.drop(db.engine)
    return '表已删除'


@blue.route('/oncodrug/delete_drugid/<int:Drug_combination_ID>', methods=['GET', 'POST'])
def delete_user(Drug_combination_ID):
    user_id = DrugComdb.query.get(Drug_combination_ID)
    if user_id:
        db.session.delete(user_id)
        db.session.commit()
        return '用户已删除'
    else:
        return '用户不存在'

@blue.route('/oncodrug/createdb')
def create_db():
    ####create the table schema in the database
    db.create_all()
    return 'DB success'

####数据库插入
@blue.route('/oncodrug/adddrugcom')
def adddrugcom():
    ###对数据进行操作并存储
    drugcomb_data=pd.read_excel('/h/tianyi/oncodrug/upload/drug_comdb_data/drug_combination_database.xlsx')
    for pos in range(drugcomb_data.shape[0]):
        drugcomb_data_eg=drugcomb_data.iloc[pos:(pos+1),]
        Drug_combination_ID= ''.join(drugcomb_data_eg['Drug combination ID'].astype(str))
        Evidence_level= ''.join(drugcomb_data_eg['Evidence level'].astype(str))
        Prioritization_score= int(drugcomb_data_eg['Prioritization score'].astype(int))
        Targeted_drug= ''.join(drugcomb_data_eg['Targeted drug'].astype(str))
        Non_targeted_drug= ''.join(drugcomb_data_eg['Non-targeted drug'].astype(str))
        Cancer_type_oncotree_Level2= ''.join(drugcomb_data_eg['Cancer type (oncotree Level2)'].astype(str))
        Action_mutations= ''.join(drugcomb_data_eg['Action mutations'].astype(str))
        Response= ''.join(drugcomb_data_eg['Response'].astype(str))
        Adverse_effect= ''.join(drugcomb_data_eg['Adverse effect'].astype(str))
        Drug_dosage= ''.join(drugcomb_data_eg['Drug dosage'].astype(str))
        Evidence_level_score= int(drugcomb_data_eg['Evidence level score'].astype(int))
        Actionable_mutation_precision_score= int(drugcomb_data_eg['Actionable mutation precision score'].astype(int))
        Action_mutation_level_score= int(drugcomb_data_eg['Action mutation level score'].astype(int))
        FDA_evidence_score= int(drugcomb_data_eg['FDA evidence score'].astype(int))
        Response_score=int(drugcomb_data_eg['Response score'].astype(int))

        
        DrugComdb_content=DrugComdb(Drug_combination_ID=Drug_combination_ID,
        Evidence_level=Evidence_level,
        Prioritization_score=Prioritization_score,
        Targeted_drug=Targeted_drug,
        Non_targeted_drug=Non_targeted_drug,
        Cancer_type_oncotree_Level2=Cancer_type_oncotree_Level2,
        Action_mutations=Action_mutations,
        Response=Response,
        Adverse_effect=Adverse_effect,
        Drug_dosage=Drug_dosage,
        Evidence_level_score=Evidence_level_score,
        Action_mutation_level_score=Action_mutation_level_score,
        FDA_evidence_score=FDA_evidence_score,
        Actionable_mutation_precision_score=Actionable_mutation_precision_score,
        Response_score=Response_score)
        db.session.add(DrugComdb_content)
    db.session.commit()

    return 'drugcomb Success'


@blue.route('/oncodrug/browse_levela')
def browse_levela():
    ####对应models中person()
    DrugComdb_content_levela=DrugComdb.query.filter(DrugComdb.Evidence_level == 'Level A').all()
    return render_template('browse_levela.html',DrugComdb_content=DrugComdb_content_levela)

@blue.route('/oncodrug/browse_levela_test')
def browse_levela_test():
    ####对应models中person()
    DrugComdb_content_levela=DrugComdb.query.filter(DrugComdb.Evidence_level == 'Level A').all()
    return render_template('browse_levela_test.html',DrugComdb_content=DrugComdb_content_levela)

@blue.route('/oncodrug/browse_levelb')
def browse_levelb():
    ####对应models中person()
    DrugComdb_content_levelb=DrugComdb.query.filter(DrugComdb.Evidence_level == 'Level B').all()
    return render_template('browse_levelb.html',DrugComdb_content=DrugComdb_content_levelb)

@blue.route('/oncodrug/browse_levelc')
def browse_levelc():
    ####对应models中person()
    DrugComdb_content_levelc=DrugComdb.query.filter(DrugComdb.Evidence_level == 'Level C').all()
    return render_template('browse_levelc.html',DrugComdb_content=DrugComdb_content_levelc)


@blue.route('/oncodrug/browse_leveld')
def browse_leveld():
    ####对应models中person()
    DrugComdb_content_leveld=DrugComdb.query.filter(DrugComdb.Evidence_level == 'Level D').all()
    return render_template('browse_leveld.html',DrugComdb_content=DrugComdb_content_leveld)

####将数据库id参数传入函数中，然后获取超链接网页
@blue.route('/oncodrug/item/<int:item_id>')
def item_details(item_id):

    # 根据item_id查询数据库中的数据详情
    item = DrugComdb.query.get(item_id)
    drug_name_target_all=item.Targeted_drug
    drug_name_other_all=item.Non_targeted_drug
    drug_name_all=drug_name_target_all.split(';')+drug_name_other_all.split(';')
    drug_name_all=pd.DataFrame(drug_name_all).dropna()
    drug_name_all.columns=['drug_name']
    drug_name_all['names']='name'
    drug_name_all=drug_name_all[drug_name_all['drug_name']!='nan']
    drug_name_all=drug_name_all[drug_name_all['drug_name']!='No info']
#####drug name and drugbank id####
    drugbank_target_all=item.Drugbank_id_of_targeted_drug
    drugbank_other_all=item.Drugbank_id_of_non_targeted_drug
    drugbank_all=drugbank_target_all.split(';')+drugbank_other_all.split(';')
    drugbank_all=pd.DataFrame(drugbank_all).dropna()
    drugbank_all.columns=['drugbank_id']
    drugbank_all=drugbank_all[drugbank_all['drugbank_id']!='nan']

    link_combi=pd.DataFrame(list(combinations(drug_name_all['drug_name'].tolist(),2)))
    link_combi.columns=['node1','node2']

    #####nodes position####


    x_all,y_all=x_y_drug_pos(drug_name_all)

    ###acition mutation
    mut_gene_all=item.Action_mutation_genes
    mut_event_all=item.Action_mutation_events

    mut_gene_str_all=pd.Series(mut_gene_all.split('; '))
    mut_event_str_all=pd.Series(mut_event_all.split('; '))
    action_mut_all=';'.join(mut_gene_str_all+' '+mut_event_str_all)

    node_data=[]
    for pos in range(drug_name_all.shape[0]):
        drug_name_eg=''.join(drug_name_all.iloc[pos:(pos+1),:]['drug_name'])
        drugbank_eg=''.join(drugbank_all.iloc[pos:(pos+1),:]['drugbank_id'])
        x_eg=x_all[pos]
        y_eg=y_all[pos]
        node_data_eg={
            'name': f'{drug_name_eg}',
            'x': x_eg,
            'y': y_eg,
            'tooltip': f'{drug_name_eg}: {drugbank_eg}'
            }
        node_data=node_data+[node_data_eg]

    link_data=[]
    for link_pos in range(link_combi.shape[0]):
        source_eg=''.join(link_combi.iloc[link_pos:(link_pos+1),:]['node1'])
        target_eg=''.join(link_combi.iloc[link_pos:(link_pos+1),:]['node2'])
        link_data_eg={
        'source': f'{source_eg}',
        'target': f'{target_eg}',
        'label': {
        'show': 'true','formatter': ''
        },
        'lineStyle': {
        'curveness': 0,
        },'tooltip':f'{action_mut_all}'
        }
        link_data=link_data+[link_data_eg]
    
    option = {'title': {
          'text': 'Drug combination strategies in entry','left': 'center'
          },
                  'tooltip': {},
                  'animationDurationUpdate': 1500,
                  'animationEasingUpdate': 'quinticInOut',
                  'series': [
                  {
                    'type': 'graph',
                    'layout': 'none',
                    'symbolSize': [25/2, 25/2],
                    'roam': 'true',
                    'label': {
                    'show': 'true','formatter': '{b}'
                    },
                    'edgeLabel': {
                    'fontSize': 20
                    },
                    'itemStyle': {
                    'borderColor': '#000',
                    'borderWidth': 1,
                    'borderType': 'solid',
                    'shadowBlur': 10,
                    'shadowColor': 'rgba(0, 0, 0, 0.3)'},
                    
                    'data': node_data,
                    'links': link_data,
                    'lineStyle': {
                    'opacity': 0.9,
                    'width': 2,'smooth': 'false'
                    }
                  }
                  ]
                };
    # 将 option 对象传递给模板渲染
    option_json = json.dumps(option)

    return render_template('each_drugcomb_id_info.html',item=item,option_json=option_json)



@blue.route('/oncodrug/each_drug_info/<drugname>')
def each_drug_info(drugname):
    try:
        ####对应models中person()
            ####对应models中person()
        all_combination_pairs_part1=original_database[original_database['Targeted drug'].str.contains(drugname).fillna(False)]
        all_combination_pairs_part2=original_database[original_database['Non-targeted drug'].str.contains(drugname).fillna(False)]
        all_combination_pairs=pd.concat([all_combination_pairs_part1,all_combination_pairs_part2])
        all_combination_pairs.columns=pd.Series(all_combination_pairs.columns.tolist()).str.replace(' ','_').str.replace('-','_').str.replace('(','').str.replace(')','')

        drug_name_all=np.union1d(all_combination_pairs['Targeted_drug'].astype(str).tolist(),all_combination_pairs['Non_targeted_drug'].astype(str).tolist())
        drug_name_all=pd.Series(drug_name_all).dropna()

        drug_name_all_str=[]
        for drug_name_eg in drug_name_all:
            drug_name_eg_str=drug_name_eg.split('; ')
            drug_name_all_str=drug_name_eg_str+drug_name_all_str

        drug_name_all=pd.DataFrame(drug_name_all_str).drop_duplicates()
        drug_name_all.columns=['drug_name']
        drug_name_all['names']='name'
        drug_name_all=drug_name_all[(drug_name_all['drug_name']!='nan') & (drug_name_all['drug_name']!='No info')]



        ##drugbank id###
        all_drugbank_info_all=pd.merge(drug_name_all,drugbank_id_raw_data,left_on=drug_name_all['drug_name'].str.upper(),right_on=drugbank_id_raw_data['All_Drug'].str.upper(),how='left')[['drug_name','DrugBank']]
        all_drugbank_info_all.columns=['drugbank_name','drugbank_id']
        all_drugbank_info_all=all_drugbank_info_all.dropna().drop_duplicates(['drugbank_id'])

        link_combi=all_drugbank_info_all[['drugbank_name']]
        link_combi['names']=drugname
        link_combi.columns=['node1','node2']
        link_combi=link_combi[link_combi['node1'] != link_combi['node2']]

    ####核心节点
        drugbank_core_data=all_drugbank_info_all[all_drugbank_info_all['drugbank_name']==drugname]
        if drugbank_core_data.shape[0]!=0:
            core_drug_name_eg=''.join(drugbank_core_data['drugbank_name'])
            core_drugbank_eg=''.join(drugbank_core_data['drugbank_id'])
        else:
            core_drug_name_eg=drugname
            core_drugbank_eg=''
        core_node_data={
                'name': f'{core_drug_name_eg}',
                'x': 0,
                'y': 0,'symbolSize': [25, 25],'itemStyle': { 'color': '#67C1C4' } ,
                'tooltip': f'{core_drug_name_eg}: {core_drugbank_eg}'
                }
    ###其他节点
        #####nodes position####
        if len(all_drugbank_info_all['drugbank_name'].drop_duplicates())!=1:
            all_comb_drug_name_drop=np.setdiff1d(all_drugbank_info_all['drugbank_name'],drugname)
        else:
            all_comb_drug_name_drop=drugname
        x_all,y_all=single_drug_pos(len(all_comb_drug_name_drop))
        node_data_nocore=[]
        for pos in range(len(all_comb_drug_name_drop)):
            drugname_eg=all_comb_drug_name_drop[pos]
            all_drugbank_info_eg=all_drugbank_info_all[all_drugbank_info_all['drugbank_name']==drugname_eg]
            if all_drugbank_info_eg.shape[0]==0:
                drug_name_eg=drugname_eg
                drugbank_eg='No info'
            else:
                drug_name_eg=''.join(all_drugbank_info_eg['drugbank_name'])
                drugbank_eg=''.join(all_drugbank_info_eg['drugbank_id'])
            x_eg=x_all[pos]
            y_eg=y_all[pos]
            node_data_eg={
                'name': f'{drug_name_eg}',
                'x': x_eg,
                'y': y_eg,
                'tooltip': f'{drug_name_eg}: {drugbank_eg}'
                }
            node_data_nocore=node_data_nocore+[node_data_eg]

        node_data=[core_node_data]+node_data_nocore
        link_data=[]
        for link_pos in range(link_combi.shape[0]):
            source_eg=''.join(link_combi.iloc[link_pos:(link_pos+1),:]['node1'])
            target_eg=''.join(link_combi.iloc[link_pos:(link_pos+1),:]['node2'])
            link_data_eg={
            'source': f'{source_eg}',
            'target': f'{target_eg}',
            'label': {
            'show': 'true','formatter': ''
            },
            'lineStyle': {
            'curveness': 0,
            }
            }
            link_data=link_data+[link_data_eg]
        
        option = {
                'title': {
                'text': 'Drug combination therapy regimens','left': 'center','top': 25, 
                },
                    'tooltip': {},
                    'animationDurationUpdate': 1500,
                    'animationEasingUpdate': 'quinticInOut',
                    'series': [
                    {
                        'type': 'graph',
                        'layout': 'none',
                        'symbolSize': [25/2, 25/2],
                        'roam': 'true',
                        'label': {
                        'show': 'true','formatter': '{b}'
                        },
                        'edgeLabel': {
                        'fontSize': 20
                        },
                        'itemStyle': {
                        'borderColor': '#000',
                        'borderWidth': 1,
                        'borderType': 'solid',
                        'shadowBlur': 10,
                        'shadowColor': 'rgba(0, 0, 0, 0.3)'},
                        
                        'data': node_data,
                        'links': link_data,
                        'lineStyle': {
                        'opacity': 0.9,
                        'width': 2,'smooth': 'false'
                        }
                    }
                    ]
                    }
        # 将 option 对象传递给模板渲染
        option_json = json.dumps(option)
        if core_drugbank_eg!='':
            drugbank_link=f'https://go.drugbank.com/drugs/{core_drugbank_eg}'
        else:
            drugbank_link=''
        #####all drug info detail#####
        drug_info_detail_eg=drug_info_detail_all[drug_info_detail_all['drugbank_id']==core_drugbank_eg]
        drug_info_detail_description1=drug_info_detail_eg['description'].dropna()
        if len(drug_info_detail_description1)==0:
            drug_info_detail_description2=''
        else:
            drug_info_detail_description2=''.join(drug_info_detail_eg['description'])

        drug_info_detail_target1=drug_info_detail_eg['target'].dropna()
        if len(drug_info_detail_target1)==0:
            drug_info_detail_target2=''
        else:
            drug_info_detail_target2=''.join(drug_info_detail_eg['target'])
        
        drug_info_detail_moa1=drug_info_detail_eg['moa'].dropna()
        if len(drug_info_detail_moa1)==0:
            drug_info_detail_moa2=''
        else:
            drug_info_detail_moa2=''.join(drug_info_detail_eg['moa'])

        drug_info_detail_pharmacodynamics1=drug_info_detail_eg['pharmacodynamics'].dropna()
        if len(drug_info_detail_pharmacodynamics1)==0:
            drug_info_detail_pharmacodynamics2=''
        else:
            drug_info_detail_pharmacodynamics2=''.join(drug_info_detail_eg['pharmacodynamics'])

        drug_info_detail_indication1=drug_info_detail_eg['indication'].dropna()
        if len(drug_info_detail_indication1)==0:
            drug_info_detail_indication2=''
        else:
            drug_info_detail_indication2=''.join(drug_info_detail_eg['indication'])
            
        all_cancer_type_id_all=pd.DataFrame(all_combination_pairs['Cancer_type_oncotree_Level2'].value_counts()).iloc[0:10,:]
        cancer_type_link_data=[]
        for pos in range(all_cancer_type_id_all.shape[0]):
            all_cancer_type_id_all_eg=all_cancer_type_id_all.iloc[pos:(pos+1),:]
            cancer_type_eg=str(all_cancer_type_id_all_eg.index.tolist()[0])
            count_eg=int(all_cancer_type_id_all_eg['count'])
            cancer_type_link_data_eg={
                'name':cancer_type_eg,
                'value':count_eg
            }
            cancer_type_link_data=cancer_type_link_data+[cancer_type_link_data_eg]

            cancer_type_option_json={
                'title': {
                        'text': 'Cancer types','left': 'center','top': 25, 
                        }
                ,'tooltip': {},
        'series': [ {
        'type': 'wordCloud',
        'gridSize': 2,
        'sizeRange': [10, 30],
        'rotationRange': [-90, 90],
        'shape': 'apple',
        'width': 600,
        'height': 400,
        'drawOutOfBound': 'true',
        'textStyle': {
        'fontFamily': 'sans-serif',
        'fontWeight': 'bold',
            'emphasis': {
            'shadowBlur': 10,
            'shadowColor': '#333'
        }
        },'data':cancer_type_link_data}]
                                };

        cancer_type_option_json = json.dumps(cancer_type_option_json)
            
        all_evidence_id_all=pd.DataFrame(all_combination_pairs['Evidence_level'].value_counts())
        all_evidence_id_all.columns=['count']
        all_evidence_id_all['Evidence_level']=all_evidence_id_all.index.tolist()
        all_evidence_id_all.reset_index(drop=True,inplace=True)
        evidence_link_data=[]
        for evidence_level in ['Level A','Level B','Level C','Level D']:
            all_evidence_id_all_eg=all_evidence_id_all[all_evidence_id_all['Evidence_level']==evidence_level]
            if all_evidence_id_all_eg.shape[0]!=0:
                evidence_eg=evidence_level
                count_eg=int(all_evidence_id_all_eg['count'])
                evidence_link_data_eg={
                    'name':evidence_eg,
                    'value':count_eg
                }
            else:
                evidence_link_data_eg={
                    'name':evidence_level,
                    'value':0
                }
            evidence_link_data=evidence_link_data+[evidence_link_data_eg]
        evidence_option_json={
        'title': {
        'text': 'Statistical evidence levels for data entries','left': 'center',
        'top': 25,    
        'left': 'center'
        },
            'tooltip': {
            'trigger': 'item'
            },
            'legend': {
            'bottom': '12.5%',
        'left': 'center'
            },
            'series': [
            {
                'name': 'Access From',
                'type': 'pie',
                'radius': ['30%', '60%'],
                'avoidLabelOverlap': 'false',
                'itemStyle': {
                'borderRadius': 10,
                'borderColor': '#fff',
                'borderWidth': 2
                },
                'label': {
                'show': 'false',
                'position': 'center'
                },
                'emphasis': {
                'label': {
                    'show': 'true',
                    'fontSize': 20,
                    'fontWeight': 'bold'
                }
                },
                'labelLine': {
                'show': 'false'
                },
                'data': evidence_link_data
            }
            ]
        };

        evidence_option_json = json.dumps(evidence_option_json)
        
        all_mutation_id_all=pd.DataFrame(all_combination_pairs['Action_mutations'].value_counts())
        all_mutation_id_all=all_mutation_id_all.iloc[0:10,:]
        all_mutation_id_all=all_mutation_id_all.loc[list(~(pd.Series(all_mutation_id_all.index.tolist()).str.contains('Unavailable|No info',case=False))),:]
        mutation_link_data=[]
        for pos in range(all_mutation_id_all.shape[0]):
            all_mutation_id_all_eg=all_mutation_id_all.iloc[pos:(pos+1),:]
            mutation_eg=str(all_mutation_id_all_eg.index.tolist()[0])
            count_eg=int(all_mutation_id_all_eg['count'])
            mutation_link_data_eg={
                'name':mutation_eg,
                'value':count_eg
            }
            mutation_link_data=mutation_link_data+[mutation_link_data_eg]

            mutation_option_json={
                'title': {
                        'text': 'Biomarkers predicting drug responses','left': 'center','top': 25,
                        }
                ,'tooltip': {},
        'series': [ {
        'type': 'wordCloud',
        'gridSize': 2,
        'sizeRange': [10, 30],
        'rotationRange': [-90, 90],
        'shape': 'apple',
        'width': 600,
        'height': 400,
        'drawOutOfBound': 'true',
        'textStyle': {
        'fontFamily': 'sans-serif',
        'fontWeight': 'bold',
            'emphasis': {
            'shadowBlur': 10,
            'shadowColor': '#333'
        }
        },'data':mutation_link_data}]
                                };

        mutation_option_json = json.dumps(mutation_option_json)
        return render_template('each_drug_info.html',core_drugbank_eg=core_drugbank_eg,drugbank_link=drugbank_link,drugname=drugname,drug_info_detail_description=drug_info_detail_description2,drug_info_detail_target=drug_info_detail_target2,drug_info_detail_moa=drug_info_detail_moa2,drug_info_detail_pharmacodynamics=drug_info_detail_pharmacodynamics2,drug_info_detail_indication=drug_info_detail_indication2,all_combination_pairs=all_combination_pairs.to_dict(orient='records'),option_json=option_json,cancer_type_option_json=cancer_type_option_json,evidence_option_json=evidence_option_json,mutation_option_json=mutation_option_json)
    except:
        return render_template('navbar_error.html')


#####show drug combination information#######

@blue.route('/oncodrug/each_combdrug_info/<drugname>')
def each_combdrug_info(drugname):
    ####对应models中person()
        ####对应models中person()
    drugname1=drugname
    drugname2='; '.join(sorted(drugname1.split('; ')))
    if drugname1==drugname1:
        drugname2='; '.join(sorted(drugname1.split('; '))[::-1])
    ####对应models中person()
        ####对应models中person()
    all_combination_pairs_part1=original_database[original_database['Targeted drug'].str.contains(f'{drugname1}|{drugname2}',case=False).fillna(False)]
    all_combination_pairs_part2=original_database[original_database['Non-targeted drug'].str.contains(f'{drugname1}|{drugname2}',case=False).fillna(False)]
    all_combination_pairs=pd.concat([all_combination_pairs_part1,all_combination_pairs_part2])
    all_combination_pairs.columns=pd.Series(all_combination_pairs.columns.tolist()).str.replace(' ','_').str.replace('-','_').str.replace('(','').str.replace(')','')
    drug_name_all=np.union1d(all_combination_pairs['Targeted_drug'].astype(str).tolist(),all_combination_pairs['Non_targeted_drug'].astype(str).tolist())
    drug_name_all=pd.Series(drug_name_all).dropna()
    drug_name_all_str=[]
    for drug_name_eg in drug_name_all:
        drug_name_eg_str=drug_name_eg.split('; ')
        drug_name_all_str=drug_name_eg_str+drug_name_all_str

    drug_name_all=pd.DataFrame(drug_name_all_str).drop_duplicates()
    drug_name_all.columns=['drug_name']
    drug_name_all['names']='name'
    drug_name_all=drug_name_all[(drug_name_all['drug_name']!='nan') & (drug_name_all['drug_name']!='No info')]

        
    all_cancer_type_id_all=pd.DataFrame(all_combination_pairs['Cancer_type_oncotree_Level2'].value_counts()).iloc[0:10,:]
    all_cancer_type_id_all.columns=['count']
    cancer_type_link_data=[]
    for pos in range(all_cancer_type_id_all.shape[0]):
        all_cancer_type_id_all_eg=all_cancer_type_id_all.iloc[pos:(pos+1),:]
        cancer_type_eg=str(all_cancer_type_id_all_eg.index.tolist()[0])
        count_eg=int(all_cancer_type_id_all_eg['count'])
        cancer_type_link_data_eg={
            'name':cancer_type_eg,
            'value':count_eg
        }
        cancer_type_link_data=cancer_type_link_data+[cancer_type_link_data_eg]
        cancer_type_option_json={
            'title': {
                        'text': 'Cancer types','left': 'center','top': 25, 
                        }
            ,'tooltip': {},
    'series': [ {
    'type': 'wordCloud',
    'gridSize': 2,
    'sizeRange': [10, 30],
    'rotationRange': [-90, 90],
    'shape': 'apple',
    'width': 600,
    'height': 400,
    'drawOutOfBound': 'true',
    'textStyle': {
    'fontFamily': 'sans-serif',
    'fontWeight': 'bold',
        'emphasis': {
        'shadowBlur': 10,
        'shadowColor': '#333'
    }
    },'data':cancer_type_link_data}]
                            };

    cancer_type_option_json = json.dumps(cancer_type_option_json)
        
    all_evidence_id_all=pd.DataFrame(all_combination_pairs['Evidence_level'].value_counts())
    all_evidence_id_all.columns=['count']
    all_evidence_id_all['Evidence_level']=all_evidence_id_all.index.tolist()
    all_evidence_id_all.reset_index(drop=True,inplace=True)
    evidence_link_data=[]
    for evidence_level in ['Level A','Level B','Level C','Level D']:
        all_evidence_id_all_eg=all_evidence_id_all[all_evidence_id_all['Evidence_level']==evidence_level]
        if all_evidence_id_all_eg.shape[0]!=0:
            evidence_eg=evidence_level
            count_eg=int(all_evidence_id_all_eg['count'])
            evidence_link_data_eg={
                'name':evidence_eg,
                'value':count_eg
            }
        else:
            evidence_link_data_eg={
                'name':evidence_level,
                'value':0
            }
        evidence_link_data=evidence_link_data+[evidence_link_data_eg]
    evidence_option_json={
        'title': {
        'text': 'Statistical evidence levels for data entries','left': 'center',
        'top': 25,    
        'left': 'center'
        },
            'tooltip': {
            'trigger': 'item'
        },
        'legend': {
            'bottom': '12.5%',
        'left': 'center'
        },
        'series': [
            {
            'name': 'Access From',
            'type': 'pie',
            'radius': ['30%', '60%'],
            'avoidLabelOverlap': 'false',
            'itemStyle': {
                'borderRadius': 10,
                'borderColor': '#fff',
                'borderWidth': 2
            },
            'label': {
                'show': 'false',
                'position': 'center'
            },
            'emphasis': {
                'label': {
                'show': 'true',
                'fontSize': 20,
                'fontWeight': 'bold'
                }
            },
            'labelLine': {
                'show': 'false'
            },
            'data': evidence_link_data
            }
        ]
        };

    evidence_option_json = json.dumps(evidence_option_json)

    all_mutation_id_all=pd.DataFrame(all_combination_pairs['Action_mutations'].value_counts())
    all_mutation_id_all.columns=['count']
    all_mutation_id_all=all_mutation_id_all.iloc[0:10,:]
    all_mutation_id_all=all_mutation_id_all.loc[list(~(pd.Series(all_mutation_id_all.index.tolist()).str.contains('Unavailable|No info'))),:]
    mutation_link_data=[]
    for pos in range(all_mutation_id_all.shape[0]):
        all_mutation_id_all_eg=all_mutation_id_all.iloc[pos:(pos+1),:]
        mutation_eg=str(all_mutation_id_all_eg.index.tolist()[0])
        count_eg=int(all_mutation_id_all_eg['count'])
        mutation_link_data_eg={
            'name':mutation_eg,
            'value':count_eg
        }
        mutation_link_data=mutation_link_data+[mutation_link_data_eg]

        mutation_option_json={
            'title': {
                        'text': 'Biomarkers predicting drug responses','left': 'center','top': 25,
                        }
            ,'tooltip': {},
    'series': [ {
    'type': 'wordCloud',
    'gridSize': 2,
    'sizeRange': [10, 30],
    'rotationRange': [-90, 90],
    'shape': 'apple',
    'width': 600,
    'height': 400,
    'drawOutOfBound': 'true',
    'textStyle': {
    'fontFamily': 'sans-serif',
    'fontWeight': 'bold',
        'emphasis': {
        'shadowBlur': 10,
        'shadowColor': '#333'
    }
    },'data':mutation_link_data}]
                            };

    mutation_option_json = json.dumps(mutation_option_json)
    return render_template('each_combdrug_info.html',drugname=drugname,all_combination_pairs=all_combination_pairs.to_dict(orient='records'),cancer_type_option_json=cancer_type_option_json,evidence_option_json=evidence_option_json,mutation_option_json=mutation_option_json)




@blue.route('/oncodrug/each_actmut_info/<actmut>')
def each_actmut_info(actmut):
    ####对应models中person()
        ####对应models中person()
    action_mutation_data=original_database[original_database['Action mutations'].str.contains(actmut,case=False).fillna(False)]
    action_mutation_data.columns=pd.Series(action_mutation_data.columns.tolist()).str.replace(' ','_').str.replace('-','_').str.replace('(','').str.replace(')','')

    all_action_mutation_id_all=[]
    for action_mut_eg in action_mutation_data['Action_mutations'].drop_duplicates():
        action_mut_eg_split=action_mut_eg.split('; ')
        all_action_mutation_id_all=np.union1d(action_mut_eg_split,all_action_mutation_id_all)
    all_action_mutation_id_all=pd.Series(all_action_mutation_id_all)
    all_action_mutation_id_all=all_action_mutation_id_all[~(all_action_mutation_id_all.str.contains(actmut,case=False)) ].tolist()
    ##drugbank id###
    link_act_combi=pd.DataFrame()
    link_act_combi['action_mut1']=all_action_mutation_id_all
    link_act_combi['action_mut2']=actmut
    link_act_combi=link_act_combi[link_act_combi['action_mut1'] != link_act_combi['action_mut2']]

####核心节点

    core_node_data={
            'name': f'{actmut}',
            'x': 0,
            'y': 0,'itemStyle': { 'color': '#F4B2B0' } ,'symbolSize': [25, 25]
            }
###其他节点
    #####nodes position####
    if len(all_action_mutation_id_all)==1:
        all_comb_mut_name_drop=all_action_mutation_id_all
    else:
        all_comb_mut_name_drop=np.setdiff1d(all_action_mutation_id_all,actmut)
    
    x_all,y_all=single_drug_pos(len(all_comb_mut_name_drop))
    node_data_nocore=[]
    for pos in range(len(all_comb_mut_name_drop)):
        action_mut_eg=all_comb_mut_name_drop[pos]
        x_eg=x_all[pos]
        y_eg=y_all[pos]
        node_data_eg={
            'name': f'{action_mut_eg}',
            'x': x_eg,
            'y': y_eg
            }
        node_data_nocore=node_data_nocore+[node_data_eg]
    node_data=[core_node_data]+node_data_nocore

    link_data=[]
    for link_pos in range(link_act_combi.shape[0]):
        source_eg=''.join(link_act_combi.iloc[link_pos:(link_pos+1),:]['action_mut1'])
        target_eg=''.join(link_act_combi.iloc[link_pos:(link_pos+1),:]['action_mut2'])
        link_data_eg={
        'source': f'{source_eg}',
        'target': f'{target_eg}',
        'label': {
        'show': 'true','formatter': ''
        },
        'lineStyle': {
        'curveness': 0,
        }
        }
        link_data=link_data+[link_data_eg]
    
    option = {            'title': {
                      'text': 'Biomarkers predicting drug responses','left': 'center','top':25
                      },
                  'animationDurationUpdate': 1500,
                  'animationEasingUpdate': 'quinticInOut',
                  'series': [
                  {
                    'type': 'graph',
                    'layout': 'none',
                    'symbolSize': [25/2, 25/2],
                    'roam': 'true',
                    'label': {
                    'show': 'true','formatter': '{b}'
                    },
                    'edgeLabel': {
                    'fontSize': 20
                    },
                    'itemStyle': {
                    'borderColor': '#000',
                    'borderWidth': 1,
                    'borderType': 'solid',
                    'shadowBlur': 10,
                    'shadowColor': 'rgba(0, 0, 0, 0.3)'},
                    
                    'data': node_data,
                    'links': link_data,
                    'lineStyle': {
                    'opacity': 0.9,
                    'width': 2,'smooth': 'false'
                    },'tooltip': {'show': 'false'},
                  }
                  ]
                };
    # 将 option 对象传递给模板渲染
    option_json = json.dumps(option)
    
    
    ##drug name####
    comb_drug_all=pd.DataFrame()
    for pos in range(action_mutation_data.shape[0]):
        action_mutation_data_eg=action_mutation_data.iloc[pos:(pos+1),:]
        drug_name_eg=action_mutation_data_eg['Targeted_drug'].str.split('; ').iloc[0]+action_mutation_data_eg['Targeted_drug'].str.split('; ').iloc[0]
        drug_name_eg=np.setdiff1d(drug_name_eg,['No info'])
        drug_name_sort_eg=pd.Series(drug_name_eg).drop_duplicates().sort_values()
        drug_name_sort_str_eg='; '.join(drug_name_sort_eg)
        drug_name_sort_str_df_eg=pd.DataFrame([drug_name_sort_str_eg])
        comb_drug_all=pd.concat([drug_name_sort_str_df_eg,comb_drug_all])
    
    all_drug_id_all=pd.DataFrame(comb_drug_all[0].value_counts()).sort_values(0).iloc[0:10,:]

    drug_link_data=[]
    for pos in range(all_drug_id_all.shape[0]):
        all_drug_id_all_eg=all_drug_id_all.iloc[pos:(pos+1),:]
        all_drug_id_all_eg.columns=['count']
        drug_eg=str(all_drug_id_all_eg.index.tolist()[0])
        count_eg=int(all_drug_id_all_eg['count'])
        drug_link_data_eg={
            'name':drug_eg,
            'value':count_eg
        }
        drug_link_data=drug_link_data+[drug_link_data_eg]

        drug_option_json={
            'title': {
                      'text': 'Drugs associated with biomarkers','left': 'center','top':25
                      }
            ,'tooltip': {},
    'series': [ {
    'type': 'wordCloud',
    'gridSize': 2,
    'sizeRange': [10, 30],
    'rotationRange': [-90, 90],
    'shape': 'apple',
    'width': 600,
    'height': 400,
    'drawOutOfBound': 'true',
    'textStyle': {
    'fontFamily': 'sans-serif',
    'fontWeight': 'bold',
        'emphasis': {
        'shadowBlur': 10,
        'shadowColor': '#333'
    }
    },'data':drug_link_data}]
                            };

    drug_option_json = json.dumps(drug_option_json)
    
    all_evidence_id_all=pd.DataFrame(action_mutation_data['Evidence_level'].value_counts())
    all_evidence_id_all.columns=['count']
    all_evidence_id_all['Evidence_level']=all_evidence_id_all.index.tolist()
    all_evidence_id_all.reset_index(drop=True,inplace=True)
    evidence_link_data=[]
    for evidence_level in ['Level A','Level B','Level C','Level D']:
        all_evidence_id_all_eg=all_evidence_id_all[all_evidence_id_all['Evidence_level']==evidence_level]
        if all_evidence_id_all_eg.shape[0]!=0:
            evidence_eg=evidence_level
            count_eg=int(all_evidence_id_all_eg['count'])
            evidence_link_data_eg={
                'name':evidence_eg,
                'value':count_eg
            }
        else:
            evidence_link_data_eg={
                'name':evidence_level,
                'value':0
            }
        evidence_link_data=evidence_link_data+[evidence_link_data_eg]
    evidence_option_json={
      'title': {
      'text': 'Statistical evidence levels for data entries','left': 'center',
      'top': 25,    
      'left': 'center'
      },
          'tooltip': {
          'trigger': 'item'
        },
        'legend': {
          'bottom': '12.5%',
      'left': 'center'
        },
        'series': [
          {
            'name': 'Access From',
            'type': 'pie',
            'radius': ['30%', '60%'],
            'avoidLabelOverlap': 'false',
            'itemStyle': {
              'borderRadius': 10,
              'borderColor': '#fff',
              'borderWidth': 2
            },
            'label': {
              'show': 'false',
              'position': 'center'
            },
            'emphasis': {
              'label': {
                'show': 'true',
                'fontSize': 20,
                'fontWeight': 'bold'
              }
            },
            'labelLine': {
              'show': 'false'
            },
            'data': evidence_link_data
          }
        ]
      };

    evidence_option_json = json.dumps(evidence_option_json)
    
    all_cancer_type_id_all=pd.DataFrame(action_mutation_data['Cancer_type_oncotree_Level2'].value_counts())
    all_cancer_type_id_all=all_cancer_type_id_all.iloc[0:10,:]
    cancer_type_link_data=[]
    for pos in range(all_cancer_type_id_all.shape[0]):
        all_cancer_type_id_all_eg=all_cancer_type_id_all.iloc[pos:(pos+1),:]
        cancer_type_eg=str(all_cancer_type_id_all_eg.index.tolist()[0])
        count_eg=int(all_cancer_type_id_all_eg['count'])
        cancer_type_link_data_eg={
            'name':cancer_type_eg,
            'value':count_eg
        }
        cancer_type_link_data=cancer_type_link_data+[cancer_type_link_data_eg]

        cancer_type_option_json={
            'title': {
                      'text': 'Cancer types associated with biomarkers','left': 'center','top':25
                      }
            ,'tooltip': {},
    'series': [ {
    'type': 'wordCloud',
    'gridSize': 2,
    'sizeRange': [10, 30],
    'rotationRange': [-90, 90],
    'shape': 'apple',
    'width': 600,
    'height': 400,
    'drawOutOfBound': 'true',
    'textStyle': {
    'fontFamily': 'sans-serif',
    'fontWeight': 'bold',
        'emphasis': {
        'shadowBlur': 10,
        'shadowColor': '#333'
    }
    },'data':cancer_type_link_data}]
                            };

    cancer_type_option_json = json.dumps(cancer_type_option_json)
    
    return render_template('each_action_mutation_info.html',actmut=actmut,action_mutation_data=action_mutation_data.to_dict(orient='records'),option_json=option_json,drug_option_json=drug_option_json,evidence_option_json=evidence_option_json,cancer_type_option_json=cancer_type_option_json)

######添加oncotree界面
@blue.route('/oncodrug/each_cancer_oncotree/<cancer_oncotree>')
def each_oncotree_info(cancer_oncotree):
    try:
        ####对应models中person()
            ####对应models中person()
        cancer_oncotree_data=original_database[original_database['Cancer type (oncotree Level2)'].str.contains(cancer_oncotree,case=False).fillna(False)]
        cancer_oncotree_data.columns=pd.Series(cancer_oncotree_data.columns.tolist()).str.replace(' ','_').str.replace('-','_').str.replace('(','').str.replace(')','')
        cancer_oncotree_data.columns.tolist()
        ####action mutation
        all_mutation_id_all=pd.DataFrame(cancer_oncotree_data['Action_mutations'].value_counts()).iloc[0:10,:]
        all_mutation_id_all=all_mutation_id_all.loc[list(~(pd.Series(all_mutation_id_all.index.tolist()).str.contains('Unavailable',case=False))),:]
        mutation_link_data=[]
        for pos in range(all_mutation_id_all.shape[0]):
            all_mutation_id_all_eg=all_mutation_id_all.iloc[pos:(pos+1),:]
            mutation_eg=str(all_mutation_id_all_eg.index.tolist()[0])
            count_eg=int(all_mutation_id_all_eg['count'])
            mutation_link_data_eg={
                'name':mutation_eg,
                'value':count_eg
            }
            mutation_link_data=mutation_link_data+[mutation_link_data_eg]
    
            mutation_option_json={
                'title': {
                          'text': 'Biomarkers associated with cancer type','left': 'center','top': 25, 
                          }
                ,'tooltip': {},
        'series': [ {
        'type': 'wordCloud',
        'gridSize': 2,
        'sizeRange': [10, 30],
        'rotationRange': [-90, 90],
        'shape': 'apple',
        'width': 600,
        'height': 400,
        'drawOutOfBound': 'true',
        'textStyle': {
        'fontFamily': 'sans-serif',
        'fontWeight': 'bold',
            'emphasis': {
            'shadowBlur': 10,
            'shadowColor': '#333'
        }
        },'data':mutation_link_data}]
                                };
    
        mutation_option_json = json.dumps(mutation_option_json)
        ###drug#####
        
        comb_drug_all=pd.DataFrame()
        for pos in range(cancer_oncotree_data.shape[0]):
            cancer_oncotree_data_eg=cancer_oncotree_data.iloc[pos:(pos+1),:]
            cancer_oncotree_data_eg.fillna('No info',inplace=True)
            drug_name_eg=cancer_oncotree_data_eg['Targeted_drug'].str.split('; ').iloc[0]+cancer_oncotree_data_eg['Targeted_drug'].str.split('; ').iloc[0]
            drug_name_eg=np.setdiff1d(drug_name_eg,['No info','nan'])
            drug_name_sort_eg=pd.Series(drug_name_eg).drop_duplicates().sort_values()
            drug_name_sort_str_eg='; '.join(drug_name_sort_eg)
            drug_name_sort_str_df_eg=pd.DataFrame([drug_name_sort_str_eg])
            comb_drug_all=pd.concat([drug_name_sort_str_df_eg,comb_drug_all])
    
        all_drug_id_all=pd.DataFrame(comb_drug_all[0].value_counts()).iloc[0:10,:]
        drug_link_data=[]
        for pos in range(all_drug_id_all.shape[0]):
            all_drug_id_all_eg=all_drug_id_all.iloc[pos:(pos+1),:]
            drug_eg=str(all_drug_id_all_eg.index.tolist()[0])
            count_eg=int(all_drug_id_all_eg['count'])
            drug_link_data_eg={
                'name':drug_eg,
                'value':count_eg
            }
            drug_link_data=drug_link_data+[drug_link_data_eg]
    
            drug_option_json={
                'title': {
                          'text': 'Drugs associated with cancer type','left': 'center','top': 25, 
                          }
                ,'tooltip': {},
        'series': [ {
        'type': 'wordCloud',
        'gridSize': 2,
        'sizeRange': [10, 30],
        'rotationRange': [-90, 90],
        'shape': 'apple',
        'width': 600,
        'height': 400,
        'drawOutOfBound': 'true',
        'textStyle': {
        'fontFamily': 'sans-serif',
        'fontWeight': 'bold',
            'emphasis': {
            'shadowBlur': 10,
            'shadowColor': '#333'
        }
        },'data':drug_link_data}]
                                };
    
        drug_option_json = json.dumps(drug_option_json)
    
        ####evidence####
        all_evidence_id_all=pd.DataFrame(cancer_oncotree_data['Evidence_level'].value_counts())
        all_evidence_id_all.columns=['count']
        all_evidence_id_all['Evidence_level']=all_evidence_id_all.index.tolist()
        all_evidence_id_all.reset_index(drop=True,inplace=True)
        evidence_link_data=[]
        for evidence_level in ['Level A','Level B','Level C','Level D']:
            all_evidence_id_all_eg=all_evidence_id_all[all_evidence_id_all['Evidence_level']==evidence_level]
            if all_evidence_id_all_eg.shape[0]!=0:
                evidence_eg=evidence_level
                count_eg=int(all_evidence_id_all_eg['count'])
                evidence_link_data_eg={
                    'name':evidence_eg,
                    'value':count_eg
                }
            else:
                evidence_link_data_eg={
                    'name':evidence_level,
                    'value':0
                }
            evidence_link_data=evidence_link_data+[evidence_link_data_eg]
        evidence_option_json={
          'title': {
          'text': 'Statistical evidence levels for data entries','left': 'center',
          'top': 25,    
          'left': 'center'
          },
              'tooltip': {
              'trigger': 'item'
            },
            'legend': {
              'bottom': '12.5%',
          'left': 'center'
            },
            'series': [
              {
                'name': 'Access From',
                'type': 'pie',
                'radius': ['30%', '60%'],
                'avoidLabelOverlap': 'false',
                'itemStyle': {
                  'borderRadius': 10,
                  'borderColor': '#fff',
                  'borderWidth': 2
                },
                'label': {
                  'show': 'false',
                  'position': 'center'
                },
                'emphasis': {
                  'label': {
                    'show': 'true',
                    'fontSize': 20,
                    'fontWeight': 'bold'
                  }
                },
                'labelLine': {
                  'show': 'false'
                },
                'data': evidence_link_data
              }
            ]
          };
        evidence_option_json = json.dumps(evidence_option_json)
        
        return render_template('each_cancer_oncotree_info.html',cancer_oncotree=cancer_oncotree,cancer_oncotree_data=cancer_oncotree_data.to_dict(orient='records'),mutation_option_json=mutation_option_json,drug_option_json=drug_option_json,evidence_option_json=evidence_option_json)
    except:
        return render_template('navbar_error.html' )

###添加下载功能
@blue.route('/oncodrug/download_db')
def download_db():
    ####对应models中person()
    return render_template('download_drugdb.html')

@blue.route('/oncodrug/download/<filename>')
def download_db_file(filename):
    # 定义文件路径，您可以根据实际情况进行调整
    file_path = '/h/tianyi/oncodrug/upload/drug_comdb_data/' + filename
    return send_file(file_path, as_attachment=True)


@blue.route('/oncodrug/about_document')
def about_document():
    return render_template('about_document.html')

