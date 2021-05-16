#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 11:44:13 2021

@author: josephcavataio
"""
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import requests

#%%# Read in Km Data ####
data = pd.read_csv('km_data_corrected.txt', sep='\t')

#skiprows=1,names=["entry_id","year","ec_number","organism","reaction_id",\
    #"reaction_equation","km_start_value","km_end_value","km_deviation",\
    #"unit","ph","temperature","buffer","kineticlaw_type",\
    #"kineticlaw_formula","substrate"])

#### Percent of data excluded for different categories ####
kmend = data['km_end_value']
devs =  data['km_deviation']
klt = data['kineticlaw_type']
klf = data['kineticlaw_formula']
buffer = data['buffer']
org = data['organism']

print("Percent of entries missing km_end_value")
print(round(100*len(kmend[kmend=='(null)'])/len(data),2))

print("Percent of entries missing km_deviation")
print(round(100*len(devs[(devs=='(null)') | (devs=='-')])/len(data),2))

print("Percent of entries missing kineticlaw_type")
print(round(100*len(data[(klt=='unknown') | (klt=='-')])/len(data),2))

print("Percent of entries missing kineticlaw_formula")
print(round(100*len(data[(klf=='unknown') | (klf=='-')])/len(data),2))

print("Percent of entries missing buffer")
print(round(100*len(data[(buffer=='not available') | (buffer=='-')])/len(data),2)) #not available & -
buffer_uniq = pd.DataFrame(np.unique(data['buffer'],return_counts=True)).transpose()

print("Percent of entries missing organism")
print(round(100*len(data[(org=='unknown') | (org=='-')])/len(data),2), '\n')
#%%# Clean-up: Drop km end values and deviation from entries ####  
data = data.drop(['km_end_value'],axis=1)
#data = data.drop(['km_deviation'],axis=1)

#### Clean-up: Sort entries by entry_id ####
data = data.set_index('entry_id', drop=False)
data = data.sort_index(0) #There are duplicate entry IDs for sequential ordered bi-bi entries
L0 = len(data)

# Number of measurements before Unit Clean Up
print('Number of measurements before Unit Clean Up = {}'.format(len(data)))

#### Clean-up: Units (include only M) ####
unit_uniq = pd.DataFrame(np.unique(data['unit'],return_counts=True)).transpose()
data_cor_units = data[data['unit']==unit_uniq[0]]

# Number of measurements after Unit Clean Up
print('Number of measurements after Unit Clean Up = {}'.format(len(data_cor_units)))

# select MM mechanisms
lawtypes = pd.DataFrame(np.unique(data['kineticlaw_type'],return_counts=True))
data_MM_pregroup = data_cor_units#[data_cor_units['kineticlaw_type']=='Michaelis-Menten']

# Number of measurements after Selecting for Michaelis Menten Mechanism
print('Number of measurements after Selecting for Michaelis Menten Mechanism = {}'.format(len(data_MM_pregroup)))

IDX_1 = np.where(data_MM_pregroup['km_start_value'] == 0)[0] 
IDX_2 = np.where(data_MM_pregroup['km_start_value'] < 0)[0]
IDX_3 = np.concatenate([IDX_1,IDX_2])

ENTRIES_BAD = [data_MM_pregroup['entry_id'].iloc[i] for i in IDX_3]
KM_BAD = [data_MM_pregroup['km_start_value'].iloc[i] for i in IDX_3]
bad_entries = pd.DataFrame(list(zip(ENTRIES_BAD,KM_BAD)),columns=['entry id','Km'])
#bad_entries.to_csv('bad_entries.csv')

data_MM_pregroup=data_MM_pregroup[data_MM_pregroup['km_start_value']>0]
data_MM_pregroup=data_MM_pregroup.reset_index(drop=True)

# Create column with pKm values
data_MM_pregroup['pKm'] = np.nan
ind_km_start = data_MM_pregroup.columns.get_loc('km_start_value')
ind_pkm = data_MM_pregroup.columns.get_loc('pKm')

Ltot = len(data_MM_pregroup)
for i in range(Ltot):
    data_MM_pregroup.iloc[i,ind_pkm] = -1*np.log10(data_MM_pregroup.iloc[i,ind_km_start])

# Number of measurements after removing negative measurements
print('Number of measurements after removing negative measurements = {}'.format(len(data_MM_pregroup)))

#%%# Get publications for all MM entries and Save as .csv ####
QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv'
PARAM_QUERY_URL = 'http://sabiork.h-its.org/entry/exportToExcelCustomizable'

data_MM_pregroup['Publication Title'] = np.nan
ind_pub = data_MM_pregroup.columns.get_loc('Publication Title')
ind_enid = data_MM_pregroup.columns.get_loc('entry_id')
count = 0

while Ltot != count:
    if Ltot-count > 100:
        data_field= {"entryIDs[]":data_MM_pregroup.iloc[count:count+100,ind_enid],'EnzymeType':'"wildtype"'}
        
        query = {'format':'tsv','fields[]':['Title']} 
        
        # Other possible fields: 'EntryID','Enzymename','Substrate','pH','Temperature'
        
        request = requests.post(PARAM_QUERY_URL, params=query, data=data_field)
        request.raise_for_status()
        req = request.text
        
        req = req[6:] # Eliminate 'Title' at start of string
        req = req[:-1] # Eliminate the space at end of string
        req = req.split('\n') # Separate string by rows
        
        for i in range(len(req)):
            data_MM_pregroup.iloc[count+i,ind_pub] = pd.Series(req).iloc[i]
        
        count += len(req)
        
    elif Ltot-count<=10 and count != Ltot:
        data_field= {"entryIDs[]":data_MM_pregroup.iloc[count:,ind_enid],'EnzymeType':'"wildtype"'}
        
        query = {'format':'tsv','fields[]':['Title']} 
        
        # Other possible fields: 'EntryID','Enzymename','Substrate','pH','Temperature'
        
        request = requests.post(PARAM_QUERY_URL, params=query, data=data_field)
        request.raise_for_status()
        req = request.text
        
        req = req[6:] # Eliminate 'Title' at start of string
        req = req[:-1] # Eliminate the space at end of string
        req = req.split('\n') # Separate string by rows
        
        for i in range(len(req)):
            data_MM_pregroup.iloc[count+i,ind_pub] = pd.Series(req)[i]
        count += len(req)
    
    del req,data_field

# data_MM_pregroup.to_csv('data_MM_pregroup_pubs.csv')

#%% Retrieve Authors
ind_auth = data_MM_pregroup.columns.get_loc('Publication Author')

while Ltot != count:
    if Ltot-count > 100:
        data_field= {"entryIDs[]":data_MM_pregroup.iloc[count:count+100,ind_enid],'EnzymeType':'"wildtype"'}
        
        query = {'format':'tsv','fields[]':['Author']} 
        
        # Other possible fields: 'EntryID','Enzymename','Substrate','pH','Temperature'
        
        request = requests.post(PARAM_QUERY_URL, params=query, data=data_field)
        request.raise_for_status()
        req = request.text
        
        req = req[7:] # Eliminate 'Title' at start of string
        req = req[:-1] # Eliminate the space at end of string
        req = req.split('\n') # Separate string by rows
        
        for i in range(len(req)):
            data_MM_pregroup.iloc[count+i,ind_auth] = pd.Series(req).iloc[i]
        
        count += len(req)
        
    elif Ltot-count<=10 and count != Ltot:
        data_field= {"entryIDs[]":data_MM_pregroup.iloc[count:,ind_enid],'EnzymeType':'"wildtype"'}
        
        query = {'format':'tsv','fields[]':['Title']} 
        
        # Other possible fields: 'EntryID','Enzymename','Substrate','pH','Temperature'
        
        request = requests.post(PARAM_QUERY_URL, params=query, data=data_field)
        request.raise_for_status()
        req = request.text
        
        req = req[7:] # Eliminate 'Title' at start of string
        req = req[:-1] # Eliminate the space at end of string
        req = req.split('\n') # Separate string by rows
        
        for i in range(len(req)):
            data_MM_pregroup.iloc[count+i,ind_auth] = pd.Series(req)[i]
        count += len(req)
    
    del req,data_field

data_MM_pregroup.to_csv('data_MM_pregroup_pubs2.csv')

# Create column with only first author
count=0
data_MM_pregroup['First Author'] = ''
ind_1a = data_MM_pregroup.columns.get_loc('First Author')
ind_a = data_MM_pregroup.columns.get_loc('Publication Author')
for i in range(len(data_MM_pregroup)):
    temp = data_MM_pregroup.iloc[i,ind_a]
    if type(temp) == float:
        data_MM_pregroup.iloc[i,ind_1a] = temp
    else:
        count = count+1
        ind = temp.find(',')
        data_MM_pregroup.iloc[i,ind_1a] = temp[0:ind]
        
# data_MM_pregroup.to_csv('data_MM_pregroup_pubs_auth.csv')
