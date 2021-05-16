#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:05:28 2020

@author: josephcavataio
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy.ma as ma
import requests
from scipy import stats
from itertools import chain
from itertools import combinations


#### Function for creating groups of entries ####
def organize_index(data):
    uq_counts_1= np.unique(data['count'])
    uq_counts = uq_counts_1[::-1]
    sorted_index = []
    j = 0
    for i in uq_counts:
        a1 = data['ec_number'][data['count'] == i]#data_MM_log10 = data_all_log10.reindex(nw_idx)
        a11 = data['substrate'][data['count'] == i]
        idx_a4 = np.argsort(a11)
        if len(a1) == 1:
            sorted_index.append(j)
            j = j+1
        else:
            k2 = j + idx_a4
            j = max(k2) + 1
            for l2 in k2:
                sorted_index.append(l2)
    return sorted_index  

def do_km_stats(x):
    otpt = pd.Series(x).describe() 
    return otpt 

def select_outs(x,l,u):
    X = np.asarray(x)
    opts = X[(X < l) | (X > u)]
    nopts1 = np.argwhere( (X < l) | (X > u) )
    
    if nopts1.size > 0:
        nopts = np.concatenate(nopts1,axis=0)
    else:
        nopts = nopts1
   
    nints =  np.concatenate(np.argwhere( (X >= l) & (X <= u) ),axis=0)
    return opts,nopts,nints

def data_stats(X,mthd):
    L = len(X)
    dataX = pd.DataFrame([do_km_stats(X["Km"][i]) for i in range(L)])
    
    #Outliers: Quartile definition or One standard deviation away
    if mthd == 'quartile':
        low_lim = dataX['25%'] - 1.5*(dataX['75%']-dataX['25%'])
        up_lim = dataX['75%'] + 1.5*(dataX['75%']-dataX['25%'])
        outs = pd.Series([select_outs(X['Km'][i],low_lim[i],up_lim[i])[0] for i in range(L)])
        idxs_outs = [select_outs(X['Km'][i],low_lim[i],up_lim[i])[1] for i in range(L)]
        idxs_fair = [select_outs(X['Km'][i],low_lim[i],up_lim[i])[2] for i in range(L)]
        len_outs = pd.Series([len(outs[i]) for i in range(L)])
        perc_outs = 100*(len_outs/dataX['count'])
        dataX['LL'] = low_lim
        dataX['UL'] = up_lim
        dataX['%outliers'] = perc_outs        
        
    elif mthd == 'std':
        low_lim = dataX['mean'] - dataX['std']
        up_lim = dataX['mean'] + dataX['std']
    
        outs = pd.Series([select_outs(X['Km'][i],low_lim[i],up_lim[i])[0] for i in range(L)])
        idxs_outs = [select_outs(X['Km'][i],low_lim[i],up_lim[i])[1] for i in range(L)]
        idxs_fair = [select_outs(X['Km'][i],low_lim[i],up_lim[i])[2] for i in range(L)]
        len_outs = pd.Series([len(outs[i]) for i in range(L)])
        perc_outs = 100*(len_outs/dataX['count'])
        dataX['LL'] = low_lim
        dataX['UL'] = up_lim
        dataX['%outliers'] = perc_outs 
        
    else:
        raise Exception('invalid outlier selection method')
        
    #distribution test    
    Q = np.asarray([stats.shapiro(X['Km'][i]) for i in range(L) ])
    pvals = pd.Series(Q[:,1])
    dataX['p-value'] = pvals
    return dataX,idxs_outs,idxs_fair


import matplotlib.ticker as mticker
class MathTextSciFormatter(mticker.Formatter):
    def __init__(self, fmt="%1.2e"):
        self.fmt = fmt
    def __call__(self, x, pos=None):
        s = self.fmt % x
        decimal_point = '.'
        positive_sign = '+'
        tup = s.split('e')
        significand = tup[0].rstrip(decimal_point)
        sign = tup[1][0].replace(positive_sign, '')
        exponent = tup[1][1:].lstrip('0')
        if exponent:
            exponent = '10^{%s%s}' % (sign, exponent)
        if significand and exponent:
            s =  r'%s{\times}%s' % (significand, exponent)
        else:
            s =  r'%s%s' % (significand, exponent)
        return "${}$".format(s)

def make_groups(data,MRC,IsLogTransf,P):
    ## Group the entries according to criteria in the groups variable ##
    data0 = data.groupby(MRC)
    group_index = list(data0.indices)
    
    km_vals = []
    km_dev = []
    entry_id = []
    formulas = []

    L = len(group_index)
    print('%d is the no of sets under MRCs.'%L)
    for i in range(L):
        example = data0.get_group(group_index[i])
        entry_id.append(list(example['entry_id']))
        formulas.append(list(example['kineticlaw_formula']))
        
        if IsLogTransf == True:
            km_vals.append(list(example['Km log10']))
        else:
            km_vals.append(list(example['km_start_value']))
            
        if np.all(np.unique(list(example['km_deviation'])) == '-'):
            km_dev.append(['-'])
        elif np.all(np.unique(list(example['km_deviation'])) == '(null)'):
            km_dev.append(['-'])
        else:
            km_dev.append(list(example['km_deviation']))
    
    grp = pd.DataFrame({'count': data0.size()}).reset_index()
    grp['Km'] = pd.Series(km_vals)
    grp['Km dev'] = pd.Series(km_dev)
    grp['ec_id'] = pd.Series(entry_id)
    grp['formula'] = pd.Series(formulas)
    
    grp  = grp.reindex(np.argsort(grp['count'])[::-1])
    grp  = grp.reset_index(drop=True)
    dat = grp[grp['count']>=P]
    L2 = len(dat)
    print('{} is the no of sets under MRCs with {} or more entries.\n'.format(L2,P))

    # Record only one formula for each set
    ind_form = dat.columns.get_loc('formula')

    for i in range(L2):
        # Lf = len(np.unique(dat.iloc[(i,ind_form)]))
        # if Lf > 1:
        #     dat.iloc[i,ind_form] = [[] for _ in range(Lf)]
        #     for j in range(Lf):
        #         dat.iloc[i,ind_form][j] = np.unique(dat.iloc[(i,ind_form)])[j]
        # else:
        dat.iloc[i,ind_form] = np.unique(dat.iloc[(i,ind_form)])[0]          
    return dat 


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
print(round(100*len(data[(buffer=='not available') | (buffer=='-')])/len(data),2), '\n') #not available & -
buffer_uniq = pd.DataFrame(np.unique(data['buffer'],return_counts=True))

print("Percent of entries missing organism")
print(round(100*len(data[(org=='unknown') | (org=='-')])/len(data),2))

#%%# Clean-up: Drop km end values and deviation from entries ####  
data = data.drop(['km_end_value'],axis=1)
#data = data.drop(['km_deviation'],axis=1)


#### Clean-up: Sort entries by entry_id ####
data = data.set_index('entry_id', drop=False)
data = data.sort_index(0) #There are duplicate entry IDs for sequential ordered bi-bi entries
L0 = len(data)


#### Clean-up: Units (include only M) ####
units = list(data["unit"].unique())
f_un_units = [len(data[data['unit']==units[i]]) for i in range(len(units))]
data_cor_units = data[data['unit']==units[0]]

#select NM mechanisms and look at the distribution of the year
lawtypes = pd.DataFrame(np.unique(data['kineticlaw_type'],return_counts=True))
data_MM_pregroup = data_cor_units[data_cor_units['kineticlaw_type']=='Michaelis-Menten']


IDX_1 = np.where(data_MM_pregroup['km_start_value'] == 0)[0] 
IDX_2 = np.where(data_MM_pregroup['km_start_value'] < 0)[0]
IDX_3 = np.concatenate([IDX_1,IDX_2])


ENTRIES_BAD = [data_MM_pregroup['entry_id'].iloc[i] for i in IDX_3]
KM_BAD = [data_MM_pregroup['km_start_value'].iloc[i] for i in IDX_3]
bad_entries = pd.DataFrame(list(zip(ENTRIES_BAD,KM_BAD)),columns=['entry id','Km'])
bad_entries.to_csv('bad_entries.csv')

data_MM_pregroup=data_MM_pregroup[data_MM_pregroup['km_start_value']>0]
data_MM_pregroup=data_MM_pregroup.reset_index(drop=True)

#%%# Intro Graph ####
ind_km = data_MM_pregroup.columns.get_loc('km_start_value')
values = data_MM_pregroup.iloc[:,ind_km].dropna(axis=0)
km_uniq, km_freq = np.unique(values,return_counts = True)
uniq_km = pd.DataFrame(list(zip(km_uniq,km_freq)),columns=['values','frequency'])

ind_zero = values.index[values == 0].tolist()
totstats = np.log10(values.drop(ind_zero)).describe()

plt.figure(figsize=(16,12))
plt.scatter(np.log10(uniq_km['values']),uniq_km['frequency'])
plt.plot([totstats['mean']]*31,range(0,155,5),color = 'r',label = '{} ± {}'.format(round(totstats['mean'],2),round(totstats['std'],2)))
plt.yticks(range(0,175,25), fontsize = 16)
plt.xticks(range(-11,6,1), fontsize = 16)
plt.xlim(-11,2)
plt.xlabel('$K_M$ (log$_{10}$)', fontsize = 16)
plt.ylabel('Frequency', fontsize = 16)
plt.title('Log$_{10}$ Histogram of Entire Dataset Using Michaelis-Menten Law Type', fontsize = 24, y=1.01)
legend = plt.legend(bbox_to_anchor=(1.0, 1), loc='upper right', borderaxespad=0.8,\
           title='Mean ± Stdev',fontsize = 16)
plt.rcParams['legend.title_fontsize'] = 'xx-small'

#%%# Group entries by the following groups ####
groups = ['ec_number','kineticlaw_type','substrate','organism','ph','temperature','buffer'] 

data_MM_pregroup['Km log10'] = np.nan
ind_km = data_MM_pregroup.columns.get_loc('km_start_value')
ind_kmlog10 = data_MM_pregroup.columns.get_loc('Km log10')

for i in range(len(data_MM_pregroup)):
    data_MM_pregroup.iloc[i,ind_kmlog10] = np.log10(data_MM_pregroup.iloc[i,ind_km])

data_all = make_groups(data_MM_pregroup,groups,False,10)
data_all=data_all.reset_index(drop=True)

data_all_log10 = make_groups(data_MM_pregroup,groups,True,10)
data_all_log10=data_all_log10.reset_index(drop=True)

nw_idx = organize_index(data_all)
data_MM = data_all.reindex(nw_idx)
data_MM=data_MM.reset_index(drop=True)

data_MM_log10 = data_all_log10.reindex(nw_idx)
data_MM_log10 = data_MM_log10.reset_index(drop=True)

#data_MM.to_csv('data_MM.csv')

#%%# Descriptive statistics for MM ####
L_MM = len(data_MM)
stats_MM = data_stats(data_MM,'std')
stats_data_MM = stats_MM[0]
out_idx_MM_std = stats_MM[1]
in_idx_MM_std = stats_MM[2]


stats_MM_log10 = data_stats(data_MM_log10,'std')
stats_data_MM_log10 = stats_MM_log10[0]
out_idx_MM_std_log10 = stats_MM_log10[1]
in_idx_MM_std_log10 = stats_MM_log10[2]


normMM = len(stats_data_MM['p-value'][stats_data_MM['p-value']>0.05])
normMM_log10 = len(stats_data_MM_log10['p-value'][stats_data_MM_log10['p-value']>0.05])

print('Estimated percentage of normally-distributed sets: %.2f' %(np.round(100*normMM/L_MM,2)))
print('Estimated percentage of log-normally-distributed sets: %.2f \n' %(np.round(100*normMM_log10/L_MM,2)))

### Getting entry ID's of outliers and inliers ###

ec_id_out_list_std = []
ec_id_in_list_std = []

buff_out_list_std = []
buff_in_list_std = []

ind_buff = data_MM.columns.get_loc('buffer')
ind_count = data_MM.columns.get_loc('count')
ind_ecid = data_MM_log10.columns.get_loc('ec_id')

for i in range(L_MM):
    if len(out_idx_MM_std_log10[i]) > 0:

        entry_ids_out_std = [data_MM_log10.iloc[i,ind_ecid][j] for j in out_idx_MM_std_log10[i]]
        entry_ids_in_std = [data_MM_log10.iloc[i,ind_ecid][j] for j in in_idx_MM_std_log10[i]]
        buf_out_std = [data_MM.iloc[i,ind_buff]] * len(out_idx_MM_std_log10[i])
        buf_in_std = [data_MM.iloc[i,ind_buff]] * (data_MM['count'][i]-len(out_idx_MM_std_log10[i]))
        

        ec_id_out_list_std.append(entry_ids_out_std)
        ec_id_in_list_std.append(entry_ids_in_std)  
        buff_out_list_std.append(buf_out_std)
        buff_in_list_std.append(buf_in_std)
    else:

        buf_in_std = [data_MM.iloc[i,ind_buff]] * data_MM.iloc[i,ind_count]
        
        ec_id_out_list_std.append([])
        ec_id_in_list_std.append(data_MM['ec_id'][i])
        
        buff_out_list_std.append([])
        buff_in_list_std.append(buf_in_std)


all_out_entryids_std = [y for x in ec_id_out_list_std for y in x]
all_in_entryids_std = [y for x in ec_id_in_list_std for y in x]
all_out_buff_std = [y for x in buff_out_list_std for y in x]
all_in_buff_std = [y for x in buff_in_list_std for y in x]

all_out_std = pd.DataFrame(list(zip(all_out_entryids_std,all_out_buff_std)),columns=['entry id','buffer'])
all_in_std = pd.DataFrame(list(zip(all_in_entryids_std,all_in_buff_std)),columns=['entry id','buffer'])

km_out_std = all_out_std
km_in_std = all_in_std

#km_out_std.to_csv('km_out_std.csv')
#km_in_std.to_csv('km_in_std.csv')

#%%$ Set analysis for range vs. publication data  ####
maxmin_MM1 = np.log10(stats_data_MM["max"]/stats_data_MM["min"])
maxmin_MM = np.asarray(maxmin_MM1)
proc_OK_MM = np.round(100*len(maxmin_MM[maxmin_MM<=1])/L_MM,2)
proc_nOK_MM  = 100 - proc_OK_MM
maxmin_MM2 = np.round(maxmin_MM,4)

fdata_MM = pd.DataFrame(list(zip(data_MM['count'],data_MM['ec_id'],data_MM['Km'],
                data_MM_log10['Km'],data_MM['ec_number'],maxmin_MM2)),
                columns = ['count','entry id','Km','Km (log10)','ec_number','range'])
#fdata_MM.to_csv('fdata_MM.csv')
data_MM['range'] = maxmin_MM2


#%%# For distribution based on years ####
MM_years = data_MM_pregroup['year']
unique_MM_years,counts_MM_years=np.unique(MM_years,return_counts = 'True')
L_MM_years = len(unique_MM_years)
max_years = np.max(counts_MM_years) 
idx_max_years = np.where(counts_MM_years==max_years)[0][0]
what_year_max = unique_MM_years[idx_max_years]

plt.figure(1)
plt.figure(figsize=(11,4), constrained_layout=False) ## constrained_layout=True
plt.rcParams.update({'font.size':22})
plt.bar(unique_MM_years[1:-1],counts_MM_years[1:-1],0.75,color='olive',align='edge')
plt.yticks([0,500,1000])
plt.xlim(1955,2020)
plt.xlabel('Publication year',labelpad=10)
plt.ylabel('No. publications',labelpad=10)
plt.title('All publications', fontsize=22)
#plt.savefig('figures_new/keep_the_title.pdf',bbox_inches='tight')

#%%# PLOTS FOR THE FLOWCHART ####

NOW = 3 # Which set to do an analysis on 
ind_km = data_MM.columns.get_loc('Km')
ind_kmlog10 = data_MM_log10.columns.get_loc('Km')

A = data_MM.iloc[NOW,ind_km]
B = data_MM_log10.iloc[NOW,ind_kmlog10]


def gauss_func(X):  
    stdX = np.std(X)
    meanX = np.mean(X)
        
    up  = np.mean(X)+np.std(X)
    dwn = np.mean(X)-np.std(X)
    
    rng_max = np.mean(X)+3.25*np.std(X)
    rng_min = np.mean(X)-3.25*np.std(X)
    x = np.linspace(rng_min,rng_max,100)
    out = (1/stdX * np.sqrt(np.pi)) * np.exp(-(((x - meanX)/stdX)**2)/2)
    Q = [out,x,dwn,up]
    return Q

q = gauss_func(B)


plt.figure(3)
plt.rcParams.update({'font.size':28})
plt.hist(A,16,color='black',label='N = %d'%len(A))
plt.xlabel('$K_M$ [M]',labelpad=10)
plt.ylabel('Frequency',labelpad=10) 
#plt.ylim(0,8)
#plt.xlim(-0.00025,0.0045)
#plt.yticks([0,4,8])
#plt.xticks([0,0.002,0.004])
plt.legend(loc='best',frameon=False)
plt.gca().xaxis.set_major_formatter(MathTextSciFormatter("%1.0e"))
#plt.savefig('figures/hist_wannabe_28.pdf',bbox_inches='tight')

IDX_LOW = np.where(q[1] == max(q[1][q[1]<q[2]]))[0][0]
IDX_HIGH = np.where(q[1] == min(q[1][q[1]>q[3]]))[0][0]

part1 = q[1][0:(IDX_LOW+1)]
y_part1 = q[0][0:(IDX_LOW+1)]

part2 = q[1][IDX_LOW:(IDX_HIGH+1)]
y_part2 = q[0][IDX_LOW:(IDX_HIGH+1)]

part3 = q[1][IDX_HIGH:-1]
y_part3 = q[0][IDX_HIGH:-1]

plt.figure(4)
plt.rcParams.update({'font.size':28})
plt.plot(q[1],q[0],color='black',linewidth=3)
plt.hist(B,16,color='black')
plt.fill_between(part1,y_part1,color='red',alpha=0.5)
plt.fill_between(part2,y_part2,color='forestgreen',alpha=0.5)
plt.fill_between(part3,y_part3,color='red',alpha=0.5)
plt.xlim(min(q[1]),max(q[1]))
plt.yticks([0,3,6])
plt.ylim(0,1+int(2*max(q[0]))) 
plt.xlabel('log$_{10}$ $K_M$ [M]',labelpad=10)
plt.ylabel('Frequency',labelpad=10)
#plt.title('N = %d' %len(A), fontsize=22)
#plt.savefig('figures_flowchart/logtransf_hist_4.pdf',bbox_inches='tight')

#%%# Plotting range for each set ####
fig = plt.figure(figsize=(11,8))
ax1 = plt.subplot2grid((1,7),(0,0),rowspan=1,colspan=6)
ax2 = plt.subplot2grid((1,7),(0,6),rowspan=1,colspan=1)
#plt.subplots(ax1,ax2,sharey = True,facecolor = 'w')

#Set ticks
lower_xticks = np.arange(10,32,2)
ax1.set_xticks(lower_xticks)
upper_xticks = np.arange(61,63,1)
ax2.set_xticks(upper_xticks)

#Plot on each subplot
#ax1.plot(data_MM_no_group['count'],maxmin_MM_no_group,'ro',markersize=12,\
         #label = 'No Buffer')
ax1.plot(data_MM['count'],maxmin_MM2,'bo',markersize=12,label = 'All MRCs')
#ax2.plot(data_MM_no_group['count'],maxmin_MM_no_group,'ro',markersize=12,\
         #label = 'No Buffer')
ax2.plot(data_MM['count'],maxmin_MM2,'bo',markersize=12,label = 'All MRCs')

ax1.set_xlim(8,32)
ax2.set_xlim(61,63)
ax1.set_ylim(0,4)
ax2.set_ylim(0,4)

#Plot horizontal line
linex = np.arange(8,64,1)
liney = [1]*56
ax1.plot(linex,liney,color='red',linestyle='--',linewidth=1)
ax2.plot(linex,liney,color='red',linestyle='--',linewidth=1)

#Formatting
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax1.yaxis.tick_left()
ax2.set_yticks([])
#Adjust space btw subplots
plt.subplots_adjust(wspace=0.1)

# Make break clips
d = .02
kwargs = dict(transform=ax1.transAxes, color='r', clip_on=False)
ax1.plot((1-d,1+d), (-d,+1.5*d), **kwargs)
ax1.plot((1-d,1+d),(1-d,1+1.5*d), **kwargs)

kwargs.update(transform=ax2.transAxes)
ax2.plot((-5*d,+5*d), (-d,+d), **kwargs)   
ax2.plot((-5*d,+7*d), (1-d,1+1.5*d), **kwargs)

# Plot % of each and labels
plt.text(63,2.1,'$\leftarrow$ %.0f %%' %proc_nOK_MM)
plt.text(63,0.4,'$\leftarrow$ %.0f %%' %proc_OK_MM)
ax1.set_ylabel('$K_M$ range (log$_{10}$)',labelpad=15)
ax1.set_xlabel('No. replicates',labelpad=20)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.8,\
           fontsize = 14)
#plt.xlabel("No. replicates",labelpad=25)
#plt.axhline(y=1,color='red',linestyle='--',linewidth=1)
plt.savefig('figures_new/log10_range_MM_new.pdf',bbox_inches='tight') 

#%%# Graphing Entry ID's ####
L = len(data_MM)
ind = range(L)

lists_ec_id = [[] for _ in range(L)]
ind_ecid = data_MM.columns.get_loc('ec_id')
plt.figure(figsize=(20,10), constrained_layout=False)

for i in range(L):
    lists_ec_id[i] = data_MM.iloc[(i,ind_ecid)]
    x = [i]*len(lists_ec_id[i])
    plt.plot(x,lists_ec_id[i], linewidth=5)
        

plt.ylabel('Entry ID Number',fontsize=16)
plt.xlabel('Sets',fontsize=16)
plt.xticks(ind, fontsize=14)
plt.yticks(np.arange(4000, 60000, 5000),size=14)
plt.title("Entry ID's for Each Set",fontsize=24,y=1.02)
plt.show()

#%%# Find out how many publications are in each group ####
QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv'
PARAM_QUERY_URL = 'http://sabiork.h-its.org/entry/exportToExcelCustomizable'

title_lists_pub = [[] for _ in range(L)]
auth_lists_pub = [[] for _ in range(L)]
pub_count = [0]*L
count = 0

# Other possible fields: 'EntryID','Enzymename','Substrate','pH','Temperature'
query1 = {'format':'tsv','fields[]':['Title']} 
query2 = {'format':'tsv','fields[]':['Author']} 

for j in range(2):
    for i in range(L):
        if j == 0:
            data_field= {"entryIDs[]":lists_ec_id[i],'EnzymeType':'"wildtype"'}
                        
            request = requests.post(PARAM_QUERY_URL, params=query1, data=data_field)
            request.raise_for_status()
            req = request.text
        
            req = req[6:] # Eliminate 'Title' at start of string
            req = req[:-1] # Eliminate the space at end of string
            req = req.split('\n') # Separate string by rows
            
            req_uniq = np.unique(req,return_counts = True)
            title_lists_pub[i] = req_uniq[0]
            pub_count[i] = req_uniq[1]
            #print(req_uniq[0])
            
            if len(req_uniq[1]) > 1:
                count = count + 1
        
        if j == 1:
            data_field= {"entryIDs[]":lists_ec_id[i],'EnzymeType':'"wildtype"'}
                        
            request = requests.post(PARAM_QUERY_URL, params=query2, data=data_field)
            request.raise_for_status()
            req = request.text
        
            req = req[7:] # Eliminate 'Title' at start of string
            req = req[:-1] # Eliminate the space at end of string
            req = req.split('\n') # Separate string by rows
            
            req_uniq = np.unique(req)
            auth_lists_pub[i] = req_uniq
            #print(req_uniq[0])
        

data_MM['Publication Titles'] = title_lists_pub
data_MM['Pub Count'] = pub_count
data_MM['Authors'] = auth_lists_pub

#fdata_MM['Publication Titles'] = title_lists_pub
#fdata_MM['Pub Count'] = pub_count

print('\n{} of {} Sets have entries from more than one publication'.format(count,L))

#%%###### Correlation between type of study and range of Km values ###########
ind_pub = data_MM.columns.get_loc('Publication Titles')
pubs = []
for i in range(len(data_MM)):
    pubs.append(data_MM.iloc[i,ind_pub][0])
pub_uniq, pub_freq= np.unique(pubs,return_counts = True)
uniq_pubs = pd.DataFrame(list(zip(pub_uniq,pub_freq)),columns=['names','frequency'])

#%%# Plotting pairs of Measured Values ####
#### Control graph ####
import random

data_MM_pregroup['Km log10'] = np.nan
ind_km = data_MM_pregroup.columns.get_loc('km_start_value')
ind_kmlog10 = data_MM_pregroup.columns.get_loc('Km log10')

for i in range(len(data_MM_pregroup)):
    data_MM_pregroup.iloc[i,ind_kmlog10] = -1*np.log10(data_MM_pregroup.iloc[i,ind_km])

# Create data for pair analysis
P=3
data_all2 = make_groups(data_MM_pregroup,groups,True,P)
data_all2=data_all2.reset_index(drop=True)

nw_idx2 = organize_index(data_all2)
data_MM2 = data_all2.reindex(nw_idx2)
data_MM2=data_MM2.reset_index(drop=True)

ind_km = data_MM2.columns.get_loc('Km')
Ltot = len(data_MM2)
km_list1 = [0]*Ltot
km_list2 = [0]*Ltot

for i in range(Ltot):
    temp = random.sample(data_MM2.iloc[i,ind_km],2)
    km_list1[i] = temp[0]
    km_list2[i] = temp[1]
    
df_pairs = pd.DataFrame({'value 1': km_list1, 'value 2': km_list2})

plt.scatter(df_pairs['value 1'],df_pairs['value 2'],color='k',alpha=0.3,edgecolor=None)
plt.xticks(range(0,12,2), fontsize=12)
plt.yticks(range(0,12,2), fontsize=12)
plt.xlabel('$pK_{M1}$',fontsize=12)
plt.ylabel('$pK_{M2}$',fontsize=12)
plt.title('Plot of Pairs of Km Values Same enzyme&sub (min {}//{} sets)'.format(P,Ltot),
          fontsize = 16, y=1.02)
plt.savefig('figures/Pair_plot_3vals_all_MRCs.pdf',bbox_inches='tight')

#%%# Pairs of Values with Sets regardless of number of publications ####
data_MM_pregroup['Km log10'] = np.nan
ind_km = data_MM_pregroup.columns.get_loc('km_start_value')
ind_kmlog10 = data_MM_pregroup.columns.get_loc('Km log10')

for i in range(len(data_MM_pregroup)):
    data_MM_pregroup.iloc[i,ind_kmlog10] = -1*np.log10(data_MM_pregroup.iloc[i,ind_km])

# Create data for pair analysis
P=3
data_all2 = make_groups(data_MM_pregroup,groups,True,P)
data_all2=data_all2.reset_index(drop=True)

nw_idx2 = organize_index(data_all2)
data_MM2 = data_all2.reindex(nw_idx2)
data_MM2=data_MM2.reset_index(drop=True)

ind_km = data_MM2.columns.get_loc('Km')
Ltot = len(data_MM2)
km_list1 = [0]*Ltot
km_list2 = [0]*Ltot

for i in range(Ltot):
    temp = random.sample(data_MM2.iloc[i,ind_km],2)
    km_list1[i] = temp[0]
    km_list2[i] = temp[1]
    
df_pairs = pd.DataFrame({'value 1': km_list1, 'value 2': km_list2})

plt.scatter(df_pairs['value 1'],df_pairs['value 2'],color='k',alpha=0.3,edgecolor=None)
plt.xticks(range(0,12,2), fontsize=12)
plt.yticks(range(0,12,2), fontsize=12)
plt.xlabel('$pK_{M1}$',fontsize=12)
plt.ylabel('$pK_{M2}$',fontsize=12)
plt.title('Plot of Pairs of Km Values Same enzyme&sub (min {}//{} sets)'.format(P,Ltot),
          fontsize = 16, y=1.02)
plt.savefig('figures/Pair_plot_3vals_all_MRCs.pdf',bbox_inches='tight')

#%%## Find out publications in sets for larger dataset and plot graph####
title_lists_pub2 = [[] for _ in range(Ltot)]
pub_count2 = [0]*Ltot
count = 0

# Other possible fields: 'EntryID','Enzymename','Substrate','pH','Temperature'
lists_ec_id2 = [[] for _ in range(Ltot)]
ind_ecid = data_MM2.columns.get_loc('ec_id')
multpubs = []

for i in range(Ltot):
    lists_ec_id2[i] = data_MM2.iloc[(i,ind_ecid)]
    data_field= {"entryIDs[]":lists_ec_id2[i],'EnzymeType':'"wildtype"'}
                
    request = requests.post(PARAM_QUERY_URL, params=query1, data=data_field)
    request.raise_for_status()
    req = request.text

    req = req[6:] # Eliminate 'Title' at start of string
    req = req[:-1] # Eliminate the space at end of string
    req = req.split('\n') # Separate string by rows
    
    req_uniq = np.unique(req,return_counts = True)
    title_lists_pub2[i] = req_uniq[0]
    pub_count2[i] = req_uniq[1]
    #print(req_uniq[0])
    
    if len(req_uniq[1]) > 1:
        count = count + 1
        multpubs.append(i)

data_MM2['Publication Titles'] = title_lists_pub2
data_MM2['Pub Count'] = pub_count2
print('\n{} of {} Sets have entries from more than one publication'.format(count,Ltot))
print('\n{}% of sets in this analysis have replicates from the same paper'.format(round((Ltot-count)/Ltot,2)*100))

df_multpubs = pd.DataFrame(columns = data_MM2.columns, index = multpubs) 
for i in range(len(multpubs)):
    df_multpubs.iloc[i,:] = data_MM2.iloc[multpubs[i],:]

ind_km = df_multpubs.columns.get_loc('Km')
L2 = len(df_multpubs)
km_list3 = [0]*L2
km_list4 = [0]*L2

for i in range(L2):
    temp = random.sample(df_multpubs.iloc[i,ind_km],2)
    km_list3[i] = temp[0]
    km_list4[i] = temp[1]
    
df_pairs_2 = pd.DataFrame({'value 1': km_list3, 'value 2': km_list4})

plt.scatter(df_pairs_2['value 1'],df_pairs_2['value 2'],color='k',alpha=0.3,edgecolor=None)
plt.xticks(range(0,10,2), fontsize=12)
plt.yticks(range(0,10,2), fontsize=12)
plt.xlabel('$pK_{M1}$',fontsize=12)
plt.ylabel('$pK_{M2}$',fontsize=12)
plt.title('Plot of Pairs of Km Values Under Same MRCs in Sets with >1 Publication ({} sets)'.format(count),
          fontsize = 16, y=1.02)
plt.savefig('figures/Pair_plot_3vals_>1pubs_all_MRCs.pdf',bbox_inches='tight')

