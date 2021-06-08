#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 13:01:54 2021

@author: josephcavataio
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:05:28 2020

@author: josephcavataio
"""

# Add dependencies so that packages automatically download
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from functions import organize_index, mean_confidence_interval,sortbypubs, \
    getallpairs, residuals, measuresofquality,randomvals_and_diff, get_pairs_tot, \
        pairplot, EnzymeClassAnalysis, findpubs, findfirstauthors, make_groups

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

#%% Load in Pre-saved 'data_MM_pregroup' data ###
data_MM_pregroup = pd.read_csv('data_MM_pregroup_pubs_auth.csv',index_col=0)
Ltot = len(data_MM_pregroup)
ind_pkm = data_MM_pregroup.columns.get_loc('pKm')

# Create a column for the enzyme class
data_MM_pregroup['Enzyme Class'] = ''
ind_EC = data_MM_pregroup.columns.get_loc('Enzyme Class')
ind_EC_num = data_MM_pregroup.columns.get_loc('ec_number')

for i in range(len(data_MM_pregroup)):
    # Create a column for Enzyme Class
    temp = data_MM_pregroup.iloc[i,ind_EC_num]
    ind = temp.find('.')
    data_MM_pregroup.iloc[i,ind_EC] = temp[0:ind]

# Calculate the standard deviation of the entire dataset
stdevtot = np.std(data_MM_pregroup['pKm'])

#%%# Plotting pairs of Measured Values
#### Control graph or Get pairs from total dataset ####
df_pairs_tot,Lpairs = get_pairs_tot(data_MM_pregroup)

# Calculate standardized residuals and find outliers
idx_std_res_tot,num_out_tot,rsquared_tot,linreg_tot,model_tot = residuals(df_pairs_tot,True)

# Perform 2-Sample Kolmolgorov-Smirnov Test on pairs taken from entire dataset
# Need 'value 1', 'value 2', 'Diff', and 'Diff^2' columns
MUE_tot,MDUE_tot,stdev_tot,R2p_tot,R2pmax_tot,Dstdev_tot = measuresofquality(stdevtot,df_pairs_tot)

# Pair Plot
title = 'Plot of Pairs of pKm Values Taken Randomly \n from Entire Dataset (only MM/{} Pairs)'.format(Lpairs)
filename = 'figures/Pair_plot_Control.png'
plot_outliers = False
save = True
slope_tot = pairplot(df_pairs_tot,plot_outliers,idx_std_res_tot,R2pmax_tot,title,save,filename,linreg_tot)

#%% Plot a Histogram of the Entire Dataset
d1,freq1 = np.unique(data_MM_pregroup['pKm'],return_counts=True)

plt.figure()
plt.hist(data_MM_pregroup['pKm'], bins = np.arange(min(d1),max(d1),0.5),color= 'white',edgecolor='k')

plt.xlabel('$pK_{M}$')
plt.ylabel('Frequency')
#plt.title('Entire Dataset')
plt.tight_layout()
plt.show()

#%% Calculate mean and confidence intervals of data
mean_tot,lower_tot,upper_tot = mean_confidence_interval(MUE_tot,stdev_tot,len(df_pairs_tot))

#%%# Create data for pair analysis by retrieving pKm values and
# lowering minimum replicate requirement with the following groups
groups = ['ec_number','kineticlaw_type','substrate','organism','ph','temperature','buffer']

P=3 # Minimum number of values in each group
islogtransf = True # True for pKm and False for Km
cor_units = False # True if inputting data_cor_units (all kinetic law types)
data_all = make_groups(data_MM_pregroup,groups,islogtransf,P,cor_units)
data_all = data_all.reset_index(drop=True)

nw_idx = organize_index(data_all)
data_MM = data_all.reindex(nw_idx)
data_MM = data_MM.reset_index(drop=True)
L2 = len(data_MM)

# Number of measurements after removing negative measurements
print('Number of groups = {}, Number of measurements = {}'.format(len(data_MM),sum(data_MM['count'])))

#%%# Get pairs regardless of number of publications ####
ind_pkm = data_MM.columns.get_loc('Km')
ind_EC = data_MM.columns.get_loc('ec_number')
df_pairs, L_pairs = randomvals_and_diff(data_MM,ind_pkm,ind_EC)

#%% Calculate the measures of quality
# Need 'value 1', 'value 2', 'Diff', and 'Diff^2' columns
MUE,MDUE,stdev,R2p,R2pmax,Dstdev = measuresofquality(stdevtot,df_pairs)

#%% Calculate Standardized residuals and find outliers
idx_std_res,num_out,rsquared,linreg,model = residuals(df_pairs,False)

# Calculate confidence intervals of the difference
mean,lower,upper = mean_confidence_interval(MUE,stdev,len(df_pairs))

# Pair Plot
title = 'Sets with Min {} Replicates \nConsidering All MRCs ({} sets)'.format(P,len(data_MM))
filename = 'figures/Type 2 MRCs - Pair Plot - No publication control.png'
plot_outliers = False
save = True
slope = pairplot(df_pairs,plot_outliers,idx_std_res,R2pmax,title,save,filename,linreg)

# Enzyme Class Analysis
multpairs = False # df_pairs hasn't been split by pubs so only has 1 pair per set
isgrouped = False # For data not grouped by MRC's yet (only for data_MM_pregroup)
df_pairs, unique_EC = EnzymeClassAnalysis(df_pairs,multpairs,isgrouped)

#%% Create data for calculating how many 'independent' measurements there are
### This will be equivalent to the number of pubs per set
P2=0 # Minimum number of values in each group
islogtransf = True # True for pKm and False for Km
cor_units = False # True if inputting data_cor_units (all kinetic law types)
data_all2 = make_groups(data_MM_pregroup,groups,islogtransf,P2,cor_units)
data_all2 = data_all2.reset_index(drop=True)

nw_idx2 = organize_index(data_all2)
data_MM2 = data_all2.reindex(nw_idx2)
data_MM2 = data_MM2.reset_index(drop=True)
data_MM2 = findpubs(data_MM2)

print('Number of measurements = {}'.format(sum(data_MM2['count'])))

#%% Find out how many unique publications are in each set
data_MM = findpubs(data_MM)

#%% Histogram of how many values are in each set
l = np.unique(data_MM['count'],return_counts=True)

print('Proportion of Groups with 3, 4, and 5 entries =',l[1][0:3]/sum(l[1]))

#%% Histogram of number of publications in each set
p = np.unique(data_MM['Number of Pubs'],return_counts=True)

print('Proportion of Groups with 1, 2, and 3 publications =',p[1][0:3]/sum(p[1]))

#%% Combined histogram of how many publications are in each set
bins = np.arange(0,40,1)
plt.hist(data_MM['count'],bins = bins, color= 'white',edgecolor='k', label = 'Entries per set')
plt.hist(data_MM['Number of Pubs'],bins = bins, color= 'white',edgecolor='dodgerblue',label = 'Publications per set')
plt.ylabel('# of sets')
#plt.title('Type 2 MRC')
plt.xlim(0,15)
plt.legend()
plt.tight_layout()
plt.show()

#%% Sort data by first authors
data_MM,df_fixes = findfirstauthors(data_MM)
df_fixes = df_fixes.dropna()

diff_fauthor = []
ind_pubnum = data_MM.columns.get_loc('Number of Pubs')
ind_fauthnum = data_MM.columns.get_loc('Number of First Authors')

# Create df's to compare authors vs pubs, first authors vs pubs, & authors vs first authors
for i in range(len(data_MM)):
    if data_MM.iloc[i,ind_fauthnum] < data_MM.iloc[i,ind_pubnum]:
        diff_fauthor.append(i)

# How many groups have different number of unique authors vs. unique pubs?
print(f'There are {len(diff_fauthor)} groups where a first author contributes multiple publications\n')

# Saves all groups with multiple pubs written by first authors
data_fauthpub = data_MM[data_MM['Number of First Authors']<data_MM['Number of Pubs']]

#%%## Get sets with >1 publication and sort the data by publications/authors ####
bools = [False,True] # Switch between sorting by publication and author
counts = [0]*2
for i in range(len(bools)):
    df_multpubs,df_pairs_2 = sortbypubs(data_MM,bools[i],data_fauthpub)
    Lmp = len(df_multpubs)

    # Measures of Quality of Data from the Difference between pairs
    df_allpairs = getallpairs(df_pairs_2)

    # Calculate the measures of quality based on differences
    # Need 'value 1', 'value 2', 'Diff', and 'Diff^2' columns
    MUE_all,MDUE_all,stdev_all,R2p_all,R2pmax_all,Dstdev_all = measuresofquality(stdevtot,df_allpairs)

    # Calculate standardized residuals and find outliers
    idx_std_res_all,num_out_all,rsquared_all,linreg_all,model_all = residuals(df_allpairs,False)

    # Pair Plot
    title = 'All Pairs in Sets with >1 Publication, and \n Considering All MRCs ({} Pairs)'.format(len(df_allpairs))
    filename = 'figures/Pair_plot_3vals_all_MRCs_>1pubs_Allpairs.pdf'
    plot_outliers = False
    save = False
    slope_all = pairplot(df_allpairs,plot_outliers,idx_std_res_all,R2pmax_all,title,save,filename,linreg_all)
    counts[i] = len(df_allpairs)

print('\nThis difference should be the same as # of First authors w >1 pub -->',counts[0]-counts[1])

#%% Make Diff Ratio Column
df_pairs_tot['Diff Ratio'] = [10 ** num for num in abs(df_pairs_tot['Diff'])]
df_pairs['Diff Ratio'] = [10 ** num for num in abs(df_pairs['Diff'])]
df_allpairs['Diff Ratio'] = [10 ** num for num in abs(df_allpairs['Diff'])]

#%% Create plot to illustrate Mean and Confidence Intervals
# Create DataFrame to save Mean-CI data
cols = ['Lower CI','Upper CI']
df_mean_CI = pd.DataFrame({'Entire Dataset':[lower_tot,upper_tot]},
                          index = cols)
# Calculate confidence intervals of the difference
mean_all,lower_all,upper_all = mean_confidence_interval(MUE_all,stdev_all,len(df_allpairs))

# Add data for datasets above
df_mean_CI['No Pub Control'] = [lower,upper]
#df_mean_CI['Pub Control'] = [m_2,lower_2,upper_2] # Should always be plotting all pairs
df_mean_CI['Pub Control (All Pairs)'] = [lower_all,upper_all]

# Create columns for ratios
df_mean_CI = df_mean_CI.transpose()
mean = [MUE_tot,MUE,MUE_all]
mean_r = []
lower_CI_r = []
upper_CI_r = []

for i in range(len(df_mean_CI)):
    mean_r.append(10 ** mean[i])
    lower_CI_r.append(10 ** df_mean_CI.iloc[i,0])
    upper_CI_r.append(10 ** df_mean_CI.iloc[i,1])

df_mean_CI['MUE Ratio'] = mean_r
df_mean_CI['Lower CI Ratio'] = lower_CI_r
df_mean_CI['Upper CI Ratio'] = upper_CI_r

# Make Box Plot
# fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(14,10))
# fig.suptitle('$pK_{M}$ Differences vs Km Ratios', fontsize=22, y=0.99)

# Creat list of the data being plotted
dats=[df_pairs_tot['Diff'],df_pairs['Diff'],df_allpairs['Diff'],
      df_pairs_tot['Diff Ratio'],df_pairs['Diff Ratio'],df_allpairs['Diff Ratio']]

# Calculate the percent of outliers
list_out = [0]*len(dats)

for i in range(len(dats)):
    Q1 = dats[i].quantile(0.25)
    Q3 = dats[i].quantile(0.75)
    IQR = Q3 - Q1

    list_out[i] = round(((dats[i] < (Q1 - 1.5 * IQR)) | (dats[i] > (Q3 + 1.5 * IQR))).sum()*100/len(dats[i]),1)

# Plot the data
# b1 = ax1.boxplot(dats[0],showfliers=False,labels = [f'% Outliers = {list_out[0]}'])
# b2 = ax2.boxplot(dats[1],showfliers=False,labels = [f'% Outliers = {list_out[1]}'])
# b3 = ax3.boxplot(dats[2],showfliers=False,labels = [f'% Outliers = {list_out[2]}'])
# b4 = ax4.boxplot(dats[3],showfliers=False,labels = [f'% Outliers = {list_out[3]}'])
# b5 = ax5.boxplot(dats[4],showfliers=False,labels = [f'% Outliers = {list_out[4]}'])
# b6 = ax6.boxplot(dats[5],showfliers=False,labels = [f'% Outliers = {list_out[5]}'])

# Median and max of the different pair sets
mmp = [[round(max(dats[0]),2),round(np.median(dats[0]),2)],
       [round(max(dats[1]),2),round(np.median(dats[1]),2)],
       [round(max(dats[2]),2),round(np.median(dats[2]),2)],
       [round(max(dats[3]),2),round(np.median(dats[3]),2)],
       [round(max(dats[4]),2),round(np.median(dats[4]),2)],
       [round(max(dats[5]),2),round(np.median(dats[5]),2)]]

# bs = [b1,b2,b3,b4,b5,b6]
# axs = [ax1,ax2,ax3,ax4,ax5,ax6]
# [x.set_xticks([1]) for x in axs]
# for i in range(len(axs)):
#     axs[i].legend([bs[i]["boxes"][0]], [f'% Outliers = {list_out[i]}%'], loc='upper right',handlelength=0)

# ax1.set_xticklabels(['Entire Dataset \n (n = {}) \n (Max = {}, Med = {})'.format(len(df_pairs_tot),mmp[0][0],mmp[0][1])],fontsize=12)
# ax2.set_xticklabels(['No Pub Control \n (n = {}) \n (Max = {}, Med = {})'.format(len(df_pairs),mmp[1][0],mmp[1][1])],fontsize=12)
# ax3.set_xticklabels(['Pub Control \n (n = {}) \n (Max = {}, Med = {})'.format(len(df_allpairs),mmp[2][0],mmp[2][1])],fontsize=12)
# ax4.set_xticklabels(['Entire Dataset \n (n = {}) \n (Max = {}, Med = {})'.format(len(df_pairs_tot),mmp[3][0],mmp[3][1])],fontsize=12)
# ax5.set_xticklabels(['No Pub Control \n (n = {}) \n (Max = {}, Med = {})'.format(len(df_pairs),mmp[4][0],mmp[4][1])],fontsize=12)
# ax6.set_xticklabels(['Pub Control \n (n = {}) \n (Max = {}, Med = {})'.format(len(df_allpairs),mmp[5][0],mmp[5][1])],fontsize=12)

# ax1.set_ylim(-0.1,5)
# ax2.set_ylim(-0.1,1.4)
# ax4.set_ylim(-0.1,350)
# ax3.set_ylim(-0.1,1)
# ax5.set_ylim(-0.1,10)
# ax6.set_ylim(-0.1,6)

# ax1.set_ylabel('|$pK_{M,1}$ - $pK_{M,2}$|', fontsize=16,x=0.5)
# ax4.set_ylabel('$K_{M,1}$/$K_{M,2}$', fontsize=16,x=0.5)

# # ax2.set_title('pKm Differences',fontsize=16)
# # ax5.set_title('Km Ratios',fontsize=16)

# plt.tight_layout()
# plt.subplots_adjust(hspace=0.3)
# plt.show()

#%% MAKE A DATAFRAME OF ALL IMPORTANT STATISTICAL VALUES FOR THE THREE IMPORTANT GROUPS OF DATA
df_statistics = pd.DataFrame({'Count': [len(df_pairs_tot),len(df_pairs),len(df_allpairs)],
    'MUE': [MUE_tot,MUE,MUE_all],
    'MDUE': [MDUE_tot,MDUE,MDUE_all],
    'R2 Pearson': [R2p_tot,R2p,R2p_all],
    'R2 Pearson Max': [R2pmax_tot,R2pmax,R2pmax_all],
    'Stdev': [stdev_tot,stdev,stdev_all],
    'Dahlberg Stdev':[Dstdev_tot,Dstdev,Dstdev_all],
    'Num Outliers':[num_out_tot,num_out,num_out_all],
    'LR Rsquared':[rsquared_tot,rsquared,rsquared_all]},index=df_mean_CI.index)

df_mean_CI['Med Diff'] = [mmp[0][1],mmp[1][1],mmp[2][1]]
df_mean_CI['Max Diff'] = [mmp[0][0],mmp[1][0],mmp[2][0]]
df_mean_CI['Med Ratio'] = [mmp[3][1],mmp[4][1],mmp[5][1]]
df_mean_CI['Max Ratio'] = [mmp[3][0],mmp[4][0],mmp[5][0]]
df_mean_CI['% Outliers Diff'] = list_out[0:3]
df_mean_CI['% Outliers Ratio'] = list_out[3:6]
print(df_mean_CI)

df_statistics = df_statistics.join(df_mean_CI)
df_statistics.to_csv('Stats_tot_MRC2.csv')

#%% Plot the distribution of the mean ratio (Using df_allpairs)
plt.figure()
bp = plt.boxplot(df_allpairs['Diff Ratio'])
plt.ylabel('$K_{M}$ Ratio')
plt.xticks([1],labels=[f'Km Ratio of All Pairs (n = {len(df_allpairs)})'])
plt.ylim(0,10)

for line in bp['medians']:
    # get position data for median line
    x, y = line.get_xydata()[1]
    # overlay median value
    plt.text(x, y, '%.1f' % y, verticalalignment='center') # draw above, centered

for line in bp['boxes']:
    x, y = line.get_xydata()[0]
    plt.text(x,y, '%.1f' % y,
         verticalalignment='top',
         horizontalalignment='right')
    x, y = line.get_xydata()[3]
    plt.text(x,y, '%.1f' % y,
         verticalalignment='bottom',
         horizontalalignment='right')

plt.legend([f'Max = {round(max(df_allpairs["Diff Ratio"]),2)}'],loc='upper right')
plt.tight_layout()
plt.show()

#%% Perform Enzyme Class Analysis
#multpairs = True # df_multpubs has multiple pairs after being controlled for pubs
#isgrouped = False # For data not grouped by MRC's yet (only for data_MM_pregroup)
#df_multpubs, unique_EC_all = EnzymeClassAnalysis(df_multpubs,multpairs,isgrouped)

#%%#
""" CREATE NEW DATASET WITH LESS MRCs """

groups = ['ec_number','kineticlaw_type','substrate','organism']

P=3
islogtransf = True # True for pKm and False for Km
cor_units = False # True if inputting data_cor_units (all kinetic law types)
data_MRC1 = make_groups(data_MM_pregroup,groups,islogtransf,P,cor_units)
data_MRC1 = data_MRC1.reset_index(drop=True)

nw_idx = organize_index(data_MRC1)
data_MM_MRC1 = data_MRC1.reindex(nw_idx)
data_MM_MRC1 = data_MM_MRC1.reset_index(drop=True)
L2 = len(data_MM_MRC1)

print('Number of measurements = {}'.format(sum(data_MM_MRC1['count'])))

#%% Histogram of how many values are in each set
l_2 = np.unique(data_MM_MRC1['count'],return_counts=True)
plt.hist(data_MM_MRC1['count'],bins = np.arange(min(l_2[0]),max(l_2[0]),1), color= 'k')
plt.xlabel('# of entries per set')
plt.ylabel('# of sets')
#plt.title('Type 1 MRC')
plt.ylim(-20,600)
plt.xlim(2,35)
plt.tight_layout()
plt.show()

print('Proportion of Groups with 3, 4, and 5 entries =',l_2[1][0:3]/sum(l_2[1]))

#%% Histogram of how many publications are in each set
data_MM_MRC1 = findpubs(data_MM_MRC1)
p_2 = np.unique(data_MM_MRC1['Number of Pubs'],return_counts=True)

plt.hist(data_MM_MRC1['Number of Pubs'],np.arange(min(p_2[0]),max(p_2[0]),1), color= 'k')
plt.xlabel('# of publications per set')
plt.ylabel('# of sets')
#plt.title('Type 1 MRC')
plt.ylim(-20,750)
plt.xlim(0.5,10)
plt.tight_layout()
plt.show()

print('Proportion of Groups with 1, 2, and 3 publications =',p_2[1][0:3]/sum(p_2[1]))

#%% Combined histogram of how many publications are in each set
bins = np.arange(0,40,1)
plt.hist(data_MM_MRC1['count'],bins = bins, color= 'k', label = 'Entries per set')
plt.hist(data_MM_MRC1['Number of Pubs'],bins = bins,label = 'Publications per set')
plt.ylabel('# of sets')
#plt.title('Type 2 MRC')
plt.xlim(0,15)
plt.legend()
plt.tight_layout()
plt.show()

#%% Find out outliers, measures of quality, and mean±CI for the lower MRC dataset
# before controlling for publications

# Get pairs
ind_pkm = data_MM_MRC1.columns.get_loc('Km')
ind_EC = data_MM_MRC1.columns.get_loc('ec_number')
df_pairs_MRC1, L_pairs_MRC1 = randomvals_and_diff(data_MM_MRC1,ind_pkm,ind_EC)

# Get outliers
idx_std_res_MRC1,num_out_MRC1,rsquared_MRC1,linreg_MRC1,model_MRC1 = residuals(df_pairs_MRC1,False)

# Get stats
MUE_MRC1,MDUE_MRC1,stdev_MRC1,R2p_MRC1,R2pmax_MRC1,Dstdev_MRC1 = measuresofquality(stdevtot,df_pairs_MRC1)

# Get mean and confidence intervals
mean_MRC1,lower_MRC1,upper_MRC1 = mean_confidence_interval(MUE_MRC1,stdev_MRC1,len(df_pairs_MRC1))

# Perform Enzyme Class analysis
multpairs = False
isgrouped = False # For data not grouped by MRC's yet (only for data_MM_pregroup)
df_pairs_MRC1, unique_EC_MRC1 = EnzymeClassAnalysis(df_pairs_MRC1,multpairs,isgrouped)

# Pair Plot
title = 'Sets with Min {} Replicates Considering Enzyme, \n Substrate, and Organism ({} pairs)'.format(P,L_pairs_MRC1)
filename = 'figures/Type 1 MRCs - Pair Plot - No publication control.png'
plot_outliers = False
save = True
slope_MRC1 = pairplot(df_pairs_MRC1,plot_outliers,idx_std_res_MRC1,R2pmax_MRC1,title,save,filename,linreg_MRC1)

#%% Do same analysis after sorting dataset by publications
data_MM_MRC1,df_fixes_MRC1 = findfirstauthors(data_MM_MRC1)

# Saves all groups with multiple pubs written by first authors
data_fauthpub_MRC1 = data_MM_MRC1[data_MM_MRC1['Number of First Authors']<data_MM_MRC1['Number of Pubs']]

bools = [False,True]
counts = [0]*2
for i in range(len(bools)):
    df_multpubs_MRC1,df_pairs_MRC1_2 = sortbypubs(data_MM_MRC1,bools[i],data_fauthpub_MRC1)
    Lmp = len(df_multpubs_MRC1)

    # Get all pairs
    df_allpairs_MRC1 = getallpairs(df_pairs_MRC1_2)

    # Get stats
    MUE_MRC1_pub,MDUE_MRC1_pub,stdev_MRC1_pub,R2p_MRC1_pub,R2pmax_MRC1_pub,Dstdev_MRC1_pub = measuresofquality(stdevtot,df_allpairs_MRC1)

    # Calculate standardized residuals and find outliers
    idx_std_res_MRC1_pub,num_out_MRC1_pub,rsquared_MRC1_pub,linreg_MRC1_pub,model_MRC1_pub = residuals(df_allpairs_MRC1,False)

    # Calculate confidence intervals of the difference
    mean_MRC1_pub,lower_MRC1_pub,upper_MRC1_pub = mean_confidence_interval(MUE_MRC1_pub,stdev_MRC1_pub,len(df_allpairs_MRC1))

    # Perform Enzyme Class analysis
    # multpairs = True
    # isgrouped = False # For data not grouped by MRC's yet (only for data_MM_pregroup)
    # df_multpubs_MRC1_pub, unique_EC_MRC1_pub = EnzymeClassAnalysis(df_multpubs_MRC1,multpairs,isgrouped)

    # Pair Plot
    title = 'All Pairs in Sets with >1 Publication, and Considering \n Enzyme, Substrate, and Organism ({} of {} sets)'.format(Lmp,L2)
    filename = 'figures/Pair_plot_3vals_Enz_Sub_Org_>1pubs.pdf'
    plot_outliers = False
    save = False
    slope_MRC1_pub = pairplot(df_allpairs_MRC1,plot_outliers,idx_std_res_MRC1_pub,rsquared_MRC1_pub,title,save,filename,linreg_MRC1_pub)

#%% MAKE A DATAFRAME OF ALL IMPORTANT STATISTICAL VALUES FOR THE THREE IMPORTANT GROUPS OF DATA
# Create DataFrame to save Mean-CI data
cols = ['Lower CI','Upper CI']
df_mean_CI_MRC1 = pd.DataFrame({'Entire Dataset':[lower_tot,upper_tot]},
                          index = cols)

# Add data for datasets above
df_mean_CI_MRC1['No Pub Control'] = [lower_MRC1,upper_MRC1]
#df_mean_CI['Pub Control'] = [m_2,lower_2,upper_2] # Should always be plotting all pairs
df_mean_CI_MRC1['Pub Control (All Pairs)'] = [lower_MRC1_pub,upper_MRC1_pub]

# Create columns for ratios
df_mean_CI_MRC1 = df_mean_CI_MRC1.transpose()
mean_MRC1 = [MUE_tot,MUE_MRC1,MUE_MRC1_pub]
mean_r_MRC1 = []
lower_CI_r_MRC1 = []
upper_CI_r_MRC1 = []

for i in range(len(df_mean_CI_MRC1)):
    mean_r_MRC1.append(10 ** mean_MRC1[i])
    lower_CI_r_MRC1.append(10 ** df_mean_CI_MRC1.iloc[i,0])
    upper_CI_r_MRC1.append(10 ** df_mean_CI_MRC1.iloc[i,1])

df_mean_CI_MRC1['MUE Ratio'] = mean_r_MRC1
df_mean_CI_MRC1['Lower CI Ratio'] = lower_CI_r_MRC1
df_mean_CI_MRC1['Upper CI Ratio'] = upper_CI_r_MRC1

df_statistics_MRC1 = pd.DataFrame({'MUE': [MUE_tot,MUE_MRC1,MUE_MRC1_pub],
                    'MDUE': [MDUE_tot,MDUE_MRC1,MDUE_MRC1_pub],
                    'R2 Pearson': [R2p_tot,R2p_MRC1,R2p_MRC1_pub],
                    'R2 Pearson Max': [R2pmax_tot,R2pmax_MRC1,R2pmax_MRC1_pub],
                    'Stdev': [stdev_tot,stdev_MRC1,stdev_MRC1_pub],
                    'Dahlberg Stdev':[Dstdev_tot,Dstdev_MRC1,Dstdev_MRC1_pub],
                    'Num Outliers':[num_out_tot,num_out_MRC1,num_out_MRC1_pub],
                    'LR Rsquared':[rsquared_tot,rsquared_MRC1,rsquared_MRC1_pub],
                    'LR Slope':[slope_tot,slope_MRC1,slope_MRC1_pub]},index=df_mean_CI_MRC1.index)

# Creat list of the data being plotted
df_pairs_MRC1['Diff Ratio'] = [10 ** num for num in abs(df_pairs_MRC1['Diff'])]
df_allpairs_MRC1['Diff Ratio'] = [10 ** num for num in abs(df_allpairs_MRC1['Diff'])]

dats_MRC1=[df_pairs_tot['Diff'],df_pairs_MRC1['Diff'],df_allpairs_MRC1['Diff'],
      df_pairs_tot['Diff Ratio'],df_pairs_MRC1['Diff Ratio'],df_allpairs_MRC1['Diff Ratio']]

# Calculate the percent of outliers
list_out_MRC1 = [0]*len(dats_MRC1)

for i in range(len(dats)):
    Q1 = dats[i].quantile(0.25)
    Q3 = dats[i].quantile(0.75)
    IQR = Q3 - Q1

    list_out_MRC1[i] = round(((dats_MRC1[i] < (Q1 - 1.5 * IQR)) | (dats_MRC1[i] > (Q3 + 1.5 * IQR))).sum()*100/len(dats_MRC1[i]),1)

# Median and max of the different pair sets
mmp_MRC1 = [[round(max(dats_MRC1[0]),2),round(np.median(dats_MRC1[0]),2)],
       [round(max(dats_MRC1[1]),2),round(np.median(dats_MRC1[1]),2)],
       [round(max(dats_MRC1[2]),2),round(np.median(dats_MRC1[2]),2)],
       [round(max(dats_MRC1[3]),2),round(np.median(dats_MRC1[3]),2)],
       [round(max(dats_MRC1[4]),2),round(np.median(dats_MRC1[4]),2)],
       [round(max(dats_MRC1[5]),2),round(np.median(dats_MRC1[5]),2)]]
df_mean_CI_MRC1['Med Diff'] = [mmp_MRC1[0][1],mmp_MRC1[1][1],mmp_MRC1[2][1]]
df_mean_CI_MRC1['Max Diff'] = [mmp_MRC1[0][0],mmp_MRC1[1][0],mmp_MRC1[2][0]]
df_mean_CI_MRC1['Med Ratio'] = [mmp_MRC1[3][1],mmp_MRC1[4][1],mmp_MRC1[5][1]]
df_mean_CI_MRC1['Max Ratio'] = [mmp_MRC1[3][0],mmp_MRC1[4][0],mmp_MRC1[5][0]]
df_mean_CI_MRC1['% Outliers Diff'] = list_out_MRC1[0:3]
df_mean_CI_MRC1['% Outliers Ratio'] = list_out_MRC1[3:6]
print(df_mean_CI_MRC1)

df_statistics_MRC1 = df_statistics_MRC1.join(df_mean_CI_MRC1)
df_statistics_MRC1.to_csv('Stats_MRC1.csv')

#%% Discriptive Stats
set_num = [len(data_MM),len(data_MM_MRC1)]
num_vals = [sum(data_MM['count']),sum(data_MM_MRC1['count'])]
set_num_pub = [len(df_multpubs),len(df_multpubs_MRC1)]

print('Number of sets in type 1 and type 2 MRC, respectively',set_num)
print('Number of sets  after controlling for publications in type 1 and type 2 MRC, respectively',set_num_pub)
print('Total number of entries of all sets in type 1 and type 2 MRC, respectively',num_vals)

# Find unique enzyme classes
multpairs = False
isgrouped = True # For data not grouped by MRC's yet (only for data_MM_pregroup)
df_EC_analysis, unique_EC_tot = EnzymeClassAnalysis(data_MM_pregroup,multpairs,isgrouped)
