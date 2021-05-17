#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 16 14:58:13 2021

@author: josephcavataio
"""
# Pairplot analysis functions

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import statistics as stats
from scipy import stats as st
import statsmodels.api as sm
from statsmodels.distributions.empirical_distribution import ECDF
from sklearn.linear_model import LinearRegression
from difflib import SequenceMatcher as sq
import random


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
    
def mean_confidence_interval(m,s,L):
    sem = s/math.sqrt(L) # Figure out the z value CI = mean ± z * mean/√n
    h = sem * 1.96
    return m, m-h, m+h

# nopub == True when you want to get one pair from all sets, and multiple pairs from sets with more than 3 pubs
def sortbypubs(data,sortbyauth,data_fauthpub,nopub):
    L2 = len(data)
    multpubs = []
    num_pairs = []
    count = 0
    
    # Create columns and get indices of Km values for different pubs
    ind_km = data.columns.get_loc('Km')
    values = ['value 1','value 2','value 3','value 4','value 5','value 6']
    ind_kms = [0]*6
    for i in range(len(values)):
        data[values[i]] = np.empty((L2,0)).tolist()
        ind_kms[i] = data.columns.get_loc(values[i])

    # Consider only the sets with pairs taken from different first authors
    if sortbyauth:
        data = data.drop(list(data_fauthpub.index),axis=0)
        data = data.reset_index(drop=True)
        L2 = len(data)
    
    # Place random values from diff pubs into separate columns
    ind = data.columns.get_loc('Publication Title')   
    for i in range(L2):
        # Find unique pubs for each set, their counts, and indices
        temp_pub = data.iloc[i,ind]
        pub_uniq,pub_ind,pub_freq = np.unique(temp_pub,return_index = True,return_counts = True, )
        uniq_pubs = pd.DataFrame(list(zip(pub_uniq,pub_freq,pub_ind)),columns=['names','frequency','indices'])
        
        # Sort unique pubs starting with publication with index 0 in the set
        uniq_pubs = uniq_pubs.sort_values(by=['indices'], ascending = True)
        uniq_pubs = uniq_pubs.reset_index(drop = True)
        
        pub_num = len(uniq_pubs)
        count_pairs = 0
    
        # Find number of pairs based on number of unique pubs
        if pub_num >= 1 and pub_num < 4 and nopub == True:
            count_pairs = 1
        elif pub_num >= 2 and pub_num < 4:
            count_pairs = 1
        elif pub_num >= 4 and pub_num < 6:
            count_pairs = 2
        elif pub_num >= 6:
            count_pairs = 3
        
        # Record the number of pairs in each set
        num_pairs.append(count_pairs)
        
        # If a set contains more than one unique publication
        if pub_num == 1 and nopub == True:
            multpubs.append(i)
            temp_km = data.iloc[i,ind_km] # Gather the sets km values

            # Cycle through unique pubs
            for j in range(pub_num):
                    if j < len(ind_kms):
                        # Divide km values based on publication
                        start = uniq_pubs.iloc[j,2]
                        end = start + uniq_pubs.iloc[j,1]
                        data.iloc[i,ind_kms[j]].extend(temp_km[start:end])
                        
        elif pub_num > 1: 
            count = count + 1
            multpubs.append(i) # Save index
            temp_km = data.iloc[i,ind_km] # Gather the sets km values
            
            # Cycle through unique pubs
            for j in range(pub_num):
                    if j < len(ind_kms):
                        # Divide km values based on publication
                        start = uniq_pubs.iloc[j,2]
                        end = start + uniq_pubs.iloc[j,1]
                        data.iloc[i,ind_kms[j]].extend(temp_km[start:end])
                    
                
    print('\n{}% of sets in this analysis have entries from the same paper\n'.format(round((L2-count)/L2,2)*100))
    
    # Add column for number of pairs
    data['num pairs'] = pd.Series(num_pairs)
    
    # Create dataframe with only sets with multiple publications
    if nopub == True:
        df_multpubs = data
    else:
        df_multpubs = pd.DataFrame(columns = data.columns, index = multpubs) 
        for i in range(len(multpubs)):
            df_multpubs.iloc[i,:] = data.iloc[multpubs[i],:]
    
    Lmp = len(df_multpubs)
    ind_numpair = data.columns.get_loc('num pairs')
    km_list1 = [0]*Lmp
    km_list2 = [0]*Lmp
    km_list3 = [0]*Lmp
    km_list4 = [0]*Lmp
    km_list5 = [0]*Lmp
    km_list6 = [0]*Lmp
    
    # Take random values from different publications and create pairs
    for i in range(Lmp):
        rand_km1 = random.choice(df_multpubs.iloc[i,ind_kms[0]])

        if data.iloc[i,ind_numpair] == 1 and nopub == True:
            rand_km2 = random.choice(df_multpubs.iloc[i,ind_kms[0]])
            test = 0
            while rand_km1 == rand_km2:
                if test == 15:
                    break
                test = test + 1
                
                rand_km2 = random.choice(df_multpubs.iloc[i,ind_kms[0]])
        elif nopub == False:
            rand_km2 = random.choice(df_multpubs.iloc[i,ind_kms[1]])

        km_list1[i] = rand_km1
        km_list2[i] = rand_km2
        
        if df_multpubs.iloc[i,ind_numpair] == 2:
            rand_km3 = random.choice(df_multpubs.iloc[i,ind_kms[2]])
            rand_km4 = random.choice(df_multpubs.iloc[i,ind_kms[3]])
            km_list3[i] = rand_km3
            km_list4[i] = rand_km4
            
        elif df_multpubs.iloc[i,ind_numpair] == 3:
            rand_km3 = random.choice(df_multpubs.iloc[i,ind_kms[2]])
            rand_km4 = random.choice(df_multpubs.iloc[i,ind_kms[3]])
            rand_km5 = random.choice(df_multpubs.iloc[i,ind_kms[4]])
            rand_km6 = random.choice(df_multpubs.iloc[i,ind_kms[5]])
            km_list3[i] = rand_km3
            km_list4[i] = rand_km4
            km_list5[i] = rand_km5
            km_list6[i] = rand_km6
    
        
    pairs = pd.DataFrame({'value 1': km_list1, 'value 2': km_list2,
                               'value 3': km_list3, 'value 4': km_list4,
                               'value 5': km_list5, 'value 6': km_list6,
                               'num pairs': df_multpubs['num pairs']})
    
    # Create column for their difference
    pairs['Diff'] = 0
    for i in range(len(pairs)):
        pairs.iloc[i,-1] = pairs.iloc[i,0] - pairs.iloc[i,1]
        
    return df_multpubs, pairs

def getallpairs(dat):
    L_pairs = sum(dat['num pairs'])
    ind_numpair = dat.columns.get_loc('num pairs')
    ind_val1 = dat.columns.get_loc('value 1')
    ind_val2 = dat.columns.get_loc('value 2')
    ind_val3 = dat.columns.get_loc('value 3')
    ind_val4 = dat.columns.get_loc('value 4')
    ind_val5 = dat.columns.get_loc('value 5')
    ind_val6 = dat.columns.get_loc('value 6')

    diff_list = [0]*L_pairs
    pair_list_full1 = [0]*L_pairs
    pair_list_full2 = [0]*L_pairs
    i = 0
    
    # Find difference in pairs depending on number of pairs available
    # Create a full list of all pairs considering multiple pairs from sets that 
    # contain more than 2 publications
    
    while i < L_pairs-1:
        for j in range(len(dat)):
            k = dat.iloc[j,ind_numpair]
            if k == 1:
                diff_list[i] = dat.iloc[j,ind_val1] - dat.iloc[j,ind_val2]
                
                pair_list_full1[i] = dat.iloc[j,ind_val1]
                pair_list_full2[i] = dat.iloc[j,ind_val2]
    
            elif k == 2:
                diff_list[i] = dat.iloc[j,ind_val1] - dat.iloc[j,ind_val2]
                diff_list[i+1] = dat.iloc[j,ind_val3] - dat.iloc[j,ind_val4]
                
                pair_list_full1[i] = dat.iloc[j,ind_val1]
                pair_list_full2[i] = dat.iloc[j,ind_val2]
                pair_list_full1[i+1] = dat.iloc[j,ind_val3]
                pair_list_full2[i+1] = dat.iloc[j,ind_val4]
                
            elif k == 3:
                diff_list[i] = dat.iloc[j,ind_val1] - dat.iloc[j,ind_val2]
                diff_list[i+1] = dat.iloc[j,ind_val3] - dat.iloc[j,ind_val4]
                diff_list[i+2] = dat.iloc[j,ind_val5] - dat.iloc[j,ind_val6]
                
                pair_list_full1[i] = dat.iloc[j,ind_val1]
                pair_list_full2[i] = dat.iloc[j,ind_val2]
                pair_list_full1[i+1] = dat.iloc[j,ind_val3]
                pair_list_full2[i+1] = dat.iloc[j,ind_val4]
                pair_list_full1[i+2] = dat.iloc[j,ind_val5]
                pair_list_full2[i+2] = dat.iloc[j,ind_val6]
            
            i = i + k

    # Create figure of all pairs of values (not just one pair from every set)
    df_allpairs = pd.DataFrame({'value 1': pair_list_full1,'value 2': pair_list_full2,
                                'Diff': diff_list, 'Diff^2': [num ** 2 for num in diff_list]})
    return df_allpairs

def residuals(data,intercept):
    # Calculate standardized residuals and find outliers
    val1 = np.array(data['value 1']).reshape((-1, 1))
    val2 = np.array(data['value 2'])
    
    # Fit linear regression using scikit linear_model
    linreg = LinearRegression(fit_intercept=intercept)
    linreg = linreg.fit(val1, val2)
    
    rsquared = linreg.score(val1,val2)
    
    # Fit linear regression model using statsmodels
    val1 = sm.add_constant(val1)
    
    model = sm.OLS(val2, val1).fit()

    ''' Uncomment below to plot using statsmodels OLS '''
    # fig, ax = plt.subplots(figsize=(8, 4))
    # ax.scatter(val1,val2,color='k',alpha=0.3,edgecolor=None)
    # x_pred = np.linspace(val1.min(), val1.max(), 50)
    # x_pred1 = sm.add_constant(x_pred)
    # y_pred = model.predict(x_pred1)
    # ax.plot(x_pred, y_pred, '-', color='green', linewidth=2)
        
    #create instance of influence
    influence = model.get_influence()
    
    #obtain standardized residuals
    standardized_residuals = influence.resid_studentized_internal
    
    #display standardized residuals
    idx_std_res = np.where(standardized_residuals > 3)[0]
    idx_std_res = np.append(idx_std_res,np.where(standardized_residuals < -3)[0])
    print('The number of outliers = {}'.format(len(idx_std_res)))
    
    num_out = len(idx_std_res)

    #rsquared = model.rsquared
    
    return idx_std_res,num_out,rsquared,linreg,model

def measuresofquality(data_MM_pregroup,data):
    # Mean unsigned error
    L = len(data)
    MUE = (1/math.sqrt(2))*(sum(abs(data['Diff']))/len(data))
    MUE2 = round(pow(10,MUE),2)
    print("\n The Mean Unsigned Error is {} in pKm units (a factor of {} in Km)".format(round(MUE,3),MUE2))
    
    # Median unsigned error
    MDUE = (1/(math.sqrt(2)))* stats.median(abs(data['Diff']))
    MDUE2 = round(pow(10,MDUE),2)
    print("\n The Median Unsigned Error is {} in pKm units (a factor of {})".format(round(MDUE,3),MDUE2))
    
    stdev = math.sqrt((1/(2*(L-1)))*sum(data['Diff^2']))
    print("\n The Standard Deviation is {} in pKm units".format(round(stdev,3)))
    
    mean1 = stats.mean(data['value 1'])
    mean2 = stats.mean(data['value 2'])
    data['val1-mean1'] = data['value 1'] - mean1
    data['val2-mean2'] = data['value 2'] - mean2
    data['(val1-mean1)^2'] = data['val1-mean1'] ** 2
    data['(val2-mean2)^2'] = data['val2-mean2'] ** 2
    
    R2p = abs(sum(data['val1-mean1'] * data['val2-mean2'])\
        /(math.sqrt(sum(data['(val1-mean1)^2']))*math.sqrt(sum(data['(val2-mean2)^2']))))
    print("\n The R^2 Pearson Number is ", round(R2p,3))
    
    stdevtot = np.std(data_MM_pregroup['pKm'])
    R2pmax = 1 - (stdev/stdevtot)**2
    print("\n The Max R^2 Pearson Number is ", round(R2pmax,3))
    
    Dstdev = math.sqrt(sum(data['Diff^2'])/(2*len(data)))
    print('\n The Dahlberg Standard Deviation is {}'.format(round(Dstdev,3)))
    
    d_ks,p_ks = st.ks_2samp(data['value 1'], data['value 2'])
    print('\n 2 Sample KS test D statistic = {}'.format(round((d_ks),2)))
    print('\n 2 Sample D-critical statistic = {}'.format(round((p_ks),2)))
    
    if d_ks < p_ks:
        print('\n Value to Value pairs come from a normal distribution')
    else:
        print('Value to Value pairs do not come from a normal distribution')
    
    # DIFFERENCE BETWEEN PAIRS
    u = np.unique(data['Diff'])
    
    plt.hist(data['Diff'], bins=np.arange(min(u),max(u),0.2), color= 'white',edgecolor='k')
    plt.xlabel('Difference')
    plt.ylabel('# of pairs')
    plt.title('Difference Between Pairs')
    plt.tight_layout()
    plt.show()
    
    # ecdf = ECDF(data['Diff'])
    # plt.plot(ecdf.x, ecdf.y)
    # plt.show()
    
    # sm.qqplot(data['Diff'], line ='45')
    d_ks2,p_ks2 = st.kstest(data['Diff'],'norm')
    if d_ks2 < p_ks2:
        print('\n Differences between pairs come from a normal distribution')
    else:
        print('Differences between pairs do not come from a normal distribution')
    
    return MUE,MUE2,MDUE,MDUE2,stdev,R2p,R2pmax,Dstdev,d_ks,p_ks

def randomvals_and_diff(data,ind_pkm,ind_EC):
    L = len(data)
    km_list1 = [0]*L
    km_list2 = [0]*L
    EC_list = [0]*L

    # ind_count = data.columns.get_loc('count')
    # L2 = sum(data['count'])
    # km_list3 = [0]*L2
    # km_list4 = [0]*L2
    
    # for i in range(L):
    #     x = data.iloc[i,ind_count]
    #     for j in range(x):
    #         for k in range(j+1,x):
    #             km_list3.append(data.iloc[i,ind_pkm][j])
    #             km_list4.append(data.iloc[i,ind_pkm][k])
    
    for i in range(L):
        temp = random.sample(data.iloc[i,ind_pkm],2)
        km_list1[i] = temp[0]
        km_list2[i] = temp[1]
        EC_list[i] = data.iloc[i,ind_EC]
        
    df_pairs = pd.DataFrame({'value 1': km_list1, 'value 2': km_list2, 'ec_number': EC_list})
    
    L_pairs = len(df_pairs)
    
    # Create column for their difference
    df_pairs['Diff'] = 0
    for i in range(len(df_pairs)):
        df_pairs.iloc[i,-1] = df_pairs.iloc[i,0] - df_pairs.iloc[i,1]
    
    df_pairs['Diff^2'] = [num ** 2 for num in df_pairs['Diff']]
    
    return df_pairs, L_pairs

def pairplot(data,plot_outliers,ind_outliers,rsquared,title,save,filename,linreg):
    plt.figure()
    plt.scatter(data['value 1'],data['value 2'],color='k',alpha=0.3,edgecolor=None)
    
    if plot_outliers == True:
        plt.scatter(data.iloc[ind_outliers,0],data.iloc[ind_outliers,1],color='r')
    
    # Fit a line to data with np.polyfit
    x = np.linspace(0,10,100)
    m = linreg.coef_
    b = linreg.intercept_
    plt.plot(x, m*x+b,'g',label = 'R^2 = {}'.format(round(rsquared,2)))
    print('R^2 value = {}'.format(round(rsquared,3)))
    
    plt.xticks(range(0,12,2), fontsize=12)
    plt.yticks(range(0,12,2), fontsize=12)
    plt.ylim(0,10)
    plt.xlim(0,10)
    plt.xlabel('$pK_{M1}$',fontsize=12)
    plt.ylabel('$pK_{M2}$',fontsize=12)
    plt.title('n = {}'.format(len(data)),fontsize = 16, y=1.05)
    plt.legend()
    plt.show()
    
    if save == True:
        plt.savefig(filename,bbox_inches='tight')
    
    return m

def EnzymeClassAnalysis(data,multpairs,isgrouped):
    data['Enzyme Class'] = ''
    ind_EC = data.columns.get_loc('Enzyme Class')
    ind_EC_num = data.columns.get_loc('ec_number')
    
    for i in range(len(data)):
        # Create a column for Enzyme Class
        temp = data.iloc[i,ind_EC_num]
        ind = temp.find('.')
        data.iloc[i,ind_EC] = temp[0:ind]
        
        if multpairs == True:
            # Change pKm values from list to integer
            ind_val1 = data.columns.get_loc('value 1')
            ind_val6 = data.columns.get_loc('value 6')
            
            for j in range(ind_val1,ind_val6+1):
                if len(data.iloc[i,j]) > 0:
                    data.iloc[i,j] = data.iloc[i,j][0]
    
    # Create a boxplot of the different enzyme classes (Using pairs from data)
    EC,freq_EC = np.unique(data['Enzyme Class'],return_counts=True)
    
    unique_EC = pd.DataFrame(list(zip(EC,freq_EC)),columns=['EC','Frequency'])
    a = ['1','2','3','4','5','6']
    unique_EC = unique_EC[unique_EC['EC'].isin(a)] 
        
    unique_EC = unique_EC.reindex(np.argsort(unique_EC['EC']))
    unique_EC  = unique_EC.reset_index(drop=True)
    
    enz_classes = ['Oxidoreductases','Transferases','Hydrolases','Lyases',
                   'Isomerases','Ligases']                       
    unique_EC['EC'] = enz_classes
    
    df_EC1 = data[data['Enzyme Class']=='1']
    df_EC2 = data[data['Enzyme Class']=='2']
    df_EC3 = data[data['Enzyme Class']=='3']
    df_EC4 = data[data['Enzyme Class']=='4']
    df_EC5 = data[data['Enzyme Class']=='5']
    df_EC6 = data[data['Enzyme Class']=='6']
    
    ECs = [df_EC1,df_EC2,df_EC3,df_EC4,df_EC5,df_EC6]
    df_EC_all = [[] for _ in range(len(ECs))]
    
    for i in range(len(ECs)):
        if multpairs == True:
            d = getallpairs(ECs[i])
            d['Ratio'] = [10 ** num for num in abs(d['Diff'])]
            df_EC_all[i] = d['Ratio']
        
        elif isgrouped == True:
            ind_pkm = ECs[i].columns.get_loc('pKm')
            Lpairs = int(len(ECs[i])/2)
            km_list1 = [0]*Lpairs
            km_list2 = [0]*Lpairs
            
            for j in range(Lpairs):
                temp = random.sample(list(ECs[i].iloc[:,ind_pkm]),2)
                km_list1[j] = temp[0]
                km_list2[j] = temp[1]
                
            df_pairs_tot = pd.DataFrame({'value 1': km_list1, 'value 2': km_list2})
            
            # Create column for their difference
            df_pairs_tot['Diff'] = 0
            for k in range(Lpairs):
                df_pairs_tot.iloc[k,-1] = df_pairs_tot.iloc[k,0] - df_pairs_tot.iloc[k,1]
                
            df_pairs_tot['Ratio'] = [10 ** num for num in abs(df_pairs_tot['Diff'])]
            df_EC_all[i] = df_pairs_tot['Ratio']
            
        else:
            temp = ECs[i]['Diff']
            df_EC_all[i] = [10 ** num for num in abs(temp)]
        
    
    # Create a figure instance
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_axes([0,0,1,1])
    
    #Create the boxplot
    bp = ax.boxplot(df_EC_all,showfliers=False)
    ax.set_ylabel('Km Ratio',fontsize=12)
        
    ax.set_xticks([1,2,3,4,5,6]) # This number should match number 2 lines above and always starts from 1
    ax.set_xticklabels(['{} \n (n = {})'.format(unique_EC.iloc[0,0],len(df_EC_all[0])),
                        '{} \n (n = {})'.format(unique_EC.iloc[1,0],len(df_EC_all[1])),
                        '{} \n (n = {})'.format(unique_EC.iloc[2,0],len(df_EC_all[2])),
                        '{} \n (n = {})'.format(unique_EC.iloc[3,0],len(df_EC_all[3])),
                        '{} \n (n = {})'.format(unique_EC.iloc[4,0],len(df_EC_all[4])),
                        '{} \n (n = {})'.format(unique_EC.iloc[5,0],len(df_EC_all[5]))],
                          fontsize=12)
    
    #ax.set_title('Different Enzyme Classes', fontsize=16,x=0.5) # Change titles here
    
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
    
    plt.show()
    
    return data, unique_EC

def findpubs(data):
    L = len(data)
    
    title_lists_pub = [[] for _ in range(L)]
    pub_count = [0]*L
    count = 0
    ind_pub = data.columns.get_loc('Publication Title')
    
    for i in range(L):
        pub_uniq = np.unique(data.iloc[i,ind_pub],return_counts = True)
        title_lists_pub[i] = pub_uniq[0]
        pub_count[i] = pub_uniq[1]
        
        if len(pub_uniq[1]) > 1:
            count = count + 1
    
    data['Publication Title'] = title_lists_pub
    data['Pub Count'] = pub_count
    
    # Find how many publications are in each set
    data['Number of Pubs'] = data['Pub Count'].str.len() 
    
    return data

def findfirstauthors(data_MM):
    #Find out how many unique Authors are in each set
    L2 = len(data_MM)
    title_lists_auth = [[] for _ in range(L2)]
    auth_count = [0]*L2
    
    title_lists_fauth = [[] for _ in range(L2)]
    fauth_count = [0]*L2
    
    ind_auth = data_MM.columns.get_loc('Publication Author')
    ind_fauth = data_MM.columns.get_loc('First Author')
    fixes = [[] for _ in range(L2)]
    
    for i in range(L2):
        auth_uniq = np.unique(data_MM.iloc[i,ind_auth],return_counts = True)
        title_lists_auth[i] = auth_uniq[0]
        auth_count[i] = auth_uniq[1]
    
        fauth_uniq = pd.DataFrame(np.unique(data_MM.iloc[i,ind_fauth],return_counts = True)).transpose()
        for j in range(len(fauth_uniq)):
            for k in range(len(fauth_uniq)):
                if j != k and j < len(fauth_uniq) and k < len(fauth_uniq):
                    one = fauth_uniq.iloc[j,0]
                    two = fauth_uniq.iloc[k,0]
                    sim = sq(None,one,two).ratio()
                    if sim > 0.7:
                        fixes[i] = [j,k,fauth_uniq.iloc[k,1]]
                        fauth_uniq.iloc[j,1] = sum(fauth_uniq.iloc[[j,k],1])
                        fauth_uniq = fauth_uniq.drop(k,axis=0)
                        
        title_lists_fauth[i] = list(fauth_uniq.iloc[:,0])
        fauth_count[i] = list(fauth_uniq.iloc[:,1])
    
    data_MM['Publication Authors'] = title_lists_auth
    data_MM['Author Count'] = auth_count
    data_MM['Number of Authors'] = data_MM['Author Count'].str.len() # Find how many publications are in each set
    
    data_MM['Unique First Authors'] = title_lists_fauth
    data_MM['First Author Count'] = fauth_count
    data_MM['Number of First Authors'] = data_MM['First Author Count'].str.len() # Find how many publications are in each set
    
    df_fixes = pd.DataFrame(fixes,columns = ['Combined','Replaced','Replaced Count'])
    return data_MM, df_fixes

def make_groups(data,MRC,IsLogTransf,P):
    ## Group the entries according to criteria in the groups variable ##
    data0 = data.groupby(MRC)
    group_index = list(data0.indices)
    
    km_vals = []
    km_dev = []
    entry_id = []
    formulas = []
    pub_title = []
    author = []
    fauthor = []

    L = len(group_index)
    print('%d is the no of sets under MRCs.'%L)
    for i in range(L):
        example = data0.get_group(group_index[i])
        entry_id.append(list(example['entry_id']))
        formulas.append(list(example['kineticlaw_formula']))
        pub_title.append(list(example['Publication Title']))
        author.append(list(example['Publication Author']))
        fauthor.append(list(example['First Author']))
        
        
        if IsLogTransf == True:
            km_vals.append(list(example['pKm']))
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
    grp['Publication Title'] = pd.Series(pub_title)
    grp['Publication Author'] = pd.Series(author)
    grp['First Author'] = pd.Series(fauthor)
    
    grp  = grp.reindex(np.argsort(grp['count'])[::-1])
    grp  = grp.reset_index(drop=True)
    dat = grp[grp['count']>=P]
    L2 = len(dat)
    print('{} is the no of sets under MRCs with {} or more entries.\n'.format(L2,P))

    # Record only one formula for each set
    ind_form = dat.columns.get_loc('formula')

    for i in range(L2):
        dat.iloc[i,ind_form] = np.unique(dat.iloc[(i,ind_form)])[0]          
    return dat 
