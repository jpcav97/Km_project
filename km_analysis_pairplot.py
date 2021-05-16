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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import requests
import math
import statistics as stats
from scipy import stats as st
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from difflib import SequenceMatcher as sq

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

def measuresofquality(data):
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
    
    d_ks2,p_ks2 = st.ks_2samp(data['value 1'], data['value 2'])
    print('\n KS test D statistic = {}'.format(round((d_ks2),2)))
    print('\n D-critical statistic = {}'.format(round((p_ks2),2)))
    
    # Histogram of Difference between pairs
    u = np.unique(df_pairs_tot['Diff'])
    
    plt.hist(data['Diff'], bins=np.arange(min(u),max(u),0.2), color= 'white',edgecolor='k')
    plt.xlabel('Difference')
    plt.ylabel('# of pairs')
    plt.title('Difference Between Pairs')
    plt.tight_layout()
    plt.show()
    if d_ks2 < p_ks2:
        print('\n Measurements come from a normal distribution')
    else:
        print('Measurements do not come from a normal distribution')
    
    return MUE,MUE2,MDUE,MDUE2,stdev,R2p,R2pmax,Dstdev,d_ks2,p_ks2

def randomvals_and_diff(data,ind_pkm,ind_EC):
    L = len(data)
    km_list1 = [0]*L
    km_list2 = [0]*L
    EC_list = [0]*L
    
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

#%% Load in Pre-saved 'data_MM_pregroup' data ###
data_MM_pregroup = pd.read_csv('data_MM_pregroup_pubs_auth.csv',index_col=0)
Ltot = len(data_MM_pregroup)
ind_pkm = data_MM_pregroup.columns.get_loc('pKm')

data_MM_pregroup['Enzyme Class'] = ''
ind_EC = data_MM_pregroup.columns.get_loc('Enzyme Class')
ind_EC_num = data_MM_pregroup.columns.get_loc('ec_number')


for i in range(len(data_MM_pregroup)):
    # Create a column for Enzyme Class
    temp = data_MM_pregroup.iloc[i,ind_EC_num]
    ind = temp.find('.')
    data_MM_pregroup.iloc[i,ind_EC] = temp[0:ind]
    
#%%# Plotting pairs of Measured Values
#### Control graph ####
import random
ind_pkm = data_MM_pregroup.columns.get_loc('pKm')
Lpairs = int(len(data_MM_pregroup)/2)
km_list1 = [0]*Lpairs
km_list2 = [0]*Lpairs

for i in range(Lpairs):
    temp = random.sample(list(data_MM_pregroup.iloc[:,ind_pkm]),2)
    km_list1[i] = temp[0]
    km_list2[i] = temp[1]
    
df_pairs_tot = pd.DataFrame({'value 1': km_list1, 'value 2': km_list2})

# Create column for their difference
df_pairs_tot['Diff'] = 0
for i in range(Lpairs):
    df_pairs_tot.iloc[i,-1] = df_pairs_tot.iloc[i,0] - df_pairs_tot.iloc[i,1]

df_pairs_tot['Diff^2'] = [num ** 2 for num in df_pairs_tot['Diff']]

# Calculate standardized residuals and find outliers
idx_std_res_tot,num_out_tot,rsquared_tot,linreg_tot,model_tot = residuals(df_pairs_tot,True)

# Perform 2-Sample Kolmolgorov-Smirnov Test on pairs taken from entire dataset
# Need 'value 1', 'value 2', 'Diff', and 'Diff^2' columns
MUE_tot,MUE2_tot,MDUE_tot,MDUE2_tot,stdev_tot,R2p_tot,R2pmax_tot,Dstdev_tot,d_ks2_tot,p_ks2_tot = measuresofquality(df_pairs_tot)

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

P=3
data_all = make_groups(data_MM_pregroup,groups,True,P)
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
MUE,MUE2,MDUE,MDUE2,stdev,R2p,R2pmax,Dstdev,d_ks2,p_ks2 = measuresofquality(df_pairs)

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

# Pair Plot considering groups with only one pub and more than one pair if unique pubs > 3
#df_all,temp = sortbypubs(data_MM,False,np.nan,True)
#df_all_pairs = getallpairs(temp)

#idx_std_res,num_out,rsquared,linreg = residuals(df_all_pairs)

# Calculate confidence intervals of the difference
#mean,lower,upper = mean_confidence_interval(df_all_pairs['Diff'])

# title = 'Sets with Min {} Replicates \nConsidering All MRCs ({} sets)'.format(P,len(data_MM))
# filename = 'figures/Pair_plot_3vals_all_MRCs.pdf'
# plot_outliers = False
# save = True
# slope = pairplot(df_all_pairs,plot_outliers,idx_std_res,rsquared,title,save,filename,linreg)

# Enzyme Class Analysis
multpairs = False # df_pairs hasn't been split by pubs so only has 1 pair per set
isgrouped = False # For data not grouped by MRC's yet (only for data_MM_pregroup)
df_pairs, unique_EC = EnzymeClassAnalysis(df_pairs,multpairs,isgrouped)

#%% Create data for calculating how many 'independent' measurements there are
### This will be equivalent to the number of pubs per set
P2=0
data_all2 = make_groups(data_MM_pregroup,groups,True,P2)
data_all2 = data_all2.reset_index(drop=True)

nw_idx2 = organize_index(data_all2)
data_MM2 = data_all2.reindex(nw_idx2)
data_MM2 = data_MM2.reset_index(drop=True)
data_MM2 = findpubs(data_MM2)

print('Number of measurements = {}'.format(sum(data_MM2['count'])))

#%% Find out how many unique pubs are in each set
data_MM = findpubs(data_MM)

#%% Histogram of how many values are in each set
l = np.unique(data_MM['count'],return_counts=True)
plt.hist(data_MM['count'],bins = np.arange(min(l[0]),max(l[0]),1), color= 'white',edgecolor='k')
plt.xlabel('# of entries per set')
plt.ylabel('# of sets')
#plt.title('Type 2 MRC')
#plt.ylim(-20,600)
plt.xlim(2,35)
plt.tight_layout()
plt.show()

print('Proportion of Groups with 3, 4, and 5 entries =',l[1][0:3]/sum(l[1]))
#%% Histogram of how many publications are in each set
p = np.unique(data_MM['Number of Pubs'],return_counts=True)

plt.hist(data_MM['Number of Pubs'],np.arange(min(p[0]),max(p[0]),1), color= 'white',edgecolor='k')
plt.xlabel('# of publications per set')
plt.ylabel('# of sets')
#plt.title('Type 2 MRC')
plt.ylim(-20,800)
plt.xlim(0.5,6)
plt.tight_layout()
plt.show()

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

#%%## Get sets with >1 publication and plot that data #### 
bools = [False,True] 
counts = [0]*2
for i in range(len(bools)):
    df_multpubs,df_pairs_2 = sortbypubs(data_MM,bools[i],data_fauthpub,False)
    Lmp = len(df_multpubs)
    
    # Calculate standardized residuals and find outliers
    idx_std_res_2,num_out_2,rsquared_2,linreg_2,model_2 = residuals(df_pairs_2,False)
    
    # Measures of Quality of Data from the Difference between pairs
    df_allpairs = getallpairs(df_pairs_2)
    
    # Calculate the measures of quality based on differences
    # Need 'value 1', 'value 2', 'Diff', and 'Diff^2' columns
    MUE_all,MUE2_all,MDUE_all,MDUE2_all,stdev_all,R2p_all,R2pmax_all,Dstdev_all,d_ks2_all,p_ks2_all = measuresofquality(df_allpairs)

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
    'KS test D-stat':[d_ks2_tot,d_ks2,d_ks2_all],
    'D-critical stat':[p_ks2_tot,p_ks2,p_ks2_all],
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
""" CREATE NEW DATASET WITH LESS MRC's """

groups = ['ec_number','kineticlaw_type','substrate','organism'] 

P=3
data_MRC1 = make_groups(data_MM_pregroup,groups,True,P)
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
MUE_MRC1,MUE2_MRC1,MDUE_MRC1,MDUE2_MRC1,stdev_MRC1,R2p_MRC1,R2pmax_MRC1,Dstdev_MRC1,d_ks2_MRC1,p_ks2_MRC1 = measuresofquality(df_pairs_MRC1)

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
    df_multpubs_MRC1,df_pairs_MRC1_2 = sortbypubs(data_MM_MRC1,bools[i],data_fauthpub_MRC1,False)
    Lmp = len(df_multpubs_MRC1)
    
    # Get all pairs
    df_allpairs_MRC1 = getallpairs(df_pairs_MRC1_2)
    
    # Get stats
    MUE_MRC1_pub,MUE2_MRC1_pub,MDUE_MRC1_pub,MDUE2_MRC1_pub,stdev_MRC1_pub,R2p_MRC1_pub,R2pmax_MRC1_pub,Dstdev_MRC1_pub,d_ks2_MRC1_pub,p_ks2_MRC1_pub = measuresofquality(df_allpairs_MRC1)
    
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
                              'KS test D-stat':[d_ks2_tot,d_ks2_MRC1,d_ks2_MRC1_pub],
                              'D-critical stat':[p_ks2_tot,p_ks2_MRC1,p_ks2_MRC1_pub],
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
