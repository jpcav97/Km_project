#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 13:49:59 2021

@author: josephcavataio
"""
import numpy as np
import pandas as pd

from functions import organize_index, mean_confidence_interval,sortbypubs, \
    getallpairs, residuals, measuresofquality,randomvals_and_diff, get_pairs_tot, \
        pairplot, EnzymeClassAnalysis, findpubs, findfirstauthors, make_groups

#%% Import raw data (Michaelis-Menten)
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

"""
#%%# Import raw Data non-MM
data = pd.read_csv('km_data_corrected.txt', sep='\t')
data = data.drop(['km_end_value'],axis=1)

#### Clean-up: Sort entries by entry_id ####
data = data.set_index('entry_id', drop=False)
data = data.sort_index(0) #There are duplicate entry IDs for sequential ordered bi-bi entries

#### Clean-up: Units (include only M) ####
unit_uniq = pd.DataFrame(np.unique(data['unit'],return_counts=True)).transpose()
data_cor_units = data[data['unit']==unit_uniq[0][0]]

# Number of measurements after Unit Clean Up
print('Number of measurements after Unit Clean Up = {}'.format(len(data_cor_units)))

data_cor_units=data_cor_units[data_cor_units['km_start_value']>0]
data_cor_units=data_cor_units.reset_index(drop=True)

# Create column with pKm values
data_cor_units['pKm'] = np.nan
ind_km_start = data_cor_units.columns.get_loc('km_start_value')
ind_pkm = data_cor_units.columns.get_loc('pKm')

Ltot = len(data_cor_units)
data_cor_units['pKm'] = [-1*np.log10(x) for x in data_cor_units['km_start_value']]

#%% Create data_MM3
groups = ['ec_number','kineticlaw_type','substrate','organism','ph','temperature']

P=3 # Minimum number of values in each group
islogtransf = True # True for pKm and False for Km
cor_units = False # True if inputting data_cor_units (all kinetic law types)
data_all3 = make_groups(data_MM_pregroup,groups,islogtransf,P,cor_units)
data_all3 = data_all3.reset_index(drop=True)

nw_idx3 = organize_index(data_all3)
data_MM3 = data_all3.reindex(nw_idx3)
data_MM3 = data_MM3.reset_index(drop=True)

print('Number of measurements = {}'.format(sum(data_MM3['count'])))

#%% Create data_MM4
groups = ['ec_number','kineticlaw_type','substrate','organism','ph','buffer']

P=3 # Minimum number of values in each group
islogtransf = True # True for pKm and False for Km
cor_units = False # True if inputting data_cor_units (all kinetic law types)
data_all4 = make_groups(data_MM_pregroup,groups,islogtransf,P,cor_units)
data_all4 = data_all4.reset_index(drop=True)

nw_idx4 = organize_index(data_all4)
data_MM4 = data_all4.reindex(nw_idx4)
data_MM4 = data_MM4.reset_index(drop=True)

print('Number of measurements = {}'.format(sum(data_MM4['count'])))

#%% Create data_MM5
groups = ['ec_number','kineticlaw_type','substrate','organism','temperature','buffer']

P=3 # Minimum number of values in each group
islogtransf = True # True for pKm and False for Km
cor_units = False # True if inputting data_cor_units (all kinetic law types)
data_all5 = make_groups(data_MM_pregroup,groups,islogtransf,P,cor_units)
data_all5 = data_all5.reset_index(drop=True)

nw_idx5 = organize_index(data_all5)
data_MM5 = data_all5.reindex(nw_idx5)
data_MM5 = data_MM5.reset_index(drop=True)

print('Number of measurements = {}'.format(sum(data_MM5['count'])))

#%% Create data_MM6
groups = ['ec_number','substrate','organism','ph','temperature','buffer']

P=3 # Minimum number of values in each group
islogtransf = True # True for pKm and False for Km
cor_units = True # True if inputting data_cor_units (all kinetic law types)
data_all6 = make_groups(data_cor_units,groups,islogtransf,P,cor_units)
data_all6 = data_all6.reset_index(drop=True)

nw_idx6 = organize_index(data_all6)
data_MM6 = data_all6.reindex(nw_idx6)
data_MM6 = data_MM6.reset_index(drop=True)

print('Number of measurements = {}'.format(sum(data_MM6['count'])))
"""
#%% Create data_MM7
groups = ['ec_number','kineticlaw_type','substrate','ph','temperature','buffer']

P=3 # Minimum number of values in each group
islogtransf = True # True for pKm and False for Km
cor_units = False # True if inputting data_cor_units (all kinetic law types)
data_all7 = make_groups(data_MM_pregroup,groups,islogtransf,P,cor_units)
data_all7 = data_all7.reset_index(drop=True)

nw_idx7 = organize_index(data_all7)
data_MM7 = data_all7.reindex(nw_idx7)
data_MM7 = data_MM7.reset_index(drop=True)

print('Number of measurements = {}'.format(sum(data_MM7['count'])))

#%% Loop function but cycle through the different data_MM's
# data_MMs = [data_MM3,data_MM4,data_MM5,data_MM6,data_MM7]
data_MMs = [data_MM7]

s_log_scale = [1, 3, 6, 10, 30, 60, 100, 300, 600, 1000]
columns = ['Num Pairs','MUE','MUE std','MDUE','MDUE std','Stdev','Stdev std','R2 Pearson',\
           'R2 Pearson std','R2 Pearson Max','R2 Pearson Max std','Lower',\
           'Lower std','Upper','Upper std']
# loop_data_MRC3 = pd.DataFrame(index=s_log_scale,columns=columns)
# loop_data_MRC4 = pd.DataFrame(index=s_log_scale,columns=columns)
# loop_data_MRC5 = pd.DataFrame(index=s_log_scale,columns=columns)
# loop_data_MRC6 = pd.DataFrame(index=s_log_scale,columns=columns)
loop_data_MRC7 = pd.DataFrame(index=s_log_scale,columns=columns)

# Total stdev for entire dataset
stdevtot = np.std(data_MM_pregroup['pKm'])
# count = 0

for dat in data_MMs:
    ind_pkm = dat.columns.get_loc('Km')
    ind_EC = dat.columns.get_loc('ec_number')
    L = len(dat)
    Max = 0

    df_pairs,L_pairs = randomvals_and_diff(dat,ind_pkm,ind_EC)

    """ When you only want to consider MRC2 sets with multiple pubs
    data_MM2 = findpubs(data_MM_MRC1)
    data_MM2,df_fixes = findfirstauthors(data_MM2)
    data_fauthpub = data_MM2[data_MM2['Number of First Authors']<data_MM2['Number of Pubs']]
    df_multpubs,df_pairs = sortbypubs(data_MM2,False,data_fauthpub)
    df_pairs = getallpairs(df_pairs)
    L_pairs = len(df_pairs)
    """ ###########################################################

    MUE_evol,MDUE_evol,stdev_evol,R2p_evol,R2pmax_evol,Dstdev_evol = \
        measuresofquality(stdevtot,df_pairs)

    mean_evol,lower_evol,upper_evol = mean_confidence_interval(MUE_evol,stdev_evol,len(df_pairs))

    for i in range(len(s_log_scale)):
        Min = Max
        Max = s_log_scale[i]
        for k in range(Min,Max):

            """ When you only want to consider MRC2 sets with multiple pubs
            data_MM2 = findpubs(data_MM_MRC1)
            data_MM2,df_fixes = findfirstauthors(data_MM2)
            data_fauthpub = data_MM2[data_MM2['Number of First Authors']<data_MM2['Number of Pubs']]
            df_multpubs,df_pairs = sortbypubs(data_MM2,False,data_fauthpub)
            df_pairs = getallpairs(df_pairs)
            L_pairs = len(df_pairs)
            """ ###########################################################

            df_pairs,L_pairs = randomvals_and_diff(dat,ind_pkm,ind_EC)

            [MUE,MDUE,stdev,R2p,R2pmax,Dstdev] = measuresofquality(stdevtot,df_pairs)

            [mean,lower,upper] = mean_confidence_interval(MUE,stdev,len(df_pairs))

            """ Storing Values """
            MUE_evol = np.vstack((MUE_evol,MUE))

            MDUE_evol = np.vstack((MDUE_evol,MDUE))

            stdev_evol = np.vstack((stdev_evol,stdev))

            lower_evol = np.vstack((lower_evol,lower))

            upper_evol = np.vstack((upper_evol,upper))

            R2p_evol = np.vstack((R2p_evol,R2p))

            R2pmax_evol = np.vstack((R2pmax_evol,R2pmax))

            print('{} Loop'.format(k))

        """ Averages """
        MUE_ave = np.mean(MUE_evol)

        MDUE_ave = np.mean(MDUE_evol)

        stdev_ave = np.mean(stdev_evol)

        lower_ave = np.mean(lower_evol)

        upper_ave = np.mean(upper_evol)

        R2p_ave = np.mean(R2p_evol)

        R2pmax_ave = np.mean(R2pmax_evol)

        """ Standard Deviations """
        MUE_dev = np.std(MUE_evol)

        MDUE_dev = np.std(MDUE_evol)

        stdev_dev = np.std(stdev_evol)

        lower_dev = np.std(lower_evol)

        upper_dev = np.std(upper_evol)

        R2p_dev = np.std(R2p_evol)

        R2pmax_dev = np.std(R2pmax_evol)

        # Save looping data
        MRC = [L,MUE_ave,MUE_dev,MDUE_ave,MDUE_dev,stdev_ave,stdev_dev,R2p_ave,R2p_dev,\
                R2pmax_ave,R2pmax_dev,lower_ave,lower_dev,upper_ave,upper_dev]

        for l in range(len(MRC)):
            # if count == 0:
            #     loop_data_MRC3.iloc[i,l] = MRC[l]
            # elif count == 1:
            #     loop_data_MRC4.iloc[i,l] = MRC[l]
            # elif count == 2:
            #     loop_data_MRC5.iloc[i,l] = MRC[l]
            # elif count == 3:
            #     loop_data_MRC6.iloc[i,l] = MRC[l]
            # elif count == 4:
            loop_data_MRC7.iloc[i,l] = MRC[l]

        if i > 0:
            print('Done with looping {} times'.format(s_log_scale[i]))
        else:
            print('Done with looping 1 time')

    # count = count + 1

# loop_data_MRC3.to_csv('figures/loop_data_MRC3_rev.csv')
# loop_data_MRC4.to_csv('figures/loop_data_MRC4_rev.csv')
# loop_data_MRC5.to_csv('figures/loop_data_MRC5_rev.csv')
# loop_data_MRC6.to_csv('figures/loop_data_MRC6_rev.csv')
loop_data_MRC7.to_csv('figures/loop_data_MRC7_rev.csv')

print('3 = Buffer, 4 = Temperature, 5 = pH, 6 = Kinetic Law Type, 7 = Organism')
