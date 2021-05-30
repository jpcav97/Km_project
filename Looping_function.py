#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 15:52:13 2021

@author: josephcavataio
"""
##############################################################################
########################## Looping Function ##################################
##############################################################################
import pandas as pd
import numpy as np
import random
from functions import mean_confidence_interval,measuresofquality, \
    randomvals_and_diff, get_pairs_tot, make_groups, organize_index

#%% Import raw data
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

#%% Create data_MM
groups = ['ec_number','kineticlaw_type','substrate','organism','ph','temperature','buffer']

P=3 # Minimum number of values in each group
islogtransf = False # True for pKm and False for Km
data_all = make_groups(data_MM_pregroup,groups,islogtransf,P)
data_all = data_all.reset_index(drop=True)

nw_idx = organize_index(data_all)
data_MM = data_all.reindex(nw_idx)
data_MM = data_MM.reset_index(drop=True)

print('Number of measurements = {}'.format(sum(data_MM['count'])))

#%% Create data_MM_MRC1
groups = ['ec_number','kineticlaw_type','substrate','organism']

data_MRC1 = make_groups(data_MM_pregroup,groups,True,P)
data_MRC1 = data_MRC1.reset_index(drop=True)

nw_idx = organize_index(data_MRC1)
data_MM_MRC1 = data_MRC1.reindex(nw_idx)
data_MM_MRC1 = data_MM_MRC1.reset_index(drop=True)
L2 = len(data_MM_MRC1)

print('Number of measurements = {}'.format(sum(data_MM_MRC1['count'])))

#%% Looping
ind_pkm = data_MM.columns.get_loc('Km')
ind_EC = data_MM.columns.get_loc('ec_number')

ind_pkm_MRC1 = data_MM_MRC1.columns.get_loc('Km')
ind_EC_MRC1 = data_MM_MRC1.columns.get_loc('ec_number')

s_log_scale = [1, 3, 6, 10, 30, 60, 100, 300, 600, 1000]
columns = ['MUE','MUE std','MDUE','MDUE std','Stdev','Stdev std','R2 Pearson',\
           'R2 Pearson std','R2 Pearson Max','R2 Pearson Max std','Lower',\
               'Lower std','Upper','Upper std']

loop_data_MRC1 = pd.DataFrame(index=s_log_scale,columns=columns)
loop_data_MRC2 = pd.DataFrame(index=s_log_scale,columns=columns)
half_loop_data_MRC1 = pd.DataFrame(index=s_log_scale,columns=columns)
half_loop_data_MRC2 = pd.DataFrame(index=s_log_scale,columns=columns)


for j in range(2):

    Max = 0

    df_pairs_tot,L_pairs_tot = get_pairs_tot(data_MM_pregroup)
    df_pairs,L_pairs = randomvals_and_diff(data_MM,ind_pkm,ind_EC)
    df_pairs_MRC1,L_pairs_MRC1 = randomvals_and_diff(data_MM_MRC1,ind_pkm_MRC1,ind_EC_MRC1)

    if j == 1:
        pairs_index_tot = list(df_pairs_tot.index)
        pairs_index = list(df_pairs.index)
        pairs_index_MRC1 = list(df_pairs_MRC1.index)

        rand_tot = random.sample(pairs_index_tot,round(L_pairs_tot/2))
        rand = random.sample(pairs_index,round(L_pairs/2))
        rand_MRC1 = random.sample(pairs_index_MRC1,round(L_pairs_MRC1/2))

        df_pairs_tot = df_pairs_tot.iloc[rand_tot]
        df_pairs = df_pairs.iloc[rand]
        df_pairs_MRC1 = df_pairs_MRC1.iloc[rand_MRC1]

    MUE_tot_evol,MDUE_tot_evol,stdev_tot_evol,R2p_tot_evol,R2pmax_tot_evol,Dstdev_tot_evol = \
        measuresofquality(data_MM_pregroup,df_pairs_tot)
    MUE_evol,MDUE_evol,stdev_evol,R2p_evol,R2pmax_evol,Dstdev_evol = \
        measuresofquality(data_MM_pregroup,df_pairs)
    MUE_MRC1_evol,MDUE_MRC1_evol,stdev_MRC1_evol,R2p_MRC1_evol,R2pmax_MRC1_evol,Dstdev_MRC1_evol = \
        measuresofquality(data_MM_pregroup,df_pairs_MRC1)

    mean_tot_evol,lower_tot_evol,upper_tot_evol = mean_confidence_interval(MUE_tot_evol,stdev_tot_evol,len(df_pairs_tot))
    mean_evol,lower_evol,upper_evol = mean_confidence_interval(MUE_evol,stdev_evol,len(df_pairs))
    mean_MRC1_evol,lower_MRC1_evol,upper_MRC1_evol = mean_confidence_interval(MUE_MRC1_evol,stdev_MRC1_evol,len(df_pairs_MRC1))

    for i in  range(len(s_log_scale)):
        Min = Max
        Max = s_log_scale[i]
        for k in np.arange(Min,Max):

            # df_pairs_tot = get_pairs_tot(data_MM_pregroup)
            df_pairs,L_pairs = randomvals_and_diff(data_MM,ind_pkm,ind_EC)
            df_pairs_MRC1,L_pairs_MRC1 = randomvals_and_diff(data_MM_MRC1,ind_pkm_MRC1,ind_EC_MRC1)

            # MUE_tot,MDUE_tot,stdev_tot,R2p_tot,R2pmax_tot,Dstdev_tot = \
            #     measuresofquality(data_MM_pregroup,df_pairs_tot)
            [MUE,MDUE,stdev,R2p,R2pmax,Dstdev] = \
                measuresofquality(data_MM_pregroup,df_pairs)
            [MUE_MRC1,MDUE_MRC1,stdev_MRC1,R2p_MRC1,R2pmax_MRC1,Dstdev_MRC1] = \
                measuresofquality(data_MM_pregroup,df_pairs_MRC1)

            # mean_tot,lower_tot,upper_tot = mean_confidence_interval(MUE_tot,stdev_tot,len(df_pairs_tot))
            [mean,lower,upper] = mean_confidence_interval(MUE,stdev,len(df_pairs))
            [mean_MRC1,lower_MRC1,upper_MRC1] = mean_confidence_interval(MUE_MRC1,stdev_MRC1,len(df_pairs_MRC1))

            """ Storing Values """
            # MUE_tot_evol = np.vstack((MUE_tot_evol,MUE_tot))
            MUE_evol = np.vstack((MUE_evol,MUE))
            MUE_MRC1_evol = np.vstack((MUE_MRC1_evol,MUE_MRC1))

            # MDUE_tot_evol = np.vstack((MDUE_tot_evol,MDUE_tot))
            MDUE_evol = np.vstack((MDUE_evol,MDUE))
            MDUE_MRC1_evol = np.vstack((MDUE_MRC1_evol,MDUE_MRC1))

            # stdev_tot_evol = np.vstack((stdev_tot_evol,stdev_tot))
            stdev_evol = np.vstack((stdev_evol,stdev))
            stdev_MRC1_evol = np.vstack((stdev_MRC1_evol,stdev_MRC1))

            # lower_tot_evol = np.vstack((lower_tot_evol,lower_tot))
            lower_evol = np.vstack((lower_evol,lower))
            lower_MRC1_evol = np.vstack((lower_MRC1_evol,lower_MRC1))

            # upper_tot_evol = np.vstack((upper_tot_evol,upper_tot))
            upper_evol = np.vstack((upper_evol,upper))
            upper_MRC1_evol = np.vstack((upper_MRC1_evol,upper_MRC1))

            # R2p_tot_evol = np.vstack((R2p_tot_evol,R2p_tot))
            R2p_evol = np.vstack((R2p_evol,R2p))
            R2p_MRC1_evol = np.vstack((R2p_MRC1_evol,R2p_MRC1))

            # R2pmax_tot_evol = np.vstack((R2pmax_tot_evol,R2pmax_tot))
            R2pmax_evol = np.vstack((R2pmax_evol,R2pmax))
            R2pmax_MRC1_evol = np.vstack((R2pmax_MRC1_evol,R2pmax_MRC1))

            print('{} Loop'.format(k))

        """ Averages """
        # MUE_tot_ave = MUE_tot_evol.mean(0)
        MUE_ave = np.mean(MUE_evol)
        MUE_MRC1_ave = np.mean(MUE_MRC1_evol)

        # MDUE_tot_ave = MDUE_tot_evol.mean(0)
        MDUE_ave = np.mean(MDUE_evol)
        MDUE_MRC1_ave = np.mean(MDUE_MRC1_evol)

        # stdev_tot_ave = stdev_tot_evol.mean(0)
        stdev_ave = np.mean(stdev_evol)
        stdev_MRC1_ave = np.mean(stdev_MRC1_evol)

        # lower_tot_ave = lower_tot_evol.mean(0)
        lower_ave = np.mean(lower_evol)
        lower_MRC1_ave = np.mean(lower_MRC1_evol)

        # upper_tot_ave = upper_tot_evol.mean(0)
        upper_ave = np.mean(upper_evol)
        upper_MRC1_ave = np.mean(upper_MRC1_evol)

        # R2p_tot_ave = R2p_tot_evol.mean(0)
        R2p_ave = np.mean(R2p_evol)
        R2p_MRC1_ave = np.mean(R2p_MRC1_evol)

        # R2pmax_tot_ave = R2pmax_tot_evol.mean(0)
        R2pmax_ave = np.mean(R2pmax_evol)
        R2pmax_MRC1_ave = np.mean(R2pmax_MRC1_evol)

        """ Standard Deviations """
        # MUE_tot_dev = MUE_tot_evol.std(0)
        MUE_dev = np.std(MUE_evol)
        MUE_MRC1_dev = np.std(MUE_MRC1_evol)

        # MDUE_tot_dev = MDUE_tot_evol.std(0)
        MDUE_dev = np.std(MDUE_evol)
        MDUE_MRC1_dev = np.std(MDUE_MRC1_evol)

        # stdev_tot_dev = stdev_tot_evol.std(0)
        stdev_dev = np.std(stdev_evol)
        stdev_MRC1_dev = np.std(stdev_MRC1_evol)

        # lower_tot_dev = lower_tot_evol.std(0)
        lower_dev = np.std(lower_evol)
        lower_MRC1_dev = np.std(lower_MRC1_evol)

        # upper_tot_dev = upper_tot_evol.std(0)
        upper_dev = np.std(upper_evol)
        upper_MRC1_dev = np.std(upper_MRC1_evol)

        # R2p_tot_dev = R2p_tot_evol.std(0)
        R2p_dev = np.std(R2p_evol)
        R2p_MRC1_dev = np.std(R2p_MRC1_evol)

        # R2pmax_tot_dev = R2pmax_tot_evol.std(0)
        R2pmax_dev = np.std(R2pmax_evol)
        R2pmax_MRC1_dev = np.std(R2pmax_MRC1_evol)


        # Save looping data
        MRC2 = [MUE_ave,MUE_dev,MDUE_ave,MDUE_dev,stdev_ave,stdev_dev,R2p_ave,R2p_dev,\
                R2pmax_ave,R2pmax_dev,lower_ave,lower_dev,upper_ave,upper_dev]

        MRC1 = [MUE_MRC1_ave,MUE_MRC1_dev,MDUE_MRC1_ave,MDUE_MRC1_dev,stdev_MRC1_ave,\
                stdev_MRC1_dev,R2p_MRC1_ave,R2p_MRC1_dev,R2pmax_MRC1_ave,\
                R2pmax_MRC1_dev,lower_MRC1_ave,lower_MRC1_dev,\
                upper_MRC1_ave,upper_MRC1_dev]

        for l in range(len(MRC1)):
            if j == 0:
                loop_data_MRC1.iloc[i,l] = MRC1[l]
                loop_data_MRC2.iloc[i,l] = MRC2[l]
            if j == 1:
                half_loop_data_MRC1.iloc[i,l] = MRC1[l]
                half_loop_data_MRC2.iloc[i,l] = MRC2[l]

        # Compare deviation of current Mean unsigned error to the MUE STD taken after 3 loops
        if i > 0:
            ind_MUE_std = loop_data_MRC2.columns.get_loc('MUE std')
            MUE_dev_compare = float(loop_data_MRC2.iloc[1,ind_MUE_std])
            if MUE_dev <= MUE_dev_compare/2.718:
                magic_number = i
                print("Number of runs to reduce STD by 1/2.718 is {}".format(s_log_scale[i]))

            print('Done with looping {} times'.format(s_log_scale[i]))
        else:
            print('Done with looping 1 time')

loop_data_MRC1.to_csv('figures/loop_data_MRC1.csv')
loop_data_MRC2.to_csv('figures/loop_data_MRC2.csv')
half_loop_data_MRC1.to_csv('figures/half_loop_data_MRC1.csv')
half_loop_data_MRC2.to_csv('figures/half_loop_data_MRC2.csv')
