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
from functions import mean_confidence_interval,measuresofquality, \
    randomvals_and_diff, get_pairs_tot
    
data_MM_pregroup = pd.read_csv('data_MM_pregroup_pubs_auth.csv',index_col=0)
data_MM = pd.read_csv('data_MM.csv')
data_MM_MRC1 = pd.read_csv('data_MM_MRC1.csv')

df_pairs_tot = get_pairs_tot(data_MM_pregroup)
df_pairs = randomvals_and_diff(data_MM)
df_pairs_MRC1 = randomvals_and_diff(data_MM_MRC1)

MUE_tot_evol,MDUE_tot_evol,stdev_tot_evol,R2p_tot_evol,R2pmax_tot_evol,Dstdev_tot_evol = \
    measuresofquality(data_MM_pregroup,df_pairs_tot)
MUE_evol,MDUE_evol,stdev_evol,R2p_evol,R2pmax_evol,Dstdev_evol = \
    measuresofquality(data_MM_pregroup,df_pairs)
MUE_MRC1_evol,MDUE_MRC1_evol,stdev_MRC1_evol,R2p_MRC1_evol,R2pmax_MRC1_evol,Dstdev_MRC1_evol = \
    measuresofquality(data_MM_pregroup,df_pairs_MRC1)

mean_tot_evol,lower_tot_evol,upper_tot_evol = mean_confidence_interval(MUE_tot_evol,stdev_tot_evol,len(df_pairs_tot))
mean_evol,lower_evol,upper_evol = mean_confidence_interval(MUE_evol,stdev_evol,len(df_pairs))
mean_MRC1_evol,lower_MRC1_evol,upper_MRC1_evol = mean_confidence_interval(MUE_MRC1_evol,stdev_MRC1_evol,len(df_pairs_MRC1))

k = 0
Max = 100
while k < Max:
    
    # df_pairs_tot = get_pairs_tot(data_MM_pregroup)
    df_pairs = randomvals_and_diff(data_MM)
    df_pairs_MRC1 = randomvals_and_diff(data_MM_MRC1)
    
    # MUE_tot,MDUE_tot,stdev_tot,R2p_tot,R2pmax_tot,Dstdev_tot = \
    #     measuresofquality(data_MM_pregroup,df_pairs_tot)
    MUE,MDUE,stdev,R2p,R2pmax,Dstdev = \
        measuresofquality(data_MM_pregroup,df_pairs)
    MUE_MRC1,MDUE_MRC1,stdev_MRC1,R2p_MRC1,R2pmax_MRC1,Dstdev_MRC1 = \
        measuresofquality(data_MM_pregroup,df_pairs_MRC1)
    
    # mean_tot,lower_tot,upper_tot = mean_confidence_interval(MUE_tot,stdev_tot,len(df_pairs_tot))
    mean,lower,upper = mean_confidence_interval(MUE,stdev,len(df_pairs))
    mean_MRC1,lower_MRC1,upper_MRC1 = mean_confidence_interval(MUE_MRC1,stdev_MRC1,len(df_pairs_MRC1))
    
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
    lower_evol = np.vstack((lower_evol,stdev))
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
    
    """ Averages """
    # MUE_tot_ave = MUE_tot_evol.mean(0)
    MUE_ave = MUE_evol.mean(0)
    MUE_MRC1_ave = MUE_MRC1_evol.mean(0)
    
    # MDUE_tot_ave = MDUE_tot_evol.mean(0)
    MDUE_ave = MDUE_evol.mean(0)
    MDUE_MRC1_ave = MDUE_MRC1_evol.mean(0)
    
    # stdev_tot_ave = stdev_tot_evol.mean(0)
    stdev_ave = stdev_evol.mean(0)
    stdev_MRC1_ave = stdev_MRC1_evol.mean(0)
    
    # lower_tot_ave = lower_tot_evol.mean(0)
    lower_ave = lower_evol.mean(0)
    lower_MRC1_ave = lower_MRC1_evol.mean(0)
    
    # upper_tot_ave = upper_tot_evol.mean(0)
    upper_ave = upper_evol.mean(0)
    upper_MRC1_ave = upper_MRC1_evol.mean(0)
    
    # R2p_tot_ave = R2p_tot_evol.mean(0)
    R2p_ave = R2p_evol.mean(0)
    R2p_MRC1_ave = R2p_MRC1_evol.mean(0)
    
    # R2pmax_tot_ave = R2pmax_tot_evol.mean(0)
    R2pmax_ave = R2pmax_evol.mean(0)
    R2pmax_MRC1_ave = R2pmax_MRC1_evol.mean(0)
    
    """ Standard Deviations """
    # MUE_tot_dev = MUE_tot_evol.std(0)
    MUE_dev = MUE_evol.std(0)
    MUE_MRC1_dev = MUE_MRC1_evol.std(0)
    
    # MDUE_tot_dev = MDUE_tot_evol.std(0)
    MDUE_dev = MDUE_evol.std(0)
    MDUE_MRC1_dev = MDUE_MRC1_evol.std(0)
    
    # stdev_tot_dev = stdev_tot_evol.std(0)
    stdev_dev = stdev_evol.std(0)
    stdev_MRC1_dev = stdev_MRC1_evol.std(0)
    
    # lower_tot_dev = lower_tot_evol.std(0)
    lower_dev = lower_evol.std(0)
    lower_MRC1_dev = lower_MRC1_evol.std(0)
    
    # upper_tot_dev = upper_tot_evol.std(0)
    upper_dev = upper_evol.std(0)
    upper_MRC1_dev = upper_MRC1_evol.std(0)
    
    # R2p_tot_dev = R2p_tot_evol.std(0)
    R2p_dev = R2p_evol.std(0)
    R2p_MRC1_dev = R2p_MRC1_evol.std(0)
    
    # R2pmax_tot_dev = R2pmax_tot_evol.std(0)
    R2pmax_dev = R2pmax_evol.std(0)
    R2pmax_MRC1_dev = R2pmax_MRC1_evol.std(0)
    
    index = ['MRC2','MRC1']

    loop_data = pd.DataFrame({'Loops': [Max,Max],
    'MUE': [MUE_ave,MUE_MRC1_ave],
    'MUE std': [MUE_dev,MUE_MRC1_dev], 
    'MDUE': [MDUE_ave,MDUE_MRC1_ave],
    'MDUE std': [MDUE_dev,MDUE_MRC1_dev],
    'Stdev': [stdev_ave,stdev_MRC1_ave],
    'Stdev std': [stdev_dev,stdev_MRC1_dev],
    'R2 Pearson': [R2p_ave,R2p_MRC1_ave],
    'R2 Pearson std': [R2p_dev,R2p_MRC1_dev],
    'R2 Pearson Max': [R2pmax_ave,R2pmax_MRC1_ave],
    'R2 Pearson Max std': [R2pmax_dev,R2pmax_MRC1_dev],
    'Lower':[lower_ave,lower_MRC1_ave],
    'Lower std':[lower_dev,lower_MRC1_dev],
    'Upper':[upper_ave,upper_MRC1_ave],
    'Upper std':[upper_dev,upper_MRC1_dev]},index=index)
