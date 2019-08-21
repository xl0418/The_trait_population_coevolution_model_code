#%%
import os
import numpy as np
import platform
import pandas as pd
import matplotlib.pyplot as plt
# treeno_vec = [i for i in range(4,7)]
# test = 4
# gindicator = 4
# aindicator = 5
#

if platform.system() == 'Windows':
    dire = 'C:/Liang/Googlebox'
elif platform.system() == 'Darwin':
    dire = '/Users/dudupig/GoogleDrive'
# file = dire + '/Research/Project2/smc_newdata/test%d/smc%dg%da.npy' % (test,gindicator,aindicator)
file = dire + '/Research/Project2/BaleenWhales/Est/BW_smc_res_actualT2wsL.npy'

assert os.path.isfile(file),"%s doesn't exist!" % file
para_data = np.load(file).item()

para_data_gamma_df = pd.DataFrame(para_data['gamma'])
para_data_a_df = pd.DataFrame(para_data['a'])

# filesmcg_name = 'C:\\Liang\\PhdIntroProject2\\smcplot\\smcdatag%d.csv' % test
# filesmca_name = 'C:\\Liang\\PhdIntroProject2\\smcplot\\smcdataa%d.csv' % test
filesmcg_name = dire + '/Research/Project2/BaleenWhales/Est/est_g.csv'
filesmca_name =  dire + '/Research/Project2/BaleenWhales/Est/est_a.csv'

para_data_gamma_df.to_csv(filesmcg_name, encoding='utf-8', index=False)
para_data_a_df.to_csv(filesmca_name, encoding='utf-8', index=False)



# For model selection
generating = 'TP'
# fileMS = 'C:/Liang/Googlebox/Research/Project2/modelsele/example1/modelsele%s.npy' % generating
# fileMS = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/50itertest/BWMS_50_TPDR3.npy'
fileMS_name = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/BWMSt40000_d1_f1'
fileMS = fileMS_name + '.npy'


assert os.path.isfile(fileMS),"%s doesn't exist!" % fileMS
ms_data = np.load(fileMS).item()
modeldata = ms_data['model_data'][1:,:]
iterations = modeldata.shape[0]

fitness = ms_data['fitness']
total_population = modeldata.shape[1]
bestpercent_index = int(total_population // 2)
bestmodel = np.zeros(shape=[iterations,bestpercent_index-1])
for g in range(iterations):
    q5 = np.argsort(fitness[g, :])[-bestpercent_index]  # best 25%
    fit_index = np.where(fitness[g, :] > fitness[g, q5])[0]
    bestmodel[g,:] = modeldata[g,fit_index]


bestmodel_df = pd.DataFrame(bestmodel)

#
# modeldata_df = pd.DataFrame(ms_data['model_data'])
# fitness_df = pd.DataFrame(ms_data['fitness'])
#
# filesmcms_name = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/50itertest/BWMS_50.csv'
# filesmcfit_name = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/50itertest/BWMS_50fit.csv'
# bestmodel_name = 'C:/Liang/Googlebox/Research/Project2/modelsele/example1/bestmodel%s.csv' % generating
# modeldata_df.to_csv(filesmcms_name, encoding='utf-8', index=False)
# fitness_df.to_csv(filesmcfit_name, encoding='utf-8', index=False)


bestmodel_name = fileMS_name + '.csv'
bestmodel_df.to_csv(bestmodel_name, encoding='utf-8', index=False)


# Check the fitness
tp_fitness = fitness[iterations-1,:int(total_population/2)]
dr_fitness = fitness[iterations-1,int(total_population/2):]

# diff_norm = np.linalg.norm(Z - obsZ, axis=1)
# q5 = np.argsort(diff_norm)[int(total_population // 5)]  # best 20%
# fit_index = np.where(diff_norm < diff_norm[q5])[0]
#
# modelTPperc = len(np.where(propose_model[fit_index] == 0)[0]) / len(fit_index)
# modeldrperc = len(np.where(propose_model[fit_index] == 1)[0]) / len(fit_index)
# estimates for TP
bestTP = fit_index[np.where(fit_index<10000)[0]]
gamma_TP_mean = np.mean(ms_data['gamma_data_TP'][iterations,bestTP])
a_TP_mean = np.mean(ms_data['a_data_TP'][iterations,bestTP])
nu_TP_mean = np.mean(ms_data['nu_data_TP'][iterations,bestTP])

# estimates for DR
bestDR = fit_index[np.where(np.logical_and(fit_index>=10000, fit_index<20000))]-10000
gamma_DR_mean = np.mean(ms_data['gamma_data_DR'][iterations,bestDR])
a_DR_mean = np.mean(ms_data['a_data_DR'][iterations,bestDR])
m_DR_mean = np.mean(ms_data['m_data_DR'][iterations,bestDR])


# estimate for NH
bestNH = fit_index[np.where(fit_index>=1000)]-1000
gamma_NH_mean = np.mean(ms_data['gamma_data_NH'][iterations,bestNH])
m_NH_mean = np.mean(ms_data['m_data_NH'][iterations,bestNH])

