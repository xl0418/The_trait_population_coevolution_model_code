import os,sys
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_ms/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
import pandas as pd
import csv

def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


#full tree
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'
files = dir_path + 'treedata/'

# estimation results
fileMS_name = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/BWest_t20000_h1'
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

# Check the fitness
tp_fitness = fitness[iterations-1,:int(total_population/2)]
dr_fitness = fitness[iterations-1,int(total_population/2):]

# estimates for TP
bestTP = fit_index[np.where(fit_index<10000)[0]]
gamma_TP_mean = np.mean(ms_data['gamma_data_TP'][iterations,bestTP])
a_TP_mean = np.mean(ms_data['a_data_TP'][iterations,bestTP])
nu_TP_mean = np.mean(ms_data['nu_data_TP'][iterations,bestTP])

gamma = gamma_TP_mean
a = a_TP_mean
nu = nu_TP_mean
K = 10e8

with open(files+'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(dir_path+'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:,1]
lengthdata_array = np.array(lengthdata)
length_index = []
for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)
length = 10**logTL
meantrait = 32.50540571860487
particle_size = 1000
td = DVTreeData(path=files, scalar=20000)



param = DVParamLiang(gamma=gamma, a=a, K=K, nu=nu, r=1, theta=meantrait, Vmax=1, inittrait=meantrait, initpop=500,
                initpop_sigma = 10.0, break_on_mu=False)
params = np.tile(param, (particle_size, 1))       # duplicate


predictsim = dvcpp.DVSimLiang(td, params)

valid = np.where(predictsim['sim_time'] == td.sim_evo_time)[0]

Z = predictsim['Z'][valid]
# i, j = argsort2D(Z)
# Z = Z[i, j]
# V = pop['V'][valid][i, j]
Z = np.nan_to_num(Z)

Z_df = pd.DataFrame(Z)
Z_df.to_csv(files+'predictsimTP.csv',sep=',',index=False)

