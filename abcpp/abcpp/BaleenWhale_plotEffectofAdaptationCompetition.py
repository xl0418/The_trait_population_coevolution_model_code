import os,sys
sys.path.append('C:/Liang/abcpp_ms2/abcpp')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
import pandas as pd
import csv
import matplotlib.pyplot as plt
import seaborn as sns

def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def competitions(a, zi, nj,K):
    T = zi[:, np.newaxis] - zi  # trait-distance matrix (via 'broadcasting')
    t1 = np.exp(-a * T ** 2) * nj/K
    t2 = (2 * a) * T
    t3 = np.exp(-a * T ** 2)
    pairwis_com = t2 * t3
    sum_com = np.sum(t2 * t1, axis=1)
    return pairwis_com, sum_com,T

def adaptation(gamma,zi,theta):
    return 2*gamma * (theta - zi)


dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'
data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/'

obs_file = dir_path + 'treedata/'

with open(obs_file+'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(dir_path+'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:,1]
lengthdata_array = np.array(lengthdata)
length_index = []

distance_list = []
valid_length = []
timescaling_pair_list = []
heritability_pair_list =[]
timescaling_sum_list = []
heritability_sum_list =[]
pairwise_com_list = []
pairwise_deltaZ_list = []
a_list = []
sum_com_list = []
sum_adaptation_list = []
Z_list = []
N_list = []
count = 0
K=10e5
nu=1/(100*K)

for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)
emp_Z = 10**logTL
obsZ = sorted(emp_Z)
meantrait = np.mean(emp_Z)
timescale_vec = [20000,40000,60000,80000]
heritability_vec = [1,0.5]
for timescaling_index in range(4):
    for heritability_index in range(2):
        print(count)
        count += 1
        timescaling = timescale_vec[timescaling_index]
        heritability = heritability_vec[heritability_index]

        data_name = data_dir + 'BWest_t%i_h%i.npy' % (int(timescaling),int(heritability_index+1))
        est_data = np.load(data_name).item()
        generation = len(est_data['gamma'])
        population = len(est_data['gamma'][0])

        last_fitness = est_data['fitness'][generation-1]
        q5 = np.argsort(last_fitness)[-population // 20]  # best 5%
        fit_index = np.where(last_fitness > last_fitness[q5])[0]

        # mean of all samples
        gamma_mean = np.mean(est_data['gamma'][generation-1][fit_index])
        a_mean = np.mean(est_data['a'][generation-1][fit_index])
        nv_mean = np.mean(est_data['nu'][generation-1][fit_index])
        vm_mean = np.mean(est_data['vm'][generation-1][fit_index])
        a_list.append(a_mean)
        td = DVTreeData(path=obs_file, scalar=timescaling)

        obs_param = DVParamLiang(gamma=gamma_mean, a=a_mean, K=K, h=heritability, nu=nv_mean, r=1, theta=meantrait, V00=.5, V01=.5, Vmax=vm_mean,
                                 inittrait=meantrait, initpop=1e5,
                                 initpop_sigma=10.0, break_on_mu=False)
        for replicate in range(1000):
            print(replicate)

            simresult = dvcpp.DVSimLiang(td, obs_param)
            if simresult['sim_time'] == td.sim_evo_time:
                pic = 0
                break
            else:
                pic=1


        # mean of the top 5% samples
        Z_vec = simresult['Z']
        N_vec = simresult['N']

        pairwise_com, sum_com,delta_Z = competitions(a = a_mean,zi = Z_vec,nj = N_vec,K = K)
        sum_adaptation = adaptation(gamma = gamma_mean,zi = Z_vec,theta = meantrait)

        unique_pairwise_com = pairwise_com.flatten()
        unique_pairwise_delta_Z = delta_Z.flatten()
        pairwise_com_list.append(unique_pairwise_com)
        pairwise_deltaZ_list.append(unique_pairwise_delta_Z)
        Z_list.append(Z_vec)
        N_list.append(N_vec)

        sum_com_list.append(sum_com)
        sum_adaptation_list.append(sum_adaptation)

        timescaling_pair_list.append(np.repeat(timescaling,len(unique_pairwise_com)))
        heritability_pair_list.append(np.repeat(heritability,len(unique_pairwise_com)))
        timescaling_sum_list.append(np.repeat(timescaling,len(sum_com)))
        heritability_sum_list.append(np.repeat(heritability,len(sum_adaptation)))

pairwise_com_list_flat = [item for sublist in pairwise_com_list for item in sublist]
pairwise_deltaZ_list_flat = [item for sublist in pairwise_deltaZ_list for item in sublist]
sum_com_list_flat = [item for sublist in sum_com_list for item in sublist]
sum_adaptation_list_flat = [item for sublist in sum_adaptation_list for item in sublist]

timescaling_pair_list_flat = [item for sublist in timescaling_pair_list for item in sublist]
heritability_pair_list_flat = [item for sublist in heritability_pair_list for item in sublist]
timescaling_sum_list_flat = [item for sublist in timescaling_sum_list for item in sublist]
heritability_sum_list_flat= [item for sublist in heritability_sum_list for item in sublist]

Z_list_flat = [item for sublist in Z_list for item in sublist]
N_list_flat = [item for sublist in N_list for item in sublist]

pair_pd = {'pairwise_com':pairwise_com_list_flat,'pairwise_distance':pairwise_deltaZ_list_flat,
           'timescale':timescaling_pair_list_flat,'heritability':heritability_pair_list_flat}
sum_pd = {'sum_com':sum_com_list_flat,'sum_ada':sum_adaptation_list_flat,
           'timescale':timescaling_sum_list_flat,'heritability':heritability_sum_list_flat}
Z_N_pd = {'Z':Z_list_flat,'N':N_list_flat,
           'timescale':timescaling_sum_list_flat,'heritability':heritability_sum_list_flat}

pair_df = pd.DataFrame(pair_pd)
sum_df = pd.DataFrame(sum_pd)
Z_N_df = pd.DataFrame(Z_N_pd)

sns.set(style="ticks")
def reg_com(a, x):
    return 2*a*x*np.exp(-a*x**2)

xlim = [-1000,1000]

# empirical strength of competition

def emp_com(a, zi):
    T = zi[:, np.newaxis] - zi  # trait-distance matrix (via 'broadcasting')
    t2 = (2 * a) * T
    t3 = np.exp(-a * T ** 2)
    pairwis_com = t2 * t3
    return pairwis_com,T


# pairwise competition against trait distance plot
comdis = np.arange(xlim[0], xlim[1], 0.1)

fig_pair_com = sns.FacetGrid(pair_df, col="heritability", row="timescale",margin_titles=True,xlim=(xlim[0], xlim[1]))
# fig_pair_com.map(plt.scatter, "pairwise_distance", "pairwise_com", alpha=.5)
fig_pair_com.set_ylabels('Competition strength')
fig_pair_com.set_xlabels('Trait distance')

fig_count = 0
for ax in fig_pair_com.axes.flat:
    comstre = reg_com(a_list[fig_count],comdis)
    emp_pair_strength, emp_pair_dis = emp_com(a_list[fig_count], zi=emp_Z)
    ax.scatter(emp_pair_dis.flatten(),emp_pair_strength.flatten(),c='red',alpha = 0.5)
    ax.plot(comdis, comstre, c=".2", ls="--")
    fig_count +=1
fig_pair_com.add_legend();


# total adaptation against competition plot
com_ada_xlim = [np.min(sum_df['sum_com']),np.max(sum_df['sum_com'])]
fig_sum_ada_com = sns.FacetGrid(sum_df, col="heritability", row="timescale",margin_titles=True,xlim=(com_ada_xlim[0], com_ada_xlim[1]))
fig_sum_ada_com.map(plt.scatter, "sum_com", "sum_ada", alpha=.5)
fig_sum_ada_com.set_ylabels('Adaptation strength')
fig_sum_ada_com.set_xlabels('Competition distance')

for ax in fig_sum_ada_com.axes.flat:
    x_data = np.arange(com_ada_xlim[0],com_ada_xlim[1],(com_ada_xlim[1]-com_ada_xlim[0])/100)
    y_data = -x_data
    ax.plot(x_data, y_data, c=".2", ls="--")
    fig_count +=1
fig_sum_ada_com.add_legend();



# trait against population plot
Z_N_xlim = [np.min(Z_N_df['Z']),np.max(Z_N_df['Z'])]
fig_Z_N = sns.FacetGrid(Z_N_df, col="heritability", row="timescale",margin_titles=True,xlim=(Z_N_xlim[0], Z_N_xlim[1]))
fig_Z_N.map(plt.scatter, "Z", "N", alpha=.5)
fig_Z_N.set_ylabels('Traits')
fig_Z_N.set_xlabels('Abundances')


fig_Z_N.add_legend();