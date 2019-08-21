import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import seaborn as sns

def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err


dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'
data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/results_0529/'

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
gamma_list = []
a_list = []
nu_list = []
vm_list = []
theta_list = []
distance_list = []
ratio_dis = []
valid_length = []
timescaling_list = []
heri_list = []
count = 0
particle_size = 100
K=1e6
nu=1/(100*K)
alphavalue = 0.8
cmapvalue = 'OrRd'
for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)

timescale_vec = [20000,40000,60000,80000,100000]
heritability_vec = [0.5,1]


for timescaling_index in range(5):
    for heritability_index in range(2):
        print(count)
        count += 1
        timescaling = timescale_vec[timescaling_index]
        heritability = heritability_vec[heritability_index]

        data_name = data_dir + 'BWest_t%i_h%i.npy' % (int(timescaling), int(heritability_index))
        est_data = np.load(data_name).item()
        generation = len(est_data['gamma'])
        population = len(est_data['gamma'][0])
        fitness = est_data['fitness'][-1]
        valid = est_data['valid']
        q5 = np.argsort(fitness)[-int(population // 20+1)]  # best 5%
        fit_index = np.where(fitness > fitness[q5])[0]
        pick_fitness = fitness[fit_index]

        gamma_vec = est_data['gamma'][-1][fit_index] #*1e8
        a_vec = est_data['a'][-1][fit_index]#*1e5
        nu_vec = est_data['nu'][-1][fit_index]#*1e4
        vm_vec = est_data['vm'][-1][fit_index]
        if 'theta' in est_data:
            theta_vec = est_data['theta'][-1][fit_index]

        gamma_vec = (gamma_vec-np.mean(gamma_vec))/np.sqrt(np.var(gamma_vec))
        a_vec = (a_vec-np.mean(a_vec))/np.sqrt(np.var(a_vec))
        nu_vec = (nu_vec-np.mean(nu_vec))/np.sqrt(np.var(nu_vec))
        vm_vec = (vm_vec-np.mean(vm_vec))/np.sqrt(np.var(vm_vec))
        # pick_fitness = (pick_fitness-np.mean(pick_fitness))/np.sqrt(np.var(pick_fitness))
        if 'theta' in est_data:
            theta_vec = (theta_vec-np.mean(theta_vec))/np.sqrt(np.var(theta_vec))
        inflow_mutation = 2*K*nu_vec*1e-5 * vm_vec/(1+4*K*nu_vec*1e-5)

        print('s=%i h2=%f 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f min(Vm) = %f' % (timescaling,
                                                                                        heritability,
            np.mean(gamma_vec),
            np.mean(a_vec), np.mean(nu_vec),
            np.mean(vm_vec), np.min(vm_vec)))
        print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f inflow = %f' % (
            np.var(gamma_vec),
            np.var(a_vec), np.var(nu_vec),
            np.var(vm_vec), np.var(inflow_mutation)))

        if 'theta' in est_data:
            est_para_df = pd.DataFrame({'gamma':gamma_vec,'alpha':a_vec,'nu':nu_vec,'vm':vm_vec,
                                       'theta':theta_vec})
            scatter_fig, axes = plt.subplots(5, 5, sharex=True, sharey=True,figsize=(9,9))
            for i in range(5):
                for j in range(5):
                    axes[i][j].scatter(est_para_df[est_para_df.columns[i]], est_para_df[est_para_df.columns[j]],
                                       marker='+', c=pick_fitness,
                                       alpha=alphavalue, cmap=cmapvalue)

            axes[4][0].set_ylabel('$\\theta$')
            axes[4][0].set_xlabel('$\gamma$')
            axes[4][1].set_xlabel('$\\alpha$')
            axes[4][2].set_xlabel('$\\nu$')
            axes[4][3].set_xlabel('$V_m$')
            axes[4][4].set_xlabel('$\\theta$')

        else:
            est_para_df = pd.DataFrame({'gamma':gamma_vec,'alpha':a_vec,'nu':nu_vec,'vm':vm_vec})
            scatter_fig, axes = plt.subplots(4, 4, sharex=True, sharey=True,figsize=(9,9))
            for i in range(4):
                for j in range(4):
                    axes[i][j].scatter(est_para_df[est_para_df.columns[i]], est_para_df[est_para_df.columns[j]],
                                       marker='+', c=pick_fitness,
                                       alpha=alphavalue, cmap=cmapvalue)
            axes[3][0].set_xlabel('$\gamma$')
            axes[3][1].set_xlabel('$\\alpha$')
            axes[3][2].set_xlabel('$\\nu$')
            axes[3][3].set_xlabel('$V_m$')

        axes[0][0].set_ylabel('$\gamma$')
        axes[1][0].set_ylabel('$\\alpha$')
        axes[2][0].set_ylabel('$\\nu$')
        axes[3][0].set_ylabel('$V_m$')


        plt.suptitle('s=%i,$h^2$=%.1f' % ((int(timescaling), heritability)),
                     size=15,x=0.53,y=1)

        figsave = data_dir+ 'scatter_t%i_h%i.png' % (int(timescaling_index), int(heritability_index))
        plt.savefig(figsave)
        plt.close()