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
data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/result_0617/'

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
        fit_sorted = sorted(pick_fitness)
        order = [np.where(fit_sorted==pick_fitness[i])[0][0] for i in range(len(fit_index))]
        # last_fitness = est_data['fitness'][generation - 1]
        # q5 = np.argsort(last_fitness)[-population // 20]  # best 5%
        # fit_index = np.where(last_fitness > last_fitness[q5])[0]

        # mean of all samples
        # gamma_mean = np.mean(est_data['gamma'][generation-1])
        # a_mean = np.mean(est_data['a'][generation-1])
        # nv_mean = np.mean(est_data['nu'][generation-1])
        # mean of the top 5% samples
        fitness_level = np.repeat(['1','2','3','4'],250)
        color_fitness = fitness_level[order]
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
                                       'fitness':color_fitness,'theta':theta_vec})
        else:
            est_para_df = pd.DataFrame({'gamma':gamma_vec,'alpha':a_vec,'nu':nu_vec,'vm':vm_vec,
                                        'fitness':color_fitness})
        if 'theta' in est_data:
            pplot = sns.pairplot(est_para_df,hue='fitness',vars=["gamma", "alpha", "nu", "vm",'theta'],
        kind="reg",markers="+")
        else:
            pplot = sns.pairplot(est_para_df, hue='fitness', vars=["gamma", "alpha", "nu", "vm"],
                                 kind="reg", markers="+")

        plt.suptitle('s=%i,$h^2$=%.1f' % ((int(timescaling), heritability)),
                     size=15,x=0.53,y=1)
        #
        # corr = est_para_df.corr()
        # # Generate a mask for the upper triangle
        # mask = np.zeros_like(corr, dtype=np.bool)
        # mask[np.triu_indices_from(mask)] = True
        #
        # # Set up the matplotlib figure
        # f, ax = plt.subplots(figsize=(11, 9))
        #
        # # Generate a custom diverging colormap
        # cmap = sns.diverging_palette(220, 10, as_cmap=True)
        #
        # # Draw the heatmap with the mask and correct aspect ratio
        # sns.heatmap(corr, mask=mask, cmap=cmap, vmax=1, center=0,
        #             square=True, linewidths=.5, cbar_kws={"shrink": .5})
        # plt.title('s=%i,$h^2$=%.1f' % ((int(timescaling), heritability)))

        #
        # figsave = data_dir+ 'corr_t%i_h%i.png' % (int(timescaling_index), int(heritability_index))
        # f.savefig(figsave)
        # plt.close(f)
        figsave = data_dir+ 'corr_t%i_h%i.png' % (int(timescaling_index), int(heritability_index))
        plt.savefig(figsave)
        plt.close()