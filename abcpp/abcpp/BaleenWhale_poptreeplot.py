import sys, os
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_ms5/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_py import DVSimTVP
from dvtraitsim_shared import DVTreeData, DVParamLiang
import numpy as np
import matplotlib.pyplot as plt
import csv

# let's try to find a true simulation:


#full tree and pruned tree directory

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

count = 0
particle_size = 100
K=1e6
nu=1/(100*K)

for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)
obsZ = sorted(10**logTL)

# for timescaling_index in range(len(timescale_vec)):
#     for heritability_index in range(len(heritability_vec)):
# print(count)
# count += 1
# timescaling = timescale_vec[timescaling_index]
# heritability = heritability_vec[heritability_index]

data_name = data_dir + 'modelselec4w.npy'
est_data = np.load(data_name).item()
population = int(len(est_data['model_data'][0])/3)
fitness = est_data['fitness'][-1]

q5_TVP = np.argsort(fitness[ :population])[-int(population // 200+1)]  # best 5%
q5_TV = np.argsort(fitness[ population:2*population])[-int(population // 200+1)]+population  # best 5%
q5_TVM = np.argsort(fitness[ 2*population:])[-int(population // 200+1)]+2*population  # best 5%

fit_index_TVP = np.where(fitness[ :population] > fitness[ q5_TVP])[0]
fit_index_TV = np.where(fitness[ population:2*population] > fitness[ q5_TV])[0]+population
fit_index_TVM = np.where(fitness[ 2*population:] > fitness[ q5_TVM])[0]+2*population

previous_bestfitted_index_TVP = fit_index_TVP
previous_bestfitted_index_TV = fit_index_TV - population
previous_bestfitted_index_TVM = fit_index_TVM- 2*population

top_pop = len(previous_bestfitted_index_TVP)

# last_fitness = est_data['fitness'][generation - 1]
# q5 = np.argsort(last_fitness)[-population // 20]  # best 5%
# fit_index = np.where(last_fitness > last_fitness[q5])[0]

# mean of all samples
# gamma_mean = np.mean(est_data['gamma'][generation-1])
# a_mean = np.mean(est_data['a'][generation-1])
# nv_mean = np.mean(est_data['nu'][generation-1])
# mean of the top 5% samples

gamma_vec_TVP = est_data['gamma_data_TVP'][-2,previous_bestfitted_index_TVP]*1e8
a_vec_TVP = est_data['a_data_TVP'][-2,previous_bestfitted_index_TVP]*1e6
nu_vec_TVP = est_data['nu_data_TVP'][-2,previous_bestfitted_index_TVP]*1e4
vm_vec_TVP = est_data['vm_data_TVP'][-2,previous_bestfitted_index_TVP]
theta_vec_TVP = est_data['theta_data_TVP'][-2,previous_bestfitted_index_TVP]


gamma_vec_TV = est_data['gamma_data_TV'][-2,previous_bestfitted_index_TV]*1e8
a_vec_TV = est_data['a_data_TV'][-2,previous_bestfitted_index_TV]*1e6
nu_vec_TV = est_data['nu_data_TV'][-2,previous_bestfitted_index_TV]*1e4
vm_vec_TV = est_data['vm_data_TV'][-2,previous_bestfitted_index_TV]
theta_vec_TV = est_data['theta_data_TV'][-2,previous_bestfitted_index_TV]


gamma_vec_TVM = est_data['gamma_data_TVM'][-2,previous_bestfitted_index_TVM]*1e8
a_vec_TVM = est_data['a_data_TVM'][-2,previous_bestfitted_index_TVM]*1e6
nu_vec_TVM = est_data['nu_data_TVM'][-2,previous_bestfitted_index_TVM]*1e4
vm_vec_TVM = est_data['vm_data_TVM'][-2,previous_bestfitted_index_TVM]
theta_vec_TVM = est_data['theta_data_TVM'][-2,previous_bestfitted_index_TVM]


timegap = 20

replicates = 1
param_TVP = DVParamLiang(gamma=1, a=1, K=1e6, h=1, nu=1, r=1, theta=1,
                         V00=.5,V01=.5, Vmax=1, inittrait=1300, initpop=1e5,
                     initpop_sigma=10.0, break_on_mu=False)
# params = np.tile(param_TVM,(replicates,1))
param_TVP[0]=np.mean(gamma_vec_TVP)*1e-8
param_TVP[1]=np.mean(a_vec_TVP)*1e-6
param_TVP[4]=np.mean(nu_vec_TVP)*1e-4
param_TVP[6]=np.mean(theta_vec_TVP)
param_TVP[9] = np.mean(vm_vec_TVP)



# trait evolution plot

count = 0

scalar_vec = [10000,20000,40000,60000,80000]

f1, axes1 = plt.subplots(len(scalar_vec), 1, figsize=(9, 9)) #
axes_flatten = axes1.flatten()
for scalar in scalar_vec:
    print(count)
    td = DVTreeData(path=obs_file, scalar=scalar)
    simresult = DVSimTVP(td,param_TVP)
    if simresult['sim_time'] == td.sim_evo_time:
        pic = 0
    else:
        pic=1

    # if pic==0:
    evo_time, total_species = simresult['N'].shape
    evo_time = evo_time - 1
    trait_RI_dr = simresult['Z']
    population_RI_dr = simresult['N']
    population_RI_dr = population_RI_dr.astype(float)
    population_RI_dr[population_RI_dr==0] = np.nan
    V_dr = simresult['V']
    num_lines = total_species
    x = np.arange(evo_time/timegap)


    labels = []
    for i in range(1, num_lines + 1):
        axes1[count].plot(x, population_RI_dr[::timegap, i - 1])
    # axes[index_g, index_a].yaxis.set_major_locator(plt.NullLocator())
    axes1[count].xaxis.set_major_locator(plt.NullLocator())
    # axes[index_g, index_a].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    # axes[index_g, index_a].set_yscale('log')
    axes1[count].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    axes1[count].set_ylabel(scalar_vec[int(count)])
    axes1[count].yaxis.set_label_position("right")
    count += 1


dir_fig = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/traittree'


f1.savefig(dir_fig+'TP2w.png')
plt.close(f1)
plt.show(f1)

def show_figure(fig):

    # create a dummy figure and use its
    # manager to display "fig"

    dummy = plt.figure()
    new_manager = dummy.canvas.manager
    new_manager.canvas.figure = fig
    fig.set_canvas(new_manager.canvas)

show_figure(f1)
