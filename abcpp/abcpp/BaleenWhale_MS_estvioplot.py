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
data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster' \
           '/results_0702/'

obs_file = dir_path + 'treedata/'

with open(obs_file + 'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(dir_path + 'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:, 1]
lengthdata_array = np.array(lengthdata)
length_index = []
gamma_list = []
a_list = []
nu_list = []
vm_list = []
distance_list = []
ratio_dis = []
valid_length = []
timescaling_list = []
heri_list = []
count = 0
particle_size = 100
K = 1e6
nu = 1 / (100 * K)

for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:, 0] == label)[0][0])

logTL = lengthdata_array[length_index, 1].astype(np.float)
obsZ = sorted(10 ** logTL)

# for timescaling_index in range(len(timescale_vec)):
#     for heritability_index in range(len(heritability_vec)):
# print(count)
# count += 1
# timescaling = timescale_vec[timescaling_index]
# heritability = heritability_vec[heritability_index]

data_name = data_dir + 'modelselec2w.npy'
est_data = np.load(data_name).item()
population = int(len(est_data['model_data'][0]) / 3)
fitness = est_data['fitness'][-1]

q5_TVP = np.argsort(fitness[:population])[-int(population // 20 + 1)]  # best 5%
q5_TV = np.argsort(fitness[population:2 * population])[
            -int(population // 20 + 1)] + population  # best 5%
q5_TVM = np.argsort(fitness[2 * population:])[-int(population // 20 + 1)] + 2 * population  #
# best 5%

fit_index_TVP = np.where(fitness[:population] > fitness[q5_TVP])[0]
fit_index_TV = np.where(fitness[population:2 * population] > fitness[q5_TV])[0] + population
fit_index_TVM = np.where(fitness[2 * population:] > fitness[q5_TVM])[0] + 2 * population

previous_bestfitted_index_TVP = fit_index_TVP
previous_bestfitted_index_TV = fit_index_TV - population
previous_bestfitted_index_TVM = fit_index_TVM - 2 * population

top_pop = len(previous_bestfitted_index_TVP)

# last_fitness = est_data['fitness'][generation - 1]
# q5 = np.argsort(last_fitness)[-population // 20]  # best 5%
# fit_index = np.where(last_fitness > last_fitness[q5])[0]

# mean of all samples
# gamma_mean = np.mean(est_data['gamma'][generation-1])
# a_mean = np.mean(est_data['a'][generation-1])
# nv_mean = np.mean(est_data['nu'][generation-1])
# mean of the top 5% samples

gamma_vec_TVP = est_data['gamma_data_TVP'][-2, previous_bestfitted_index_TVP] * 1e8
a_vec_TVP = est_data['a_data_TVP'][-2, previous_bestfitted_index_TVP] * 1e6
nu_vec_TVP = est_data['nu_data_TVP'][-2, previous_bestfitted_index_TVP] * 1e4
vm_vec_TVP = est_data['vm_data_TVP'][-2, previous_bestfitted_index_TVP]
theta_vec_TVP = est_data['theta_data_TVP'][-2, previous_bestfitted_index_TVP]

gamma_vec_TV = est_data['gamma_data_TV'][-2, previous_bestfitted_index_TV] * 1e8
a_vec_TV = est_data['a_data_TV'][-2, previous_bestfitted_index_TV] * 1e6
nu_vec_TV = est_data['nu_data_TV'][-2, previous_bestfitted_index_TV] * 1e4
vm_vec_TV = est_data['vm_data_TV'][-2, previous_bestfitted_index_TV]
theta_vec_TV = est_data['theta_data_TV'][-2, previous_bestfitted_index_TV]

gamma_vec_TVM = est_data['gamma_data_TVM'][-2, previous_bestfitted_index_TVM] * 1e8
a_vec_TVM = est_data['a_data_TVM'][-2, previous_bestfitted_index_TVM] * 1e6
nu_vec_TVM = est_data['nu_data_TVM'][-2, previous_bestfitted_index_TVM] * 1e4
vm_vec_TVM = est_data['vm_data_TVM'][-2, previous_bestfitted_index_TVM]
theta_vec_TVM = est_data['theta_data_TVM'][-2, previous_bestfitted_index_TVM]

print('TVP 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.mean(gamma_vec_TVP) * 1e-8,
    np.mean(a_vec_TVP) * 1e-6, np.mean(nu_vec_TVP) * 1e-4,
    np.mean(vm_vec_TVP), np.mean(theta_vec_TVP)))
print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.var(gamma_vec_TVP) * 1e-8,
    np.var(a_vec_TVP) * 1e-6, np.var(nu_vec_TVP) * 1e-4,
    np.var(vm_vec_TVP), np.var(theta_vec_TVP)))

print('TV 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.mean(gamma_vec_TV) * 1e-8,
    np.mean(a_vec_TV) * 1e-6, np.mean(nu_vec_TV) * 1e-4,
    np.mean(vm_vec_TV), np.mean(theta_vec_TV)))
print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.var(gamma_vec_TV) * 1e-8,
    np.var(a_vec_TV) * 1e-6, np.var(nu_vec_TV) * 1e-4,
    np.var(vm_vec_TV), np.var(theta_vec_TV)))

print('TVM 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.mean(gamma_vec_TVM) * 1e-8,
    np.mean(a_vec_TVM) * 1e-6, np.mean(nu_vec_TVM) * 1e-4,
    np.mean(vm_vec_TVM), np.mean(theta_vec_TVM)))
print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.var(gamma_vec_TVM) * 1e-8,
    np.var(a_vec_TVM) * 1e-6, np.var(nu_vec_TVM) * 1e-4,
    np.var(vm_vec_TVM), np.var(theta_vec_TVM)))

est_para = ['gamma', 'alpha', 'nu', 'vm', 'theta']  # ,'vm'
model_para = ['TVP', 'TV', 'TVM']
est_array = np.concatenate([gamma_vec_TVP, a_vec_TVP, nu_vec_TVP, vm_vec_TVP, theta_vec_TVP,
                            gamma_vec_TV, a_vec_TV, nu_vec_TV, vm_vec_TV, theta_vec_TV,
                            gamma_vec_TVM, a_vec_TVM, nu_vec_TVM, vm_vec_TVM,
                            theta_vec_TVM])  # ,vm_list_flat
est_label = np.tile(np.repeat(est_para, top_pop), len(model_para))
model_label = np.repeat(model_para, top_pop * len(est_para))

ss_list = {'est': est_array, 'est_label': est_label,
           'model_label': model_label}
ss_df = pd.DataFrame(ss_list)

vioplot = sns.catplot(x="model_label", y="est", col="est_label",
                      data=ss_df, kind="box", height=5, aspect=.6, sharey=False)
# vioplot.set(ylim=(-50,400))
# vioplot.set_xticklabels(["$\gamma \cdot 10^{-8}$", "$\\alpha \cdot 10^{-5}$", "$\\nu \cdot 10^{
# -4}$","$V_{m}$"])
vioplot.set_axis_labels("", "Estimate value")
# vioplot._legend.set_title('$h^2$')
axes = vioplot.axes.flatten()
axes[0].set_title("$\gamma \cdot 10^{-8}$")
axes[1].set_title("$\\alpha \cdot 10^{-6}$")
axes[2].set_title("$\\nu \cdot 10^{-4}$")
axes[3].set_title("$V_m $")
axes[4].set_title("$\\theta $")

Z = est_data['Z']
diff_norm = np.linalg.norm(Z - obsZ, axis=1)
plt.hist(diff_norm, bins=200)
TVP_index = np.where(fitness[:population] > 0)[0]
TV_index = np.where(fitness[population:2 * population] > 0)[0] + population
TVM_index = np.where(fitness[2 * population:] > 0)[0] + 2 * population
plot_tvm_index = np.where(fitness[TVM_index] > np.min(fitness[TVP_index]))[0] + 2 * population

plt.hist(diff_norm[:len(TVP_index)])
plt.hist(diff_norm[len(TVP_index):len(TVP_index) + len(TV_index)])
plt.hist(diff_norm[plot_tvm_index])

import sys

sys.path.append('C:/Liang/abcpp_ms5/abcpp')
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp

td = DVTreeData(path=obs_file, scalar=1000)

largeini = 100
param_TVM = DVParamLiang(gamma=1, a=1, K=1e6, h=1, nu=1, r=1, theta=1,
                         V00=.5, V01=.5, Vmax=1, inittrait=1300, initpop=1e5,
                         initpop_sigma=10.0, break_on_mu=False)
params = np.tile(param_TVM, (largeini, 1))
params[:, 0] = np.mean(gamma_vec_TVP) * 1e-8
params[:, 1] = np.mean(a_vec_TVP) * 1e-6
params[:, 4] = np.mean(nu_vec_TVP) * 1e-4
params[:, 6] = np.mean(theta_vec_TVP)

vm_test = [20, 60, 100, 140, 180, 220]
diff_list = []
vm_list = []
v_list = []
n_list = []
for vm in vm_test:
    params[:, 9] = vm
    predictsim = dvcpp.DVSimTVP(td, params)
    valid = np.where(predictsim['sim_time'] == td.sim_evo_time)[0]
    Z = predictsim['Z'][valid]
    V = predictsim['V'][valid]
    N = predictsim['N'][valid]
    i, j = argsort2D(Z)
    Z_modelTVP = Z[i, j]
    N_modelTVP = N[i, j]
    V_modelTVP = V[i, j]
    diff_norm = np.linalg.norm(Z_modelTVP - obsZ, axis=1)
    diff_list.append(diff_norm)
    vm_list.append(np.repeat(vm, len(diff_norm)))
    v_list.append(V_modelTVP.flatten())
    n_list.append(N_modelTVP.flatten())

diff_flatten = np.array([item for sublist in diff_list for item in sublist])
vm_flatten = np.array([item for sublist in vm_list for item in sublist])
v_flatten = np.array([item for sublist in v_list for item in sublist])
n_flatten = np.array([item for sublist in n_list for item in sublist])

vm_v_label = np.repeat(vm_test, 15 * largeini)
species_label = np.tile(range(15), largeini * len(vm_test))
diff_df = pd.DataFrame({'diff': diff_flatten, 'vm': vm_flatten})
v_df = pd.DataFrame({'v': v_flatten, 'vm': vm_v_label})
n_df = pd.DataFrame({'n': n_flatten, 'species': species_label, 'vm': vm_v_label})

# goodness of fit
for vm in vm_test:
    sns.distplot(diff_df[diff_df["vm"] == vm]['diff'], label=str(vm))
plt.legend()

# variance
for vm in vm_test:
    sns.distplot(v_df[v_df["vm"] == vm]['v'], label=vm)
plt.legend()

# population

g = sns.catplot(x="species", y="n",
                hue="vm", col="vm",
                data=n_df, kind="box",
                height=4, aspect=.7)
for ax in g.axes.flat:
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

# n distribution across different time steps

replicates = 100
param_TVM = DVParamLiang(gamma=1, a=1, K=1e6, h=1, nu=1, r=1, theta=1,
                         V00=.5, V01=.5, Vmax=1, inittrait=1300, initpop=1e5,
                         initpop_sigma=10.0, break_on_mu=False)
params = np.tile(param_TVM, (replicates, 1))
params[:, 0] = np.mean(gamma_vec_TVP) * 1e-8
params[:, 1] = np.mean(a_vec_TVP) * 1e-6
params[:, 4] = np.mean(nu_vec_TVP) * 1e-4
params[:, 6] = np.mean(theta_vec_TVP)
params[:, 9] = np.mean(vm_vec_TVP)

t_test = [5000, 10000, 20000, 40000, 80000, 120000]  #
diff_list = []
time_list = []
v_list = []
n_list = []
for timevalue in t_test:
    td = DVTreeData(path=obs_file, scalar=timevalue)
    predictsim = dvcpp.DVSimTVP(td, params)
    valid = np.where(predictsim['sim_time'] == td.sim_evo_time)[0]
    Z = predictsim['Z'][valid]
    V = predictsim['V'][valid]
    N = predictsim['N'][valid]
    i, j = argsort2D(Z)
    Z_modelTVP = Z[i, j]
    N_modelTVP = N[i, j]
    V_modelTVP = V[i, j]
    diff_norm = np.linalg.norm(Z_modelTVP - obsZ, axis=1)
    diff_list.append(diff_norm)
    time_list.append(np.repeat(timevalue, len(diff_norm)))
    v_list.append(V_modelTVP.flatten())
    n_list.append(N_modelTVP.flatten())

diff_flatten = np.array([item for sublist in diff_list for item in sublist])
time_flatten = np.array([item for sublist in time_list for item in sublist])
v_flatten = np.array([item for sublist in v_list for item in sublist])
n_flatten = np.array([item for sublist in n_list for item in sublist])

time_v_label = np.repeat(t_test, 15 * replicates)
species_label = np.tile(range(15), replicates * len(t_test))
diff_df = pd.DataFrame({'diff': diff_flatten, 'time': time_flatten})
v_df = pd.DataFrame({'v': v_flatten, 'time': time_v_label})
n_df = pd.DataFrame({'n': n_flatten, 'species': species_label, 'time': time_v_label})

# goodness of fit
for time in t_test:
    sns.distplot(diff_df[diff_df["time"] == time]['diff'], label=str(time))
plt.legend()

# variance
for time in t_test:
    sns.distplot(v_df[v_df["time"] == time]['v'], label=time)
plt.legend()

# population

g = sns.catplot(x="species", y="n",
                hue="time", col="time",
                data=n_df, kind="box",
                height=4, aspect=.7)
for ax in g.axes.flat:
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))


# calculate growth rate
def competition_functions_Liang(a, zi, nj):
    """ competition functions, Liang's model.

    returns beta = Sum_j( exp(-a(zi-zj)^2) * Nj)
            sigma = Sum_j( 2a * (zi-zj) * exp(-a(zi-zj)^2) * Nj)
            sigmaSqr = Sum_j( 4a^2 * (zi-zj)^2 * exp(-a(zi-zj)^2) * Nj)
    """
    T = zi[:, np.newaxis] - zi  # trait-distance matrix (via 'broadcasting')
    t1 = np.exp(-a * T ** 2) * nj
    t2 = (2 * a) * T
    beta = np.sum(t1, axis=1)
    # sigma = np.sum(t2 * t1, axis=1)
    # sigmasqr = np.sum(t2 ** 2 * t1, axis=1)
    return beta  # , sigma, sigmasqr


gamma = np.mean(gamma_vec_TVP) * 1e-8
dtz = np.mean(theta_vec_TVP) - Z_modelTVP[0]
beta = competition_functions_Liang(np.mean(a_vec_TVP) * 1e-6, Z_modelTVP[0], N_modelTVP[0])
gr1q = np.exp(-gamma * dtz ** 2 + (1 - beta / 1e6))
