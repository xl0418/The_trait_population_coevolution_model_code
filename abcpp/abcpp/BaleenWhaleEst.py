import os
import sys
sys.path.append('C:/Liang/abcpp_ms2/abcpp')

import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
from scipy.stats import norm

import csv

# gamma = 0.001
# a = 0.1
#
# argsort of 2-dimensional arrays is awkward
# returns i,j so that X[i,j] == np.sort(X)
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err



#full tree
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'


td = DVTreeData(path=files, scalar=2000)

K=10e5
nu=1e-4

prior = [1e-7, 1e-4, 1e-5, 1e-2,nu,nu,100,50]
gamma_prior_mean = prior[0]
gamma_prior_var = prior[1]
a_prior_mean = prior[2]
a_prior_var = prior[3]
nu_prior_mean = prior[4]
nu_prior_var = prior[5]
vm_prior_mean = prior[6]
vm_prior_var = prior[7]

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
obsZ = length
print('trying to estimate the parameters','...')

s = np.argsort(obsZ)
obsZ = obsZ[s]
obsZ = obsZ.astype(np.float)

meantrait = np.mean(obsZ)
# let's try to find a true simulation:
obs_param = DVParamLiang(gamma=1, a=1, K=K,h=1, nu=nu, r=1, theta=meantrait,V00=.1,V01=.1, Vmax=100, inittrait=0, initpop=500,
                                initpop_sigma=10.0, break_on_mu=False)

# pop = dvcpp.DVSim(td, obs_param)

population = 2000
generations = 30
params = np.tile(obs_param, (population, 1))  # duplicate
params[:, 0] = np.random.uniform(0.0, 1e-5, params.shape[0])  # randomize 'gamma'
params[:, 1] = np.random.uniform(0.0, 1e-3, params.shape[0])  # randomize 'a'
params[:, 4] = np.random.uniform(0.0, nu*10, params.shape[0])  # randomize 'nu'
params[:, 9] = np.random.uniform(50, 150, params.shape[0])  # randomize 'Vm'


gamma_data = np.zeros(shape=(generations, population))
a_data = np.zeros(shape=(generations, population))
nu_data = np.zeros(shape=(generations, population))
vm_data = np.zeros(shape=(generations, population))

fitness= np.zeros(shape=(generations, population))
# Initialize the weights.
weight_gamma = np.zeros(population)
weight_gamma.fill(1 / population)
weight_a = np.zeros(population)
weight_a.fill(1 / population)
weight_nu = np.zeros(population)
weight_nu.fill(1 / population)
weight_vm = np.zeros(population)
weight_vm.fill(1 / population)
for g in range(generations):
    gamma_data[g, :] = params[:, 0]
    a_data[g, :] = params[:, 1]
    nu_data[g,:] = params[:,4]
    vm_data[g,:] = params[:,9]
    pop = dvcpp.DVSimLiang(td, params)

    # access fitness
    # fitness = np.zeros(population)
    valid = np.where(pop['sim_time'] == td.sim_evo_time)[0]
    if len(valid)<20:
        print("WARNING:Valid simulations are too scarce!")
    if valid.size > 0:
        Z = pop['Z'][valid]
        i, j = argsort2D(Z)
        Z = Z[i, j]
        # V = pop['V'][valid][i, j]
        Z = np.nan_to_num(Z)
        # V = np.nan_to_num(V)
        #GOF: Goodness of fit
        fitness[g,valid] += 1.0 - normalized_norm(Z, obsZ)
        # fitness[g,valid] += 1.0 - normalized_norm(np.sqrt(V), np.sqrt(obsV))

    # print something...
    q5 = np.argsort(fitness[g,:])[-len(valid)// 20]  # best 5%
    fit_index = np.where(fitness[g,:] > fitness[g,q5])[0]
    print('Valid length: %i' % len(valid))
    print('Iteration = %d 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f fitness = %f' % (g, np.mean(params[fit_index, 0]),
                                                                 np.mean(params[fit_index, 1]),np.mean(params[fit_index, 4]),
                                                                 np.mean(params[fit_index, 9]), np.mean(fitness[g,fit_index])))


    weight_gamma = weight_gamma[fit_index]/sum(weight_gamma[fit_index])
    weight_a = weight_a[fit_index]/sum(weight_a[fit_index])
    weight_nu = weight_nu[fit_index]/sum(weight_nu[fit_index])
    weight_vm = weight_vm[fit_index]/sum(weight_vm[fit_index])

    gamma_pre_mean = np.sum(params[fit_index, 0] * weight_gamma)
    gamma_pre_var = np.sum((params[fit_index, 0] - gamma_pre_mean) ** 2 * weight_gamma)
    a_pre_mean = np.sum(params[fit_index, 1] * weight_a)
    a_pre_var = np.sum((params[fit_index, 1] - a_pre_mean) ** 2 * weight_a)
    nu_pre_mean = np.sum(params[fit_index, 4] * weight_nu)
    nu_pre_var = np.sum((params[fit_index, 4] - nu_pre_mean) ** 2 * weight_nu)
    vm_pre_mean = np.sum(params[fit_index, 9] * weight_vm)
    vm_pre_var = np.sum((params[fit_index, 9] - vm_pre_mean) ** 2 * weight_vm)

    # sample parameters by the weights computed in last loop.
    sample_gamma_index = np.random.choice(fit_index, population, p=weight_gamma)
    sample_a_index = np.random.choice(fit_index, population, p=weight_a)
    sample_nu_index = np.random.choice(fit_index, population, p=weight_nu)
    sample_vm_index = np.random.choice(fit_index, population, p=weight_vm)

    # mean of the sample for gamma
    propose_gamma0 = params[sample_gamma_index, 0]
    # draw new gamma with mean and variance
    propose_gamma = abs(np.random.normal(propose_gamma0, np.sqrt(2 * gamma_pre_var)))
    # mean of the sample for a
    propose_a0 = params[sample_a_index, 1]
    # draw new a with mean and variance
    propose_a = abs(np.random.normal(propose_a0, np.sqrt(2 * a_pre_var)))
    # mean of the sample for nu
    propose_nu0 = params[sample_nu_index, 4]
    # draw new nu with mean and variance
    propose_nu = abs(np.random.normal(propose_nu0, np.sqrt(2 * nu_pre_var)))
    # mean of the sample for vm
    propose_vm0 = params[sample_vm_index, 9]
    # draw new nu with mean and variance
    propose_vm = abs(np.random.normal(propose_vm0, np.sqrt(2 * vm_pre_var)))

    extend_weight_gamma = weight_gamma[fit_index.searchsorted(sample_gamma_index)]
    extend_weight_a = weight_a[fit_index.searchsorted(sample_a_index)]
    extend_weight_nu = weight_nu[fit_index.searchsorted(sample_nu_index)]
    extend_weight_vm = weight_vm[fit_index.searchsorted(sample_vm_index)]

    # compute new weights for gamma and a
    weight_gamma_denominator = np.sum(extend_weight_gamma * norm.pdf(propose_gamma, params[:, 0],
                                                                     np.sqrt(2 * gamma_pre_var)))
    weight_gamma_numerator = norm.pdf(propose_gamma, gamma_prior_mean, gamma_prior_var)
    weight_gamma = weight_gamma_numerator / weight_gamma_denominator

    weight_a_denominator = np.sum(extend_weight_a * norm.pdf(propose_a, params[:, 1],
                                                             np.sqrt(2 * a_pre_var)))
    weight_a_numerator = norm.pdf(propose_a, a_prior_mean, a_prior_var)
    weight_a = weight_a_numerator / weight_a_denominator

    weight_nu_denominator = np.sum(extend_weight_nu * norm.pdf(propose_nu, params[:, 4],
                                                             np.sqrt(2 * nu_pre_var)))
    weight_nu_numerator = norm.pdf(propose_nu, nu_prior_mean, nu_prior_var)
    weight_nu = weight_nu_numerator / weight_nu_denominator

    weight_vm_denominator = np.sum(extend_weight_vm * norm.pdf(propose_vm, params[:, 9],
                                                             np.sqrt(2 * vm_pre_var)))
    weight_vm_numerator = norm.pdf(propose_vm, vm_prior_mean, vm_prior_var)
    weight_vm = weight_vm_numerator / weight_vm_denominator
    # normalize the weights
    # total_simulation[t] = sim_count
    weight_gamma = weight_gamma / sum(weight_gamma)
    weight_a = weight_a / sum(weight_a)
    weight_nu = weight_nu / sum(weight_nu)
    weight_vm = weight_vm / sum(weight_vm)

    params[:, 0] = propose_gamma
    params[:, 1] = propose_a
    params[:, 4] = propose_nu
    params[:, 9] = propose_vm

#
para_data = {'gamma': gamma_data, 'a': a_data, 'nu': nu_data,'vm': vm_data,'fitness': fitness}
# file='C:/Liang/Code/Pro2/abcpp/abcpp/smcndata/'
# # file = '/home/p274981/abcpp/abcpp/'
# file2 = file + 'paraT11.npy'
# np.save(file2,para_data)
#
# smc = np.load(file2).item()