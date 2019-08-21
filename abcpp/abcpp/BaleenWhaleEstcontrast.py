import sys
sys.path.append('C:/Liang/abcpp_ms5/abcpp')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
from scipy.stats import norm
import csv
import dendropy
from dendropy.model import continuous
import time
from multiprocessing import Pool
from itertools import repeat
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from pic_compute import pic_compute

# gamma = 0.001
# a = 0.1
#
# argsort of 2-dimensional arrays is awkward
# returns i,j so that X[i,j] == np.sort(X)
num_cores = Pool(8)  # the number of cores

def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err
    # calculate pic of the simulated traits

#full tree
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'
    #'/home/p274981/abcpp/'

files = dir_path + 'treedata/'



time_scalar = 20000
heri_sqr = 1

td = DVTreeData(path=files, scalar=time_scalar)


K=10e5
nu=1e-3

file_result = dir_path + 'BaleenWhales/BWest_t%i_h%f.npy' % (int(time_scalar),heri_sqr)

with open(files+'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(files+'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:,1]
lengthdata_array = np.array(lengthdata)
length_index = []
for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)
length = 10**logTL


# reorder the traits according to the simulation order
sim_species_label = [ "Balaena_mysticetus","Balaenoptera_acutorostrata" ,"Caperea_marginata" ,    "Balaenoptera_borealis" ,
  "Balaenoptera_physalus" ,     "Eschrichtius_robustus"  ,    "Balaenoptera_musculus"    ,  "Balaenoptera_omurai"    ,
 "Eubalaena_australis" ,    "Megaptera_novaeangliae" , "Balaenoptera_bonaerensis"  , "Balaenoptera_brydei"   ,
 "Balaenoptera_edeni"     ,    "Eubalaena_glacialis"   ,  "Eubalaena_japonica" ]

obsZ_ordered_sim = length[[np.where(sim_species_label[i] == extantlabels_array)[0][0] for i in range(15)]]
obsZ = obsZ_ordered_sim

# PIC calculation
taxa1 = dendropy.TaxonNamespace()
dataset_combined = dendropy.DataSet.get(path=files+"bw_char.nex",schema="nexus")
tree_emp = dataset_combined.tree_lists[0][0]
chars_emp = dataset_combined.char_matrices[0]
pic_emp = continuous.PhylogeneticIndependentConstrasts(tree=tree_emp,
char_matrix=chars_emp)
ctree_emp = pic_emp.contrasts_tree(character_index=0,
annotate_pic_statistics=True,state_values_as_node_labels=False,corrected_edge_lengths=False)
emp_pic = []
label = []
for nd in ctree_emp.postorder_internal_node_iter():
    emp_pic.append(nd.pic_contrast_standardized)
    label.append(int(nd.label))

emp_pic_orded_node = abs(np.array(emp_pic)[np.argsort(label)])


tree_sim = dendropy.Tree.get(
    path=files + "bw.nex", schema="nexus",
    taxon_namespace=taxa1)

# compute pic


print('trying to estimate the parameters','...')



prior = [1e-7, 1e-4, 1e-5, 1e-2,nu,nu,100,100,1000,1000]
gamma_prior_mean = prior[0]
gamma_prior_var = prior[1]
a_prior_mean = prior[2]
a_prior_var = prior[3]
nu_prior_mean = prior[4]
nu_prior_var = prior[5]
vm_prior_mean = prior[6]
vm_prior_var = prior[7]
theta_prior_mean = prior[8]
theta_prior_var = prior[9]


meantrait = np.mean(obsZ)
# let's try to find a true simulation:
obs_param = DVParamLiang(gamma=1, a=1, K=K,h=np.sqrt(heri_sqr), nu=nu, r=1, theta=1,V00=.1,V01=.1, Vmax=100, inittrait=meantrait, initpop=500,
                                initpop_sigma=10.0, break_on_mu=False)

# pop = dvcpp.DVSim(td, obs_param)

population = 30000
generations = 20
params = np.tile(obs_param, (population, 1))  # duplicate
params[:, 0] = np.random.uniform(0.0, 1e-5, params.shape[0])  # randomize 'gamma'
params[:, 1] = np.random.uniform(0.0, 1e-2, params.shape[0])  # randomize 'a'
params[:, 4] = np.random.uniform(0.0, nu*100, params.shape[0])  # randomize 'nu'
params[:, 9] = np.random.uniform(50, 150, params.shape[0])  # randomize 'Vm'
params[:, 6] = np.random.uniform(np.min(obsZ),np.max(obsZ), params.shape[0])  # randomize 'Vm'


gamma_data = np.zeros(shape=(generations, population))
a_data = np.zeros(shape=(generations, population))
nu_data = np.zeros(shape=(generations, population))
vm_data = np.zeros(shape=(generations, population))
theta_data = np.zeros(shape=(generations, population))

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
weight_theta = np.zeros(population)
weight_theta.fill(1 / population)
for g in range(generations):

    gamma_data[g, :] = params[:, 0]
    a_data[g, :] = params[:, 1]
    nu_data[g,:] = params[:,4]
    vm_data[g,:] = params[:,9]
    theta_data[g,:] = params[:,6]
    start_time = time.time()
    pop = dvcpp.DVSimTVP(td, params)
    print("---simulation: %s seconds ---" % (time.time() - start_time))

    # access fitness
    # fitness = np.zeros(population)
    valid = np.where(pop['sim_time'] == td.sim_evo_time)[0]
    num_valid_sims = len(valid)
    if num_valid_sims<20:
        print("WARNING:Valid simulations are too scarce!")
    if num_valid_sims > 0:
        Z = pop['Z'][valid]
        # V = pop['V'][valid][i, j]
        Z = np.nan_to_num(Z)
        # V = np.nan_to_num(V)
        start_time = time.time()

        # # calculate pic of the simulated traits
        # tree_sim = dendropy.Tree.get(
        #     path=files + "bw.nex", schema="nexus",
        #     taxon_namespace=taxa1)
        #
        # simchar_dict = {}
        # keys = ["B.mysticetus", "B.acutorostrata", "C.marginata","B.borealis",
        #             "B.physalus", "E.robustus", "B.musculus", "B.omurai",
        #             "E.australis", "M.novaeangliae", "B.bonaerensis", "B.brydei",
        #             "B.edeni", "E.glacialis", "E.japonica"]
        # values = Z
        # for i in range(15):
        #     simchar_dict[keys[i]]=values[:,i].tolist()
        # simchars = dendropy.ContinuousCharacterMatrix.from_dict(simchar_dict, taxon_namespace=taxa1)
        # simpic = continuous.PhylogeneticIndependentConstrasts(tree=tree_sim, char_matrix=simchars)
        # sim_pic_thisbatch = []
        # for pic_each in range(num_valid_sims):
        #     sim_ctree = simpic.contrasts_tree(character_index=pic_each,
        #                                       annotate_pic_statistics=True,
        #                                       state_values_as_node_labels=False,
        #                                       corrected_edge_lengths=False)
        #     sim_pic = []
        #     sim_label = []
        #     for nd in sim_ctree.postorder_internal_node_iter():
        #         sim_pic.append(nd.pic_contrast_raw)
        #         sim_label.append(int(nd.label))
        #     sim_pic_ordered = np.array(sim_pic)[np.argsort(sim_label)]
        #     sim_pic_thisbatch.append(sim_pic_ordered)


        pic_ordered_list = num_cores.starmap(pic_compute, zip(repeat(tree_sim), Z,repeat(taxa1),range(num_valid_sims)))

        # in case the parallel computation returns disordered output
        order_list = []
        contrast_list = []
        for i in range(num_valid_sims):
            order_list.append(pic_ordered_list[i][1])
            contrast_list.append(pic_ordered_list[i][0])
        ordered_contrast_list = [contrast_list[item] for item in np.argsort(order_list)]
        contrast_array = abs(np.vstack(ordered_contrast_list))
        # test_list = []
        # for i in range(len(valid)):
        #     test_list.append(pic_ordered_list[i][0])

        print("---PIC calculation: %s seconds ---" % (time.time() - start_time))

        #GOF: Goodness of fit
        # fitness[g,valid] += 1.0 - normalized_norm(Z, obsZ)
        fitness[g,valid] += 1.0 - normalized_norm(emp_pic_orded_node, contrast_array)

    # print something...
    q5 = np.argsort(fitness[g,:])[-len(valid)// 200]  # best 0.5%
    fit_index = np.where(fitness[g,:] > fitness[g,q5])[0]

    print('Iteration = %d 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f fitness = %f' % (g, np.mean(params[fit_index, 0]),
                 np.mean(params[fit_index, 1]),np.mean(params[fit_index, 4]),
                 np.mean(params[fit_index, 9]),np.mean(params[fit_index, 6]), np.mean(fitness[g,fit_index])))


    weight_gamma = weight_gamma[fit_index]/sum(weight_gamma[fit_index])
    weight_a = weight_a[fit_index]/sum(weight_a[fit_index])
    weight_nu = weight_nu[fit_index]/sum(weight_nu[fit_index])
    weight_vm = weight_vm[fit_index]/sum(weight_vm[fit_index])
    weight_theta = weight_theta[fit_index] / sum(weight_theta[fit_index])

    gamma_pre_mean = np.sum(params[fit_index, 0] * weight_gamma)
    gamma_pre_var = np.sum((params[fit_index, 0] - gamma_pre_mean) ** 2 * weight_gamma)
    a_pre_mean = np.sum(params[fit_index, 1] * weight_a)
    a_pre_var = np.sum((params[fit_index, 1] - a_pre_mean) ** 2 * weight_a)
    nu_pre_mean = np.sum(params[fit_index, 4] * weight_nu)
    nu_pre_var = np.sum((params[fit_index, 4] - nu_pre_mean) ** 2 * weight_nu)
    vm_pre_mean = np.sum(params[fit_index, 9] * weight_vm)
    vm_pre_var = np.sum((params[fit_index, 9] - vm_pre_mean) ** 2 * weight_vm)
    theta_pre_mean = np.sum(params[fit_index, 6] * weight_theta)
    theta_pre_var = np.sum((params[fit_index, 6] - theta_pre_mean) ** 2 * weight_theta)

    # sample parameters by the weights computed in last loop.
    sample_gamma_index = np.random.choice(fit_index, population, p=weight_gamma)
    sample_a_index = np.random.choice(fit_index, population, p=weight_a)
    sample_nu_index = np.random.choice(fit_index, population, p=weight_nu)
    sample_vm_index = np.random.choice(fit_index, population, p=weight_vm)
    sample_theta_index = np.random.choice(fit_index, population, p=weight_theta)

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
    # draw new vm with mean and variance
    propose_vm = abs(np.random.normal(propose_vm0, np.sqrt(2 * vm_pre_var)))
    # mean of the sample for theta
    propose_theta0 = params[sample_theta_index, 6]
    # draw new theta with mean and variance
    propose_theta = abs(np.random.normal(propose_theta0, np.sqrt(2 * theta_pre_var)))

    extend_weight_gamma = weight_gamma[fit_index.searchsorted(sample_gamma_index)]
    extend_weight_a = weight_a[fit_index.searchsorted(sample_a_index)]
    extend_weight_nu = weight_nu[fit_index.searchsorted(sample_nu_index)]
    extend_weight_vm = weight_vm[fit_index.searchsorted(sample_vm_index)]
    extend_weight_theta = weight_theta[fit_index.searchsorted(sample_theta_index)]

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

    weight_theta_denominator = np.sum(extend_weight_theta * norm.pdf(propose_theta, params[:, 6],
                                                             np.sqrt(2 * theta_pre_var)))
    weight_theta_numerator = norm.pdf(propose_theta, theta_prior_mean, theta_prior_var)
    weight_theta = weight_theta_numerator / weight_theta_denominator
    # normalize the weights
    # total_simulation[t] = sim_count
    weight_gamma = weight_gamma / sum(weight_gamma)
    weight_a = weight_a / sum(weight_a)
    weight_nu = weight_nu / sum(weight_nu)
    weight_vm = weight_vm / sum(weight_vm)
    weight_theta = weight_theta / sum(weight_theta)

    params[:, 0] = propose_gamma
    params[:, 1] = propose_a
    params[:, 4] = propose_nu
    params[:, 9] = propose_vm
    params[:, 6] = propose_theta

#
para_data = {'gamma': gamma_data, 'a': a_data, 'nu': nu_data,'vm': vm_data,'theta': theta_data,
                'fitness': fitness,'Z':Z,'valid':valid}


# np.save(file_result,para_data)
