import sys, os
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_ms2/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_py import DVSimLiang
from dvtraitsim_shared import DVTreeData, DVParamLiang


theta = 0  # optimum of natural selection
r = 1  # growth rate
Vmax = 1
scalar = 20000
K=10e5
nu=1e-4
# let's try to find a true simulation:


#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'

initrait = 0
td = DVTreeData(path=files, scalar=scalar)

# gamma, a, K, h, nu, r, theta, V00, V01, Vmax, inittrait, initpop, initpop_sigma, break_on_mu
obs_param = DVParamLiang(gamma=0.0000001, a=0.0001, K=K,h=1, nu=nu, r=r, theta=theta,V00=.1,V01=.1, Vmax=200, inittrait=initrait, initpop=500,
                    initpop_sigma=10.0, break_on_mu=False)
simresult = DVSimLiang(td,obs_param)

row = simresult['Z'].shape[0]-1
np.max(simresult['Z'][row,:])-np.min(simresult['Z'][row,:])