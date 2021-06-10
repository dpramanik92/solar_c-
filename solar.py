import numpy as np
import solchi2 as sol
import pyswarms as ps

from tqdm import *
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

import scipy.interpolate as interp
from scipy.optimize import minimize as mini

def prior(val,bf,sigma):
    return ((val-bf)/sigma)**2.0


class data_set:
    def __init__(self,_data,_level):
        self.data = _data
        self.level = _level

        self.prepare_data()


    def prepare_data(self):
        x = self.data[0]
        y = self.data[1]
        z = self.data[2]

        N = 1000
        _x = np.log10(x)

        chi2_min = np.min(z)

        z = z - chi2_min

        xi = np.linspace(_x.min(),_x.max(),N)
        yi = np.linspace(y.min(),y.max(),N)

        zi = interp.griddata((_x,y),z,(xi[None,:],yi[:,None]))

        cs = plt.contour(10**xi,yi,zi,[self.level])
        p = cs.collections[0].get_paths()[0]

        v = p.vertices
        self.X = v[:,0]
        self.Y = v[:,1]

        return 0

class Optimize:
    def __init__(self,exp,which,chan):
        self.s = sol.solchi2()
        self.s.init_exp(exp,which,chan)

    def take_cord(self,x,y):
        self.X = x
        self.Y = y

    def locFunction(self,params):
        _params = [params[0]*(np.pi/180.0),params[1],0.0,self.Y,params[2]*1e-5,
                                2.5e-3,self.X,self.X]
        _params = np.array(_params)


        val = self.s.calChi2sing(_params)+prior(params[0],33.44,0.78)+prior(params[1],
                                0.02221,0.0006)+prior(params[2],7.42,0.21)

        return val

    def function(self,params):

        _params = []
        for i in range(len(params)):
            _param = [params[i][0]*(np.pi/180.0),params[i][1],0.0,self.Y,
                   params[i][2]*1e-5,2.5e-3,self.X,self.X]
            _param = np.array(_param)

            _params.append(_param)

        _params = np.array(_params)

        val = []
        val = self.s.calChi2(_params.T)+prior(params.T[0],33.44,0.78)+prior(params.T[1],
                                0.02221,0.0006)+prior(params.T[2],7.42,0.21)
        f = np.transpose(val)

        return f


def get_curves(exp,which,chan):

    Opt = Optimize(exp,which,chan)

    min_bound = np.array([31.0,0.02,6.8])
    max_bound = np.array([36.5,0.024,8.2])

    initial_guess = np.array([33.44,0.022,7.42])

    min_b = (31.0,0.02,6.8)
    max_b = (36.5,0.024,8.2)

    bound_ = ((0.0,90.0),(0.0,1),(0,15))

    bounds = (min_bound,max_bound)

    options = {'c1':0.5,'c2':0.3,'w':0.9}


    tau = np.arange(-6.0,-1.0,0.2)
    delt = np.arange(0.01,0.99,0.02)

    x = []
    y = []
    chi = []

    for t in tqdm(tau):
        for d in delt:
          #  optimizer = ps.single.GlobalBestPSO(n_particles=200,dimensions=3,options=options
          #                                    ,bounds=bounds,ftol=0.0001)
            Opt.take_cord(10**t,d)

         #   best_cost,best_ops = optimizer.optimize(Opt.function,iters=1000,verbose=False)

            res = mini(Opt.locFunction,initial_guess,tol=1e-3,method='Powell',bounds=bound_)
            x.append(10**t)
            y.append(d)
            chi.append(res.fun)

    x = np.array(x)
    y = np.array(y)
    chi = np.array(chi)

    data = np.array([x,y,chi])

    return data


kam_data_sc4 = get_curves("Kamland","Scalar",4)
bor_data_sc4 = get_curves("Borexino","Scalar",4)

np.savetxt('Kamland_pseudo_margi_4',kam_data_sc4)
np.savetxt('Borexino_pseudo_margi_4',bor_data_sc4)

cs_kam_sc4 = data_set(kam_data_sc4,11.83)
cs_bor_sc4 = data_set(bor_data_sc4,11.83)



kam_data_sc3 = get_curves("Kamland","Scalar",3)
bor_data_sc3 = get_curves("Borexino","Scalar",3)

np.savetxt('Kamland_pseudo_margi_3',kam_data_sc3)
np.savetxt('Borexino_pseudo_margi_3',bor_data_sc3)

cs_kam_sc3 = data_set(kam_data_sc3,11.83)
cs_bor_sc3 = data_set(bor_data_sc3,11.83)


#kam_data_sc5 = get_curves("Kamland","Scalar",5)
#bor_data_sc5 = get_curves("Borexino","Scalar",5)

#np.savetxt('Kamland_pseudo_margi_5',kam_data_sc5)
#np.savetxt('Borexino_pseudo_margi_5',bor_data_sc5)

#cs_kam_sc5 = data_set(kam_data_sc5,11.83)
#cs_bor_sc5 = data_set(bor_data_sc5,11.83)

fig,ax = plt.subplots()
plt.text(0.8,0.1,r"3 $\sigma$ C.L.",font="Times New Roman",size=20,transform=ax.transAxes)

ax.set_xlabel(r"$\tau/m$ [s/eV]",font="Times New Roman",size=20)
ax.set_ylabel(r"$\delta$",font="Times New Roman",size=20)
ax.set_title("Scalar; Marginalised",font="Times New Roman",size=20)

ax.set_xlim(1e-6,0.1)
ax.set_xscale('log')
ax.set_ylim(0,1)
ax.yaxis.set_minor_locator(AutoMinorLocator(10))

ax.plot(cs_kam_sc3.X,cs_kam_sc3.Y,color='r',
        label=r'Kamland;$\nu_{3}\rightarrow\bar{\nu}_{1}$',linewidth=2.0)
ax.plot(cs_kam_sc4.X,cs_kam_sc4.Y,'--',color='r',
        label=r'Kamland;$\nu_{2}\rightarrow\bar{\nu}_{3}$',linewidth=2.0)
#ax.plot(cs_kam_sc5.X,cs_kam_sc5.Y,'-.',color='r',
#        label=r'Kamland;$\nu_{1}\rightarrow\bar{\nu}_{3}$',linewidth=2.0)


ax.plot(cs_bor_sc3.X,cs_bor_sc3.Y,color='g',
        label=r'Borexino;$\nu_{3}\rightarrow\bar{\nu}_{1}$',linewidth=2.0)
ax.plot(cs_bor_sc4.X,cs_bor_sc4.Y,'--',color='g',
        label=r'Borexino;$\nu_{2}\rightarrow\bar{\nu}_{3}$',linewidth=2.0)
#ax.plot(cs_bor_sc5.X,cs_bor_sc5.Y,'-.',color='g',
#        label=r'Borexino;$\nu_{1}\rightarrow\bar{\nu}_{3}$',linewidth=2.0)

plt.legend()
#plt.show()
plt.savefig("../plots/compare_scalar_margi_diff_chan2.eps")

