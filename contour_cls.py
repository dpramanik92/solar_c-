import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

import scipy.interpolate as interp

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


Kam_sc_data3 = np.loadtxt("novar_3_Kamland_Scalar.dat",unpack=True)
Kam_sc_data4 = np.loadtxt("novar_4_Kamland_Scalar.dat",unpack=True)
Kam_sc_data5 = np.loadtxt("test7_5_Kamland_Scalar.dat",unpack=True)


Bor_sc_data3 = np.loadtxt("novar_3_Borexino_Scalar.dat",unpack=True)
Bor_sc_data4 = np.loadtxt("novar_4_Borexino_Scalar.dat",unpack=True)
Bor_sc_data5 = np.loadtxt("test7_5_Borexino_Scalar.dat",unpack=True)


Kam_ps_data3 = np.loadtxt("test7_3_Kamland_Pseudo.dat",unpack=True)
Kam_ps_data4 = np.loadtxt("test7_4_Kamland_Pseudo.dat",unpack=True)
Kam_ps_data5 = np.loadtxt("test7_5_Kamland_Pseudo.dat",unpack=True)


Bor_ps_data3 = np.loadtxt("test7_3_Borexino_Pseudo.dat",unpack=True)
Bor_ps_data4 = np.loadtxt("test7_4_Borexino_Pseudo.dat",unpack=True)
Bor_ps_data5 = np.loadtxt("test7_5_Borexino_Pseudo.dat",unpack=True)


cs_kam_sc3 = data_set(Kam_sc_data3,11.83)
cs_kam_sc4 = data_set(Kam_sc_data4,11.83)
#cs_kam_sc5 = data_set(Kam_sc_data5,11.83)

cs_bor_sc3 = data_set(Bor_sc_data3,11.83)
cs_bor_sc4 = data_set(Bor_sc_data4,11.83)
#cs_bor_sc5 = data_set(Bor_sc_data5,11.83)


cs_kam_ps3 = data_set(Kam_ps_data3,11.83)
cs_kam_ps4 = data_set(Kam_ps_data4,11.83)
cs_kam_ps5 = data_set(Kam_ps_data5,11.83)

cs_bor_ps3 = data_set(Bor_ps_data3,11.83)
cs_bor_ps4 = data_set(Bor_ps_data4,11.83)
cs_bor_ps5 = data_set(Bor_ps_data5,11.83)

fig,ax = plt.subplots()
plt.text(0.8,0.1,r"3 $\sigma$ C.L.",font="Times New Roman",size=20,transform=ax.transAxes)

ax.set_xlabel(r"$\tau/m$ [s/eV]",font="Times New Roman",size=20)
ax.set_ylabel(r"$\delta$",font="Times New Roman",size=20)
ax.set_title("Pseudo",font="Times New Roman",size=20)

ax.set_xlim(1e-6,0.1)
ax.set_xscale('log')
ax.set_ylim(0,1)
ax.yaxis.set_minor_locator(AutoMinorLocator(10))

ax.plot(cs_kam_sc3.X,cs_kam_sc3.Y,color='r',
        label=r'Kamland;$\nu_{3}\rightarrow\bar{\nu}_{1}$',linewidth=2.0)
ax.plot(cs_kam_sc4.X,cs_kam_sc4.Y,'--',color='r',
        label=r'Kamland;$\nu_{2}\rightarrow\bar{\nu}_{3}$',linewidth=2.0)
#ax.plot(cs_kam_ps5.X,cs_kam_ps5.Y,'-.',color='r',
#        label=r'Kamland;$\nu_{1}\rightarrow\bar{\nu}_{3}$',linewidth=2.0)


ax.plot(cs_bor_sc3.X,cs_bor_sc3.Y,color='g',
        label=r'Borexino;$\nu_{3}\rightarrow\bar{\nu}_{1}$',linewidth=2.0)
ax.plot(cs_bor_sc4.X,cs_bor_sc4.Y,'--',color='g',
        label=r'Borexino;$\nu_{2}\rightarrow\bar{\nu}_{3}$',linewidth=2.0)
#ax.plot(cs_bor_ps5.X,cs_bor_ps5.Y,'-.',color='g',
#        label=r'Borexino;$\nu_{1}\rightarrow\bar{\nu}_{3}$',linewidth=2.0)

plt.legend()
#plt.show()
plt.savefig("../plots/compare_Scl_diff_chan.eps")
