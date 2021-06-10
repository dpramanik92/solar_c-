import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

'''
data = np.loadtxt("Prob_nu3nu1_1e-4.dat",unpack=True)
data1 = np.loadtxt("Prob_nu3nu1_1e-6.dat",unpack=True)


fig,ax = plt.subplots()
ax.set_xlabel("E [MeV]",font="Times New Roman",size=20)
ax.set_ylabel(r"$\Phi$",font="Times New Roman",size=20)
ax.set_title(r"$\nu_{3}\rightarrow\bar{\nu}_{1}}$;$\delta$=0.9;Scalar",font="Times New Roman",size=20)
ax.set_yscale('log')
ax.set_ylim(1e-9,0.2)
ax.set_xlim(0.1,15.0)
ax.xaxis.set_minor_locator(AutoMinorLocator(10))
ax.plot(data[0],data[1],label=r"$\tau_{3}/m_{3}=1\times10^{-4} [s/eV]$",c='r',linewidth=2.0)
ax.plot(data1[0],data1[1],label=r"$\tau_{3}/m_{3}=1\times10^{-6} [s/eV]$",c='b',linewidth=2.0)
plt.legend()
plt.savefig("../plots/test_nu3nu1_prob.eps")
#plt.show()
'''

'''
data = np.loadtxt("nu3nu1_Kamland_Scalar.dat",unpack=True)
data1 = np.loadtxt("nu3nu1_Borexino_Scalar.dat",unpack=True)
data2 = np.loadtxt("nu3nu1_SK-IV_Scalar.dat",unpack=True)

chi2_min = np.min(data[1])
chi2_min1 = np.min(data1[1])
chi2_min2 = np.min(data2[1])

data[1] = data[1] - chi2_min
data1[1] = data1[1] - chi2_min1
data2[1] = data2[1] - chi2_min2

fig,ax = plt.subplots()
ax.set_xlabel(r"$\tau_{3}/m_{3}$ [s/eV]",font="Times New Roman",size=20)
ax.set_ylabel(r"$\Delta \chi^{2}$",font="Times New Roman",size=20)
ax.set_title(r"$\nu_{3}\rightarrow\bar{\nu}_{1}}$;$\delta$=0.1;Scalar",font="Times New Roman",size=20)
ax.set_xscale('log')
ax.set_xlim(1e-6,0.1)
ax.set_ylim(0.0,20.0)
ax.yaxis.set_minor_locator(AutoMinorLocator(10))
ax.plot(data[0],data[1],label="Kamland",c='r',linewidth=2.0)
ax.plot(data1[0],data1[1],label="Borexino",c='b',linewidth=2.0)
ax.plot(data2[0],data2[1],label="SK-IV",c='g',linewidth=2.0)
plt.legend()
plt.savefig("../plots/test_nu3nu1_chi2nosys.eps")
#plt.show()

'''
'''
fig,ax = plt.subplots()
x,y,z = np.loadtxt("test4_Borexino_Scalar.dat",unpack=True)

import scipy.interpolate

N = 1000

x = np.log10(x)
chi2_min = np.min(z)
z = z - chi2_min


xi = np.linspace(x.min(),x.max(),N)
yi = np.linspace(y.min(),y.max(),N)
zi = scipy.interpolate.griddata((x,y),z,(xi[None,:],yi[:,None]))

plt.text(0.6,0.8,"Borexino",font="Times New Roman",size=20,transform=ax.transAxes)
ax.set_xlabel(r"$\tau_{3}/m_{3}$ [s/ev]",font="times new roman",size=20)
ax.set_ylabel(r"$\delta$",font="times new roman",size=20)
ax.set_title(r"$\nu_{3}\rightarrow\bar{\nu}_{1}$; Scalar",font="Times New Roman",size=20)
ax.set_xlim(1e-6,0.1)
ax.set_xscale('log')
ax.set_ylim(0,1)
ax.yaxis.set_minor_locator(AutoMinorLocator(10))
levels=[2.3,6.18,11.83]
cs=ax.contour(10**xi,yi,zi,levels,colors=['r','g','b'])
cs.levels=[r"1$\sigma$",r"2$\sigma$",r"3$\sigma$"]
ax.clabel(cs,cs.levels)
plt.savefig("../plots/delt_tau3_Borexino_Scalar.eps")
'''

'''
data = np.loadtxt("marg_nu3nu1_Kamland_Pseudo.dat",unpack=True)
data1 = np.loadtxt("marg_pr_nu3nu1_Kamland_Pseudo.dat",unpack=True)
data2 = np.loadtxt("marg_nu3nu1_Borexino_Pseudo.dat",unpack=True)
data3 = np.loadtxt("marg_pr_nu3nu1_Borexino_Pseudo.dat",unpack=True)

chi2_min = np.min(data[1])
chi2_min1 = np.min(data1[1])
chi2_min2 = np.min(data2[1])
chi2_min3 = np.min(data3[1])

data[1] = data[1] - chi2_min
data1[1] = data1[1] - chi2_min1
data2[1] = data2[1] - chi2_min2
data3[1] = data3[1] - chi2_min3

fig,ax = plt.subplots()
ax.set_xlabel(r"$\tau_{3}/m_{3}$ [s/eV]",font="Times New Roman",size=20)
ax.set_ylabel(r"$\Delta \chi^{2}$",font="Times New Roman",size=20)
ax.set_title(r"$\nu_{3}\rightarrow\bar{\nu}_{1}}$;$\delta$=0.1;Pseudo",font="Times New Roman",size=20)
ax.set_xscale('log')
ax.set_xlim(1e-6,0.1)
ax.set_ylim(0.0,20.0)
ax.yaxis.set_minor_locator(AutoMinorLocator(10))
ax.plot(data[0],data[1],label="Kamland wo prior",c='r',linestyle='-',linewidth=2.0)
ax.plot(data1[0],data1[1],label="Kamland w prior",c='r',linestyle='--',linewidth=2.0)
ax.plot(data2[0],data2[1],label="Borexino wo prior",c='b',linestyle='-',linewidth=2.0)
ax.plot(data3[0],data3[1],label="Borexino w prior",c='b',linestyle='--',linewidth=2.0)
plt.legend()
plt.savefig("../plots/marg_nu3nu1_chi2nosys_pseudo.eps")
#plt.show()
'''
'''
plt.yscale('log')
data = np.loadtxt("Prob_nu3nu1_1e-4.dat",unpack=True)
flux = np.loadtxt("flux/b8spectrum.txt",unpack=True)

plt.ylim(1e-9,1)
plt.plot(data[0],data[1])
plt.plot(data[0],data[2])
plt.plot(flux[0],flux[1])
plt.show()
'''

data = np.loadtxt("test_fl.dat",unpack=True)
plt.xlim(1e-6,1e-1)
plt.xscale('log')
plt.plot(data[0],data[1])
plt.plot(data[0],data[2])
plt.show()




