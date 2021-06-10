import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)



x_kam,y_kam,z_kam = np.loadtxt("test4_Kamland_Scalar.dat",unpack=True)
x_bor,y_bor,z_bor = np.loadtxt("test4_Borexino_Scalar.dat",unpack=True)


x_kam1,y_kam1,z_kam1 = np.loadtxt("var3_pr_Kamland_Scalar.dat",unpack=True)
x_bor1,y_bor1,z_bor1 = np.loadtxt("var3_pr_Borexino_Scalar.dat",unpack=True)

z_comb = z_kam + z_bor
z_comb1 = z_kam1 + z_bor1

import scipy.interpolate

N = 1000

x = np.log10(x_kam)

chi2_min_kam = np.min(z_kam)
chi2_min_bor = np.min(z_bor)
chi2_min_comb = np.min(z_comb)

chi2_min_kam1 = np.min(z_kam1)
chi2_min_bor1 = np.min(z_bor1)
chi2_min_comb1 = np.min(z_comb1)

z_kam = z_kam - chi2_min_kam
z_bor = z_bor - chi2_min_bor
z_comb = z_comb - chi2_min_comb

z_kam1 = z_kam1 - chi2_min_kam1
z_bor1 = z_bor1 - chi2_min_bor1
z_comb1 = z_comb1 - chi2_min_comb1

xi = np.linspace(x.min(),x.max(),N)
yi = np.linspace(y_kam.min(),y_kam.max(),N)

zi_kam = scipy.interpolate.griddata((x,y_kam),z_kam,(xi[None,:],yi[:,None]))
zi_bor = scipy.interpolate.griddata((x,y_kam),z_bor,(xi[None,:],yi[:,None]))
zi_comb = scipy.interpolate.griddata((x,y_kam),z_comb,(xi[None,:],yi[:,None]))

zi_kam1 = scipy.interpolate.griddata((x,y_kam),z_kam1,(xi[None,:],yi[:,None]))
zi_bor1 = scipy.interpolate.griddata((x,y_kam),z_bor1,(xi[None,:],yi[:,None]))
zi_comb1 = scipy.interpolate.griddata((x,y_kam),z_comb1,(xi[None,:],yi[:,None]))

cs_kam = plt.contour(10**xi,yi,zi_kam,[11.83])
p_kam = cs_kam.collections[0].get_paths()[0]

v_kam = p_kam.vertices
X_kam = v_kam[:,0]
Y_kam = v_kam[:,1]

cs_bor = plt.contour(10**xi,yi,zi_bor,[11.83])
p_bor = cs_bor.collections[0].get_paths()[0]

v_bor = p_bor.vertices
X_bor = v_bor[:,0]
Y_bor = v_bor[:,1]

cs_com = plt.contour(10**xi,yi,zi_comb,[11.83])
p_com = cs_com.collections[0].get_paths()[0]

v_com = p_com.vertices
X_com = v_com[:,0]
Y_com = v_com[:,1]

cs_kam1 = plt.contour(10**xi,yi,zi_kam1,[11.83])
p_kam1 = cs_kam1.collections[0].get_paths()[0]

v_kam1 = p_kam1.vertices
X_kam1 = v_kam1[:,0]
Y_kam1 = v_kam1[:,1]

cs_bor1 = plt.contour(10**xi,yi,zi_bor1,[11.83])
p_bor1 = cs_bor1.collections[0].get_paths()[0]

v_bor1 = p_bor1.vertices
X_bor1 = v_bor1[:,0]
Y_bor1 = v_bor1[:,1]

cs_com1 = plt.contour(10**xi,yi,zi_comb1,[11.83])
p_com1 = cs_com1.collections[0].get_paths()[0]

v_com1 = p_com1.vertices
X_com1 = v_com1[:,0]
Y_com1 = v_com1[:,1]
fig,ax = plt.subplots()

plt.text(0.8,0.1,r"3 $\sigma$ C.L.",font="Times New Roman",size=20,transform=ax.transAxes)

ax.set_xlabel(r"$\tau_{3}/m_{3}$ [s/eV]",font="Times New Roman",size=20)
ax.set_ylabel(r"$\delta$",font="Times New Roman",size=20)
ax.set_title(r"$\nu_{3}\rightarrow\bar{\nu}_{1}$; Scalar",font="Times New Roman",size=20)

ax.set_xlim(1e-6,0.1)
ax.set_xscale('log')
ax.set_ylim(0,1)
ax.yaxis.set_minor_locator(AutoMinorLocator(10))

plt.xscale('log')
plt.plot(X_kam,Y_kam,color='r',label='Kamland',linewidth=2.0)
plt.plot(X_bor,Y_bor,color='g',label='Borexino',linewidth=2.0)
plt.plot(X_com,Y_com,color='k',label='Combined',linewidth=2.0)
plt.plot(X_kam1,Y_kam1,'--',color='r',label='Kamland w prior',linewidth=2.0)
plt.plot(X_bor1,Y_bor1,'--',color='g',label='Borexino w prior',linewidth=2.0)
plt.plot(X_com1,Y_com1,'--',color='k',label='Combined w prior',linewidth=2.0)
#plt.fill_between(X_kam,Y_kam,0,facecolor="orange",color="orange",alpha=0.1)
#plt.fill_between(X_bor,Y_bor,0,facecolor="red",color="red",alpha=0.1)
#plt.fill_between(X_com,Y_com,0,facecolor="green",color="green",alpha=0.1)
#plt.plot(X_bor,Y_bor)
#plt.plot(X_com,Y_com)lt.plot(X_com,Y_com)
plt.legend()
plt.savefig("../plots/delt_tau3_all_Scalar_pr.eps")
#plt.show()
