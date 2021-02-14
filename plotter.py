import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

def plot_probability():
    fg,ax = plt.subplots()
    ax.set_xlabel("E (MeV)",size=20)
    ax.set_ylabel(r"$P_{ee}$",size=20)
    data3 = np.loadtxt("prob_test_elec_Scalar_1e-5_0.96.dat",unpack=True)
    data= np.loadtxt("prob_test_elec_Scalar_0_0.96.dat",unpack=True)
    data2 = np.loadtxt("prob_test_elec_Scalar_1e-4_0.96.dat",unpack=True)
    data1 = np.loadtxt("prob_test_elec_Scalar_1e-3_0.96.dat",unpack=True)
    ax.plot(data[0],data[1],c='k',linewidth=1.0,
            label=r'$\tau_{2}/m_{2} \rightarrow \infty$')
    ax.plot(data1[0],data1[1],'--',c='r',linewidth=1.0,
            label=r'$\tau_{2}/m_{2}=10^{-3} s.eV^{-1}$')
    ax.plot(data2[0],data2[1],'-.',c='b',linewidth=1.0,
            label=r'$\tau_{2}/m_{2} = 10^{-4} s.eV^{-1}$')
    ax.plot(data3[0],data3[1],':',c='k',linewidth=1.0,
            label=r'$\tau_{2}/m_{2} = 10^{-5} s.eV^{-1}$')
    ax.set_xscale('log')
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax.set_ylim(0.0,1.0)
    ax.set_xlim(0.1,20)
    plt.legend()
    plt.show()



def plot_flux():
    fig,ax = plt.subplots()
    ax.set_xlabel("E (MeV)",size=20,font="Times New Roman")
    ax.set_ylabel(r"$\Phi$",size=20,font="Times New Roman")
    flux = np.loadtxt("flux/b8spec-2006.dat",unpack=True)
    flux[1] = 8.5e-4*flux[1]
    data = np.loadtxt("flux_test_elec_Scalar_1e-2_0.3.dat",unpack=True)
    data1 =np.loadtxt("flux_test_anti_Scalar_1e-2_0.3.dat",unpack=True)
    ax.plot(flux[0],flux[1],':',c='k',linewidth=1.0)
    ax.plot(data[0],data[1],'-',c='b',linewidth=1.0)
    ax.plot(data1[0],data1[1],'--',c='b',linewidth=1.0)
    ax.set_ylim(1e-9,0.2)
    ax.set_yscale('log')
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    plt.tick_params(axis='y',which='minor')
    ax.set_xlim(1.0,17.0)
#    ax.set_ylim(1e-9,1)
    plt.show()

def plot_w():
    data = np.loadtxt("Scalar_05_w.dat",unpack=True)
    plt.plot(data[0],data[1])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.002,1.0)
    plt.ylim(1e-4,3)
    plt.show()

plot_flux()
