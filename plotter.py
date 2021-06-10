import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

import os

def make_e_flux_plot(file_path,file_prefix,which_type,val,clr,style,legend,delt):
    fig,ax = plt.subplots(figsize=(8,6))
    ax.set_xlabel("E (MeV)",size=20)
    ax.set_ylabel(r"$\Phi$",size=20)
    Title = which_type + "; $\delta=$" + str(delt)
    plt.title(Title)
#    _file_pre = file_path + file_prefix

    _plot_pre = file_path + "/"+file_prefix


    flux = np.loadtxt("flux/b8spectrum.txt",unpack=True)
    flux[1] = flux[1];
    ax.plot(flux[0],flux[1],c='k',label='Unoscillated flux')
    for i in range(0,len(val)):

        cmd = "./prob_main "+ str(_plot_pre)+ " " +str(which_type) +" "+str(val[i]) + " "+ str(delt)
        print(cmd)

        os.system(cmd)
        file_name = _plot_pre + "_elec_"+str(which_type)+"_"+str(val[i])+"_"+str(delt)+".dat"
        file_name1 = _plot_pre + "_anti_"+str(which_type)+"_"+str(val[i])+"_"+str(delt)+".dat"
        data =np.loadtxt(file_name,unpack=True)
        ax.plot(data[0],data[1],style[i],c=clr[i],linewidth=1.0,label=legend[i])
        data1 =np.loadtxt(file_name1,unpack=True)
        ax.plot(data1[0],data1[1],style[i],c=clr[i],linewidth=1.0)

    ax.set_ylim(1e-9,0.2)
    ax.set_yscale('log')
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.tick_params(axis='y',which='minor')
    ax.set_xlim(0.1,15.0)
    plt.legend()
    output = "../plots/"+ which_type +"_"+"delta_"+str(delt)+".png"
    plt.savefig(output)






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

    ax.set_ylim(1e-9,1)
    ax.set_yscale('log')
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    plt.tick_params(axis='y',which='minor')
    ax.set_xlim(0.1,15)
    plt.legend()
    plt.show()



def plot_flux():
    fig,ax = plt.subplots()
    ax.set_xlabel("E (MeV)",size=20,font="Times New Roman")
    ax.set_ylabel(r"$\Phi$",size=20,font="Times New Roman")
    flux = np.loadtxt("flux/b8spec-2006.dat",unpack=True)
    flux[1] = 8.5e-4*flux[1]
    data = np.loadtxt("test_elec_Scalar_1e-4_0.95.dat",unpack=True)
    data1 =np.loadtxt("test_anti_Scalar_1e-4_0.95.dat",unpack=True)
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
    os.system("./match")
    data = np.loadtxt("Scalar_05_w.dat",unpack=True)
    plt.plot(data[0],data[1])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.002,1.0)
    plt.ylim(1e-4,3)
    plt.show()

def plot_test():
    os.system("./match")
    data = np.loadtxt("test_flux.dat",unpack=True)
    plt.plot(data[0],data[1],c='r')
    data1 = np.loadtxt("test_anti_Scalar_7.5e-3_0.dat",unpack=True)
    plt.plot(data1[0],data1[1],c='k')
    plt.yscale('log')
    plt.ylim(1e-9,0.2)
    plt.xlim(1.0,15.0)
    plt.show()

def plot_event_kamland():
    fig,ax = plt.subplots(figsize=(8,5))
    data0 = np.loadtxt("../kamland.dat",unpack=True,delimiter=',')
    data = np.loadtxt("exp_data/KamLAND_bkg.dat",unpack=True)
    os.system("./match")
    test = np.loadtxt('test_event_kamland.dat',unpack=True)
    ax.set_xlabel("E(MeV)",size=20,font="Times New Roman")
    ax.set_ylabel("Events/MeV")
    ax.plot(data0[0],data0[1],c='b',label='kamland expected')
    ax.plot(data[0],data[2],'--',c='k',label='Background')
    ax.plot(test[0],test[1],'-.',c='r',label=r'$\tau_{2}/m_{2}=1.2\times10^{-3}$ s/eV')
    ax.set_xlim(7.5,30.0)
    ax.set_ylim(0,12)
    plt.title(r"Scalar, $\delta=0$, SK")
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    plt.legend()
    plt.show()
  #  plt.savefig("../plots/SK_IV_match.png")




def plot_event_borexino():
    fig,ax = plt.subplots(figsize=(8,5))
    data0 = np.loadtxt("../kamland.dat",unpack=True,delimiter=',')
    data = np.loadtxt("exp_data/KamLAND_bkg.dat",unpack=True)
    os.system("./match")
    test = np.loadtxt('test_event_kamland.dat',unpack=True)
    ax.set_xlabel("E(MeV)",size=20,font="Times New Roman")
    ax.set_ylabel("Events/MeV")
    ax.plot(data0[0],data0[1],c='b',label='kamland expected')
    ax.plot(data[0],data[2],'--',c='k',label='Background')
    ax.plot(test[0],test[1],'-.',c='r',label=r'$\tau_{2}/m_{2}=1.2\times10^{-3}$ s/eV')
    ax.set_xlim(7.5,30.0)
    ax.set_ylim(0,12)
    plt.title(r"Scalar, $\delta=0$, SK")
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    plt.legend()
    plt.show()


def plot_event_SK():
    fig,ax = plt.subplots(figsize=(8,5))
    data0 = np.loadtxt("../sk-fit.dat",unpack=True,delimiter=',')
    data = np.loadtxt("exp_data/SuperK-IV_bkg.dat",unpack=True)
    os.system("./match")
    test = np.loadtxt('test_event_SK.dat',unpack=True)
    ax.set_xlabel("E(MeV)",size=20,font="Times New Roman")
    ax.set_ylabel("Events/MeV")
    ax.step(data0[0],data0[1],c='b',label='Super-K expected',where='mid')
    ax.step(data[0],data[2],'--',c='k',label='Background',where='mid')
    ax.step(test[0],test[1],'-.',c='r',label=r'$\tau_{2}/m_{2}=1.2\times10^{-3}$ s/eV',where='mid')
    ax.set_xlim(7.5,30.0)
    ax.set_ylim(0,50)
    plt.title(r"Scalar, $\delta=0$, SK")
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    plt.legend()
    plt.show()


def plot_test_flux():
    fig,ax = plt.subplots(figsize=(8,6))
    data0 = np.loadtxt("test_anti_Scalar_7.5e-3_0.dat",unpack=True)
    data1 = np.loadtxt("test_anti_Scalar_6.8e-3_0.3.dat",unpack=True)
    data2 = np.loadtxt("test_anti_Scalar_4.1e-3_0.7.dat",unpack=True)
    data3 = np.loadtxt("test_anti_Scalar_1.2e-3_0.9.dat",unpack=True)
    data4 = np.loadtxt("test_anti_Scalar_2.3e-4_0.96.dat",unpack=True)
    data5 = np.loadtxt("test_anti_Scalar_1e-8_0.99.dat",unpack=True)

    test1 = np.loadtxt("test_flux.dat",unpack=True)


    ax.set_xlabel("E[MeV]",size=20,font="Times New Roman")
    ax.set_ylabel("Normalized Flux",size=20,font="Times New Roman")

    ax.plot(data0[0],data0[1],c='r',label="7.5e-3;0")
    ax.plot(data1[0],data1[1],c='g',label="6.8e-3;0.3")
    ax.plot(data2[0],data2[1],c='b',label="4.1e-3;0.7")
    ax.plot(data3[0],data3[1],c='c',label="1.2e-3;0.9")
    ax.plot(data4[0],data4[1],c='m',label="2.3e-4;0.96",linewidth=2.0)
    ax.plot(data5[0],data5[1],c='y',label="1e-8;0.99",linewidth=2.0)
    ax.plot(test1[0],test1[1],'--',c='k',label="conversion=5.3e-5",linewidth=2.0)

    ax.set_xlim(1.0,15)
    ax.set_ylim(1e-9,0.2)
    ax.set_yscale('log')


    plt.legend()
    plt.savefig("../plots/test_flux.png")




def plot_chi2():
    os.system("./chi2_main")
    fig,ax = plt.subplots(figsize=(8,6))
    data0 = np.loadtxt("test_chi2.dat",unpack=True)

    chi2_min = np.min(data0[1])

    ax.plot(data0[0],data0[1]-chi2_min,c='r')
    ax.set_xlabel(r'$\tau_{2}/m_{2}$ [eV$^{2}$]')
    ax.set_ylabel(r'$\Delta\chi^{2}$')

    ax.set_xscale('log');
    ax.set_ylim(0,50);

    plt.show()

'''
val = np.array([1e-4,1e-3,1e-2])
clr = ['r','g','b']
style = ["--",":","-."]
legend = [r"$\tau_{2}/m_{2}=10^{-4} s/eV$",r"$\tau_{2}/m_{2}=10^{-3} s/eV$",r"$\tau_{2}/m_{2}=10^{-2} s/eV$"]

#os.system("./event_main")
os.system("make")

#plot_flux()
delt = [0,0.5,0.96,0.99,0.999]
which_type = ["Scalar","Pseudo","Mixed"]

for d in delt:
    for wh in which_type:
        make_e_flux_plot("Probability_data","test_new",wh,val,clr,style,legend,d)
'''
#os.system("make")
#plot_event_kamland()
#plot_test()
#plot_test_flux()
#plot_event_SK()
#plot_chi2()

#data = np.loadtxt("test.dat",unpack=True)
#data1 = np.loadtxt("test1.dat",unpack=True)
#data2 = np.loadtxt("test2.dat",unpack=True)

#plt.plot(data[0],data[1],c='r')
#plt.plot(data1[0],data1[1],c='g')
#plt.plot(data2[0],data2[1],c='k')

'''
import scipy.interpolate

N = 1000

x,y,z = np.loadtxt('test_SK-IV_Scalar.dat',unpack=True)

x1,y1,z1 = np.loadtxt('test_SK-IV_Pseudo.dat',unpack=True)

chi2_min = np.min(z)
z = z-chi2_min



chi2_min1 = np.min(z1)
z1 = z1-chi2_min1





print(chi2_min)
print(chi2_min1)

xi = np.linspace(x.min(),x.max(),N)
yi = np.linspace(y.min(),y.max(),N)
zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

xi1 = np.linspace(x1.min(),x1.max(),N)
yi1 = np.linspace(y1.min(),y1.max(),N)
zi1 = scipy.interpolate.griddata((x1, y1), z1, (xi1[None,:], yi1[:,None]), method='cubic')

fig = plt.figure()
plt.xlabel(r'$\delta$',font="Times New Roman",size=18)
plt.ylabel(r'$\tau_{2}/m_{2}$ s/eV',font="Times New Roman",size=18)
plt.yscale('log')
plt.xlim(0.0,1.0)
plt.ylim(1e-5,1.0)
plt.contour(xi,10**yi,zi,levels=[5.99],colors='r',clabel='Scalar')
plt.contour(xi1,10**yi1,zi1,levels=[5.99],colors='g',clabel='Pseudo scalar')
#plt.contour(xi,yi,zi)
plt.legend()
#plt.show()
plt.savefig('../test1.png')
'''

data1 = np.loadtxt('corr_SK-IV_Scalar.dat',unpack=True)
data2 = np.loadtxt('uncorr_SK-IV_Scalar.dat',unpack=True)
data3 = np.loadtxt('nosys_SK-IV_Scalar.dat',unpack=True)

chi2_min = np.min(data1[1])
chi2_min1 = np.min(data2[1])
chi2_min2 = np.min(data3[1])

data1[1] = data1[1]- chi2_min
data2[1] = data2[1] - chi2_min1
data3[1] = data3[1] - chi2_min2


plt.xlim(1e-6,0.1)
plt.xscale('log')
plt.xlabel(r"$\tau_{2}/m_{2}$[$eV^{2}}$]",size=20,font="Times New Roman")
plt.ylabel(r"$\Delta \chi^{2}$",size=20,font="Times New Roman")
plt.ylim(0,10)
plt.title(r'Scalar; $\delta=0.9$',size=20,font="Times New Roman")
plt.plot(data1[0],data1[1],c='r',label='Correlated systematics',linewidth=2.0)
plt.plot(data2[0],data2[1],c='k',label='Uncorrelated systematics',linewidth=2.0)
plt.plot(data3[0],data3[1],c='g',label='No systematics',linewidth=2.0)
plt.legend()
plt.savefig("../plots/sys_sk_iv.png")
#plt.show()
