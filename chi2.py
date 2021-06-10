import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

import os

def calculate_lifetime_limit(file_path,file_prefix,which_type,delta):
    _prefix = file_path +"/"+file_prefix
    os.system("make")
    for Ex in ["SK-IV","Borexino","Kamland"]:
        cmd = "./chi2_main " + _prefix + " " + Ex + " " + which_type + " " + str(delta)
        os.system(cmd)


    fig,ax = plt.subplots(figsize=(8,6))
    ax.set_xlabel(r"$\tau_{2}/m_{2}$ [s/eV]",size=20,font="Times New Roman")
    ax.set_ylabel(r"$\Delta \chi^{2}$",size=20,font="Times New Roman")

    ax.set_xscale('log')
    ax.set_xlim(1e-5,1e-1)
    ax.set_ylim(0,20)
    Title = which_type + r";$\delta = $" + str(delta)
    ax.set_title(Title)
#    ax.axhline(9,c='k',linewidth=0.5)

    clr = ['r','b','g']

    Exp = ["Kamland","Borexino","SK-IV"]
    combo = []
    X = []
    for i in range(0,len(Exp)):
        _exp = Exp[i]
        input_file = _prefix + "_" + _exp + "_" + which_type + "_" + str(delta) +".dat"
        print(input_file)
        data = np.loadtxt(input_file,unpack=True)
        combo.append(data)
        X = data[0]
        chi2_min = np.min(data[1])
        data[1] = data[1] - chi2_min

        ax.plot(data[0],data[1],c=clr[i],label=Exp[i],linewidth=1.0)

    combine = combo[0][1] + combo[1][1] + combo[2][1]

    comb_min = np.min(combine)

    combine = combine - comb_min
 #   ax.plot(combo[0][0],combine,'-.',c='k',label="Combined",linewidth=2.0)

    ax.yaxis.set_minor_locator(AutoMinorLocator(10))

    plt.legend()
    plt.show()




calculate_lifetime_limit("chi2_data","test","Scalar",0.9)
