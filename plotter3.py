import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                 AutoMinorLocator)

class plot_2d:
    def __init__(self,datafile):
        self.data = np.loadtxt(datafile)
        print("[Messeage]: data loaded!!")

    def set_canvas(self):
        self.fig,self.ax = plt.subplots()



