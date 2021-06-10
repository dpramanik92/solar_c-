import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("chi2_data/marg_test_Kamland_Scalar_0.9_col_0.dat",unpack=True)

#plt.xscale('log')
plt.plot(data[0],data[1],c='r')
plt.show()
