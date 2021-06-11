import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


x = np.random.normal(0,1,1000)

plt.xlabel(r'$\alpha$',size=20 )
plt.hist(x)
plt.savefig("test.eps")


