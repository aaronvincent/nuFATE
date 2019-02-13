""" Example of showing the maximum difference betwen nuFATE and ODE solver at different zenith
"""

import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(6,5))
for i in range(1,4):
    diff = np.loadtxt('../../resources/g'+str(i)+'.txt')
    plt.plot(diff[:,0],diff[:,1])
plt.xlabel(r"Zenith Angle")
plt.ylabel(r"Max Diffrence between nuFATE and ODE")
plt.legend(["g=1 NuTauBar", "g=2 NuTauBar", "g=3 NuTauBar"])
plt.grid()
plt.show()
