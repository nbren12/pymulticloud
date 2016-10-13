#!/usr/bin/env python
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.insert(0, "../../")

from python.read import read_data
from python.cmt import calc_du

data = read_data("./data")


# u plot
fig, axs= plt.subplots(1,3)
for i in range(3):
    axs[i].pcolormesh(data['u'][:,i,:], cmap='bwr')
    axs[i].set_title("CMT={0}".format(i))
    plt.savefig('u.png')

plt.figure()
umean = data['u'][:,:,2:-2].mean(axis=0)
plt.plot(umean.T)
plt.legend([0,1,2])
plt.savefig("umean.png")

plt.figure()
dul, duhi = calc_du(umean)
plt.plot(dul, label='dulow')
plt.plot(duhi, label='duhi')
plt.legend()
plt.savefig("du.png")


plt.figure()
plt.plot(np.sqrt(data['u'][:,:,2:-2]**2).mean(axis=0).T)
plt.legend([0,1,2])
plt.savefig("urms.png")


plt.figure()
plt.pcolormesh(data['lmd'])
plt.colorbar()
plt.savefig("lmd.png")
