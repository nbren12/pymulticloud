import sys

import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.insert(0, "../../")

from python.read import read_data

data = read_data("./data")

from IPython import embed; embed()

fig, axs= plt.subplots(1,3)

for i in range(3):
    axs[i].pcolormesh(data['u'][:,i,:], cmap='Greys')
    axs[i].set_title("CMT={0}".format(i))
    plt.savefig('scmt.png')
