import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys


G = sys.argv[1]

df= pd.read_table(G, index_col = 0, sep = ' *',
        names=['u1', 'u2', 'theta1', 'theta2',
            'theta_eb', 'q', 'hs', 'hc', 'hd',
            'fcls', 'fdls','fsls']
        )

ax = plt.subplot(311)
df[['fcls', 'fdls', 'fsls' ]].plot(ax = ax,  style='-', color=['b', 'k', 'g'])

ax = plt.subplot(312, sharex=ax)
df[['u1', 'u2']].plot(ax = ax, style='-', color=['b', 'k', 'g'])

ax = plt.subplot(313, sharex=ax)
df[['theta1', 'theta2', 'theta_eb']].plot(ax = ax, style='-', color=['b', 'k', 'g'])

plt.xlabel('Days')

# these two numbers should match if nstoch = 110
print 1/df.fcls[df.fcls>0].min()
print 110**2

plt.show()
