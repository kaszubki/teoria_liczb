import numpy
import pylab as pl
import matplotlib.pyplot as plt
import seaborn
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

pl.figure(figsize=(7, 7))  # Don't create a humongous figure
seaborn.set(style='ticks')

x = []
f = []

with open("values_of_f_function.out") as plik:
    for i, s in enumerate(plik):
        a = [int(z) for z in s.strip('\n').split(" ")]
        x.append(i + 1)
        f.append(a[0])


fig, ax = plt.subplots()
ax.set_ylabel('F(n)')
ax.set_xlabel('n')
ax.plot(x, f, linewidth=1)
pl.xlim([0,max(x) + 1])
plt.xticks(numpy.arange(0, max(x)+1, 10))
plt.yticks(numpy.arange(1, min(f) - 1, -5))
ax.set_aspect('equal')
ax.grid(True, which='both')
seaborn.despine(ax=ax, offset=0) # the important part here
plt.savefig("graph_of_f.pdf")
print f
