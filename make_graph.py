import numpy
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.use('Agg')
import pylab as pl
import matplotlib.pyplot as plt
import glob

pl.figure(figsize=(7, 7))

print_every = int(raw_input())

all_lists = []
x = []

for file_name in glob.glob('values_of_s*.out'):
    with open(file_name) as f:
        l = []
        append_to_x = (len(x) == 0)
        for i, s in enumerate(f):
            if append_to_x is True:
                x.append(print_every * i + 1)
            l.append(int(s))
        all_lists.append(l)

fig, ax = plt.subplots()
for l in all_lists:
    plt.plot(x, l, linewidth=0.1)
pl.xlim([0,max(x) + 1])
ax.grid(True, which='both')
plt.savefig("graph_of_l.pdf")
