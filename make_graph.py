import numpy
import pylab as pl
import matplotlib.pyplot as plt

pl.figure(figsize=(7, 7))  # Don't create a humongous figure

x = []
k = []
l = []

with open("liczby.out") as f:
    for i, s in enumerate(f):
        a = [int(z) for z in s.strip('\n').split(" ")]
        x.append(500 * i + 1)
        k.append(a[0])
        l.append(a[1])


plt.plot(x, l, linewidth=0.1)
plt.plot(x, k, linewidth=0.1)
plt.savefig("graph.pdf")
