import numpy as np
import matplotlib.pyplot as plt
from plot_cluster import read_file
import sys

N = 100
dt = 0.001
tcoll = 3.0
eps = 0.2
eqtime = 75


filename = "../benchmarks/pos_N%d_dt%1.6f_tcoll_%1.6f_eps%1.6f.xyz" %(N, dt, tcoll, eps)

pos, mass, potential, kinetic = read_file(filename)

total = np.array(potential) + np.array(kinetic)

tot_potEnergyPerPart = []
tot_kinEnergyPerPart = []
tot_EnergyPerPart = []

for i in range(0, len(potential)):
    tot_potEnergyPerPart.append(sum(potential[i]))
    tot_kinEnergyPerPart.append(sum(kinetic[i]))
    tot_EnergyPerPart.append(sum(total[i]))

tot_potEnergy = sum(tot_potEnergyPerPart)
tot_kinEnergy = sum(tot_kinEnergyPerPart)
tot_Energy = sum(tot_EnergyPerPart)

x = np.linspace(0, tcoll, 100)

print len(x)
print len(total)
print len(total[i])



counterarray = np.zeros(len(x))
for i in range(0, len(x)):
    for j in range(0, len(total)):
        if total[j][i] > 0:
            counterarray[i] += 1

plt.plot(x, counterarray/N)
plt.show()
