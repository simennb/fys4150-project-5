import numpy as np
import matplotlib.pyplot as plt
from plot_cluster import read_file

N = 300
dt = 0.01
tcoll = 1.0
eps = 0.2
eqtime = 75


filename = "../benchmarks/pos_N%d_dt%1.6f_tcoll_%1.6f_eps%1.6f.xyz" %(N, dt, tcoll, eps)



pos, mass = read_file(filename)

r = np.zeros(N)

counter = 0
for i in pos:
    r[counter] = np.sqrt(i[eqtime, 0]**2 + i[eqtime, 1]**2 + i[eqtime, 2]**2)
    counter += 1


#print r
#print np.mean(r)


dr = 70
radii = np.linspace(0, 20, dr)
dV = np.zeros(len(radii)-1)
number_density = np.zeros(len(dV))
mass_density = np.zeros(len(dV))

for i in range(0, len(radii)-1):
    dV[i] = 4./3*np.pi*radii[i+1]**3 - 4./3*np.pi*radii[i]**3
    for j in range(0, len(r)):
        if r[j] < radii[i+1] and r[j] >= radii[i]:
            mass_density[i] += mass[j]
            number_density[i] += 1

#print mass_density
print number_density

print mass_density

n, bins, patches = plt.hist(radii[:-1], bins = radii, weights = number_density)

plt.show()
