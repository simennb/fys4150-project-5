import numpy as np
import matplotlib.pyplot as plt
from plot_cluster import read_file
import sys

def radDist(N, dt, tcoll, eps):

    filename = "../benchmarks/pos_N%d_dt%1.6f_tcoll_%1.6f_eps%1.6f.xyz" %(N, dt, tcoll, eps)

    pos, mass, potential, kinetic = read_file(filename)

    r = np.zeros(N)

    counter = 0
    for i in pos:
        r[counter] = np.sqrt(i[eqtime, 0]**2 + i[eqtime, 1]**2 + i[eqtime, 2]**2)
        counter += 1


    #print np.std(r)
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

    return radii, number_density/dV, mass_density/dV, np.mean(r), np.std(r), mass

def nnn(r, r0, n0):
    return n0/(1 + (r/float(r0))**4)

def ppp(r, r0, p0):
    return p0/(r/float(r0)*(1 + (r/r0))**2)

N = [100, 200, 300, 400, 500, 600, 1000]
dt = 0.001
tcoll = 5.0
eps = 0.1
eqtime = -1

r0 = np.linspace(0.0001, 10, 1000)
n_sums = np.zeros(len(r0))
#p_sums = np.zeros((len(r0),3))

meanr_array = np.zeros(len(N))
stdr_array = np.zeros(len(N))

for i in range(0, len(N)):
    radii, number_density, mass_density, meanr, stdr, mass = radDist(N[i], dt, tcoll, eps)
    meanr_array[i] = meanr
    stdr_array[i] = stdr
    
    plt.figure(figsize = (9,7))
    n, bins, patches = plt.hist(radii[:-1], bins = radii, weights = number_density, normed = True)

    # Least squares method searching for r0
    n0 = max(n)
    n0_index = np.argmax(n)
    r0real = N[i]**(1./3)
    alpha = np.linspace(0.01,2,1000)
    for alpha_ in range(len(alpha)):
        n_array = nnn(bins[:-1]-bins[n0_index], alpha[alpha_]*r0real, n0)
        n_sums[alpha_] = sum((n_array - n)**2)

    min_index = np.argmin(n_sums)
    r0real *= alpha[min_index]

    plt.plot(bins+bins[np.argmax(n)], nnn(bins, r0real, n0),'k',lw=3)
    plt.title('Number density per radii for %d particles' %N[i], fontsize = 18)
    plt.xlabel('Radius [light years]', fontsize = 18)
    plt.ylabel('Number density [1/ly^3]', fontsize = 18)
    plt.xlim([0,10])
    plt.savefig("../figures/taskf/radialDens_N%d.png" %N[i])

    ############################### Mass density
    plt.figure(figsize = (9,7))
    n, bins, patches = plt.hist(radii[:-1], bins = radii, weights = mass_density, normed = True)

    # Least squares method searching for r0
    p0 = sum(mass)/(4./3*np.pi*radii[-1]**3)
    print p0
    p0_index = np.argmax(n)
    r0real = N[i]**(1./3)
    rs = np.linspace(1,20,1000)
    for rs_index in range(len(rs)):
        n_array = ppp(bins[:-1], rs[rs_index], p0)
        n_sums[rs_index] = sum((n_array - n)**2)

    min_index = np.argmin(n_sums)
    r0real = rs[min_index]

#    plt.plot(bins+bins[np.argmax(n)], ppp(bins, r0real, p0),'k',lw=3)
    plt.title('Mass density per radii for %d particles' %N[i], fontsize = 18)
    plt.xlabel('Radius [light years]', fontsize = 18)
    plt.ylabel('Mass density [M_sun/ly^3]', fontsize = 18)
    plt.xlim([0,10])
    plt.savefig("../figures/taskf/massDens_N%d.png" %N[i])

    if len(sys.argv) > 1 and sys.argv[1] == "show":
        plt.show()

'''
ylog = np.log(meanr_array)

first = np.polyfit(N, meanr_array, 1)
second = np.polyfit(N, 1./meanr_array, 2)
exponential = np.polyfit(N, ylog, 1)

x = np.linspace(1, N[-1] + 100, 1000)


yfirst = np.polyval(first, x)
ysecond = np.polyval(second, x)
yexponential = np.polyval(exponential, x)

ylog2 = np.log(stdr_array)

first2 = np.polyfit(N, stdr_array, 1)
second2 = np.polyfit(N, 1./stdr_array, 2)
exponential2 = np.polyfit(N, ylog, 1)

yfirst2 = np.polyval(first2, x)
ysecond2 = np.polyval(second2, x)
yexponential2 = np.polyval(exponential2, x)


plt.figure(3, figsize = (9,7))
plt.plot(N, meanr_array, '*', label = 'Mean r')
plt.plot(x, 1/ysecond, '-k',label = 'Second order inverse fit for mean')
plt.plot(N, stdr_array, 'o', label = 'std(r)')
plt.plot(x, 1/ysecond2, label = 'Second order inverse fit for std')
plt.legend(loc = 'best')
plt.ylabel('[light years]', fontsize = 18)
plt.xlabel('Number of particles', fontsize = 18)
plt.title('Mean radius and standard deviation as a function of N', fontsize = 18)
plt.savefig('../figures/taskf/meanstd.png')
plt.show()
'''
