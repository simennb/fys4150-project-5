import numpy as np
import matplotlib.pyplot as plt
from plot_cluster import read_file
import sys

def energyfinder(N, dt, tcoll, eps, eqtime):
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

    x = np.linspace(0, tcoll, 101)

    print len(x)
    print len(total)
    print len(total[i])

    b = 0

    tot_EnergyPerTime = np.zeros(len(total[1]))
    tot_kinEnergyPerTime = np.zeros(len(total[1]))
    tot_potEnergyPerTime = np.zeros(len(total[1]))


    for j in range(0, len(total[0])):   #time
        for i in range(0, len(total)):  #particle
        #if j == 0:
            #print total[i][j]
            #b += total[i][j]
            tot_EnergyPerTime[j] += total[i][j]
            tot_kinEnergyPerTime[j] += kinetic[i][j]
            tot_potEnergyPerTime[j] += potential[i][j]
    print
    print tot_EnergyPerTime[0]  #At time zero all the particles combined has this energy
    print b

    counterarray = np.zeros(len(x))
    kinP = np.zeros(len(x))
    kinTot = np.zeros(len(x))
    potTot = np.zeros(len(x))

    for i in range(0, len(x)):   #time
        for j in range(0, len(total)):  #particle
            if total[j][i] > 0:
                kinP[i] += kinetic[j][i]
                counterarray[i] += 1
            kinTot[i] += kinetic[j][i]

        
    kinN = kinTot - kinP


    ######################################
    ####   The part that is from m_   ####
    ######################################

    filenameEnergy = "../benchmarks/energy_N%d_dt%1.6f_tcoll_%1.6f_eps%1.6f.xyz" %(N, dt, tcoll, eps)
    data = np.loadtxt(filenameEnergy)
    m_kinEnergy = data[:, 0]
    m_potEnergy = data[:, 1]
    m_totEnergy = m_kinEnergy + m_potEnergy

    viral = 2*kinN/float(N) + m_potEnergy/float(N)

    return m_totEnergy, x, counterarray, viral, kinP, kinTot, kinN, pos
    
#plt.plot(x, kinTot/abs(m_potEnergy[0]), label = 'Total')
#plt.plot(x, kinN/abs(potTot[0]), label = 'Bound')
#plt.plot(x, kinP/abs(m_potEnergy[0]), label = 'Unbound')
#plt.legend()
#plt.plot(x, counterarray/N)
#plt.plot(x, viral)
#plt.ylim([-1, 1])
#plt.show()

######################################
####            Task b            ####
######################################
if sys.argv[1] == "b":
    N = [1000]
    dt = 0.001
    tcoll = 5.0
    eps = 0
    eqtime = 0

    for i in range(0, len(N)):
        m_totEnergy, x, counterarray, viral, kinP, kinTot, kinN, pos = energyfinder(N[i], dt, tcoll, eps, eqtime)
        plt.figure()
        max_time = len(pos[0][:,0])
        for time in range(max_time):
            plt.clf()
            for j in range(len(pos)):
                plt.plot(pos[j][time,0],pos[j][time,1], '*b',ms=12)
            plt.xlim([-22,22])
            plt.ylim([-22,22])
            plt.xlabel('x [ly]',size=18)
            plt.ylabel('y [ly]',size=18)
            plt.title('Open cluster at %1.2f collision times'%(time/float(max_time)*tcoll),size=18)
            plt.grid('on')
            plt.savefig('../figures/taskb/N%d/plot_N%d_time%1.2f.png'%(N[i], N[i], time))
            #plt.show()

######################################
####            Task c            ####
######################################
if sys.argv[1] == "c":
    N = [100, 200, 300, 400, 500, 1000]
    dt = 0.001
    tcoll = 5.0
    eps = 0
    eqtime = 0

    for i in range(0, len(N)):
        m_totEnergy, x, counterarray, viral, kinP, kinTot, kinN, pos = energyfinder(N[i], dt, tcoll, eps, eqtime)
        plt.figure(figsize = (9,10))
        plt.subplot(3,1,1)
        plt.plot(x, m_totEnergy)
        plt.title(r'%d particles' %(N[i]), fontsize = 18)
        #plt.xlabel("Time", fontsize = 18)
        plt.ylabel("Energy", fontsize = 18)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        plt.subplot(3,1,2)
        plt.plot(x, counterarray/N[i], label = 'Fraction of particles with positive total energy')
        #plt.title('Fraction of particles with positive total energy', fontsize = 18)
        plt.ylabel(r'f$^{p}$', fontsize = 24)
        #plt.xlabel('Time', fontsize = 18)
        #plt.legend(loc = 'best')

        plt.subplot(3,1,3)
        plt.plot(x, kinP, label = 'Kinetic energy of escaped particles')
        plt.ylabel(r'K$^{P}$', fontsize = 20)
        plt.xlabel('Time', fontsize = 18)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.savefig("../figures/taskc/totenergy_N%d.png" %(N[i]))
        plt.show()


######################################
####            Task d            ####
######################################

if sys.argv[1] == "d1":
    N = [200, 200, 200, 200, 200, 200]
    dt = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    tcoll = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    eps = [5, 2, 0.5, 0.2, 0.05, 0.02]
    eqtime = [75, 75, 75, 75, 75, 75]

    for i in range(0, len(N)):
        m_totEnergy, x, counterarray, viral, kinP, kinTot, kinN, pos = energyfinder(N[i], dt[i], tcoll[i], eps[i], eqtime[i])
        plt.figure(figsize = (9, 7))
        plt.subplot(2,1,1)
        plt.plot(x, m_totEnergy)
        plt.title(r'Total Energy with $\epsilon = %1.2e$ and %d particles' %(eps[i], N[i]), fontsize = 18)
        #plt.xlabel("Time", fontsize = 18)
        plt.ylabel("Energy", fontsize = 18)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        #plt.savefig("../figures/taskd/totenergy_eps%1.2e_N%d.png" %(eps[i], N[i]))
        #plt.show()

        plt.subplot(2,1,2)
        plt.plot(x, counterarray/N[i], label = 'Fraction of particles with positive total energy')
        #plt.title('Fraction of particles with positive total energy', fontsize = 18)
        plt.ylabel('High energy fraction', fontsize = 18)
        plt.xlabel('Time', fontsize = 18)
        #plt.legend(loc = 'best')
        plt.savefig("../figures/taskd/totenergy_eps%1.2e_N%d.png" %(eps[i], N[i]))
        #plt.show()

if sys.argv[1] == "d2":
    N = [100, 200, 300, 400, 500, 600, 1000]
    dt = 0.001
    tcoll = 5.0
    eps = 0.1
    eqtime = 75

    for i in range(0, len(N)):
        m_totEnergy, x, counterarray, viral, kinP, kinTot, kinN, pos = energyfinder(N[i], dt, tcoll, eps, eqtime)
        plt.figure(figsize = (9,10))
        plt.subplot(3,1,1)
        plt.plot(x, m_totEnergy)
        plt.title(r'%d particles' %(N[i]), fontsize = 18)
        #plt.xlabel("Time", fontsize = 18)
        plt.ylabel("Energy", fontsize = 18)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        plt.subplot(3,1,2)
        plt.plot(x, counterarray/N[i], label = 'Fraction of particles with positive total energy')
        #plt.title('Fraction of particles with positive total energy', fontsize = 18)
        plt.ylabel(r'f$^{p}$', fontsize = 24)
        #plt.xlabel('Time', fontsize = 18)
        #plt.legend(loc = 'best')

        plt.subplot(3,1,3)
        plt.plot(x, kinP, label = 'Kinetic energy of escaped particles')
        plt.ylabel(r'K$^{P}$', fontsize = 20)
        plt.xlabel('Time', fontsize = 18)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.savefig("../figures/taskd/totenergy_N%d.png" %(N[i]))
        plt.show()
        
######################################
####            Task e            ####
######################################

if sys.argv[1] == "e":
    N = [100, 200, 300, 400, 500, 600, 1000]
    dt = 0.001
    tcoll = 5.0
    eps = 0.1
    eqtime = 75

    for i in range(0, len(N)):
        m_totEnergy, x, counterarray, viral, kinP, kinTot, kinN, pos = energyfinder(N[i], dt, tcoll, eps, eqtime)
        plt.figure(figsize = (9, 7))
        plt.plot(x, viral, label = '2*<K> - <V>')
        plt.plot(x, np.zeros(len(x)), '-k')
        plt.title(r'Viral theorem test for %d particles' %(N[i]), fontsize = 20)
        plt.xlabel("Time", fontsize = 20)
        plt.ylabel("Energy", fontsize = 20)
        plt.legend(loc = 'best')
        #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.savefig("../figures/taske/viral_N%d.png" %(N[i]))
        plt.show()
