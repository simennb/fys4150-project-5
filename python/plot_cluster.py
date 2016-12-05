from pylab import *

def read_file(filename):
    data = open(filename,'r')
    totlines = 0
    for line in data:
        totlines += 1
    data.close()
    data = open(filename,'r')
    # I am most definitely overcomplicating this
    i = 0
    j = 0
    k = 0
    m = 0
    for line in data:
        if i == 0:
            num_bodies = int(line)
            positions = [zeros((totlines/(num_bodies+2),3)) for a in range(num_bodies)]
            j += 1
        elif i == j*(num_bodies+2):
            k = 0
            j += 1
        elif not line.startswith('Comment'):
            a = line.split(' ')

            positions[k][int(m)] = (array([float(a[1]),float(a[2]),float(a[3])]))
            k += 1
            m += 1./num_bodies
        i += 1
    data.close()
    return positions

if __name__=='__main__':
    filename = '../benchmarks/pos_N100_dt0.001000_tcoll_5.000000_eps0.200000.xyz'
    pos = read_file(filename)
    figure()
    max_time = len(pos[0][:,0])
    for time in range(max_time):
        clf()
        for i in range(len(pos)):
            plot(pos[i][time,0],pos[i][time,1], '*b',ms=12)
        xlim([-22,22])
        ylim([-22,22])
        xlabel('x [ly]',size=15)
        ylabel('y [ly]',size=15)
        title('Open cluster, t = %.1f t_coll'%(time/float(max_time)*5),size=15)
        grid('on')
        savefig('../figures/plot_200eps02_%3d.png'%time)
#    show()

    import os
    import subprocess
    subprocess.call('convert -delay 8 -loop 0 ../figures/plot_200eps02_*.png ../figures/anim.gif', shell=True)
