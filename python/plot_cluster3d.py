from plot_cluster import *
from mpl_toolkits.mplot3d import Axes3D

filename = '../benchmarks/pos_N152_dt0.010000_tcoll_1.000000_eps0.000000.xyz'
pos, mass, potential, kinetic = read_file(filename)

fig3d = figure()
ax3d = fig3d.add_subplot(111,projection='3d')

max_time = len(pos[0][:,0])
for time in range(max_time):
    ax3d.cla()
    for i in range(len(pos)):
        ax3d.scatter(pos[i][time,0],pos[i][time,1],pos[i][time,2])

    ax3d.set_xlabel('x [ly]',size=15)
    ax3d.set_ylabel('y [ly]',size=15)
    ax3d.set_zlabel('z [ly]',size=15)
    ax3d.set_title('Open cluster, t = %.1f t_coll'%(time/float(max_time)*10),size=15)
    ax3d.set_xlim([-22,22])
    ax3d.set_ylim([-22,22])
    ax3d.set_zlim([-22,22])

    fig3d.savefig('../figures/Eirik_plot_%3d.png'%time)

#import subprocess
#subprocess.call('convert -delay 8 -loop 0 ../figures/Eirik_plot_*.png ../figures/anim3d_long.gif', shell=True)
