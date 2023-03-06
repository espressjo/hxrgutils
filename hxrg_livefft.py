"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""
#basically export PYTHONPATH=$PYTHONPATH:/vltuser/nirps/module
import sys
sys.path.append('/vltuser/nirps/module')
from scanf import scanf
from os.path import basename
from ff import ff
from matplotlib import pyplot as plt
from matplotlib import animation
from hxrg.type.read import hxread
#to fake the readout of file.
p = '/data/NIRPS/INS_ROOT/SYSTEM/DETDATA'
fmt = "NIRPS_R%2.2d_R%2.2d.fits"
ii = 0
# First set up the figure, the axis, and the plot element we want to animate
plt.rcParams['axes.facecolor'] = '#242925'

f,ax = plt.subplots(figsize=(11,6),facecolor=(.31,.31,.31))
ax.set_xlim((0,200))
ax.set_ylim((0,40))
ax.set_xlabel('f (Hz)',fontsize=20)
ax.set_ylabel('Amplitude',fontsize=20)

ax.xaxis.label.set_color('red')
ax.yaxis.label.set_color('red')
#ttl = ax.text(10, 10, '')#ax = plt.axes(xlim=(0, 300), ylim=(0, 40))
line, = ax.plot([], [], lw=1,color='#87CEEB')
time_text = ax.text(60, 38, '')
time_text.set_color('green')
time_text.set_fontsize(22)
ax.text(60, 25, '')

ff = ff(p)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line,
# animation function.  This is called sequentially
def animate(i):
    
    print(i)
    #f1,f2 = ff.find_next()
    f1,f2 = ff.find_last()
    r1 = hxread()
    r1(f1)
    r2 = hxread()
    r2(f2)
    cds = r2-r1
    x,y = cds.fft()
    d = f1.split('/')[-2]
    ra = scanf("NIRPS_R%d_R01.fits",basename(f1))[0]
    time_text.set_text("Path: %s, ramp: %d"%(d,ra))
    
    line.set_data(x, y)
    return tuple([line]) + tuple([time_text])



# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(f, animate, init_func=init,
                               frames=None, interval=20, blit=True)

plt.show()