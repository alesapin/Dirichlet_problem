import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.interpolate import griddata
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import math

def read_file(fpath):
    result = []
    with open(fpath,'r') as f:
        for line in f:
            result.append(tuple(float(n) for n in line.split('\t')))
    return result

def func(x,y):
    return math.sqrt(4+x*y)

def plot_3d(x, y, z, xlabel, ylabel, zlabel):
    X = np.asarray(x)
    Y = np.asarray(y)
    Z = np.asarray(z)
    print 'arrayed'
    xi = np.linspace(X.min(),X.max(),100)
    print 'linspaced x'
    yi = np.linspace(Y.min(),Y.max(),100)
    print 'linspaced y'
    zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')
    print 'approximated'
    xig, yig = np.meshgrid(xi, yi)
    print 'meshed'

    surf = ax.plot_surface(xig, yig, zi, cmap=cm.rainbow_r, rstride=1, cstride=1)
    cbar = fig.colorbar(surf, aspect=10)
    cbar.remove()
    ax.set_zlim(zmin=0)
    ax.set_ylim(ymin=0)
    ax.set_xlim(xmin=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    plt.show()

def error(x, y, zreal, zcount):
    error = 0
    mesh_size = int(math.sqrt(len(x)))
    print 'mesh_size:' + str(mesh_size)
    for i in xrange(1, mesh_size-1):
        for j in xrange(1, mesh_size-1):
            val = (zreal[i*mesh_size + j]-zcount[i*mesh_size + j])
            hs1 = ((x[i*mesh_size + j + 1] - x[i*mesh_size + j]) + (x[i*mesh_size + j] - x[i*mesh_size + j - 1]))/2.
            hs2 = ((y[i*mesh_size + j + 1] - y[i*mesh_size + j]) + (y[i*mesh_size + j] - y[i*mesh_size + j - 1]))/2.
            #if hs1 <= 0 or hs2 <= 0:
            #    print i
            error += val*val*hs1*hs2
    print error
    return math.sqrt(error)

if __name__ == "__main__":
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    data = read_file('result_2000x2000_512.txt')
    print 'readed'
    x,y,z = zip(*data)
    real_z = []
    for i in xrange(0, len(x)):
        real_z.append(func(x[i],y[i]))
    print 'error: ' + str(error(x,y,real_z,z))
    #print 'splited'
    #plot_3d(x,y,real_z,'X', 'Y', 'Z')

