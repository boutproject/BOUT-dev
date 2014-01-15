from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np


def SURFACE(Z, fig, xtitle=None, ytitle=None, title=None, var=None, sub=None):

    

    if sub==None :
        ax = fig.gca(projection='3d')
    else:
        ax = fig.add_subplot(sub[0],sub[1],sub[2], projection='3d')
     
    nx=np.shape(Z)[0]
    ny=np.shape(Z)[1]
    
    zmin=np.min(Z)
    zmax=np.max(Z)
     
    X = np.arange(nx)
    Y = np.arange(ny)
    X, Y = np.meshgrid(X, Y)

    
    surf = ax.plot_surface(X.T, Y.T, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
            linewidth=0, antialiased=False)
    ax.set_zlim(zmin-1., zmax+1.)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    ax.set_zticklabels([])
    
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    if title != None : ax.set_zlabel(title)
    
    
    cbar=fig.colorbar(surf, shrink=1., aspect=15, orientation='horizontal', format='%.1e')
    cbar.ax.set_xlabel(var)
    cbar.ax.tick_params(labelsize=10)
   
     

    plt.draw()

