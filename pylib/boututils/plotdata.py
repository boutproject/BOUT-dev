# Plot a data set

try:
    import numpy as np
    import matplotlib
    import matplotlib.cm as cm
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
except ImportError:
    print "ERROR: plotdata needs numpy and matplotlib to work"
    raise

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

def plotdata(data, x=None, y=None,
             title=None, xtitle=None, ytitle=None,
             output=None, range=None,
             fill=True, mono=False, colorbar=True,
             xerr=None, yerr=None):
    """Plot 1D or 2D data, with a variety of options."""
    
    size = data.shape
    ndims = len(size)
    
    if ndims == 1:
        if (xerr != None) or (yerr != None):
            # Points with error bars
            if x == None:
                x = np.arange(size)
            errorbar(x, data, xerr, yerr)
        # Line plot
        if x == None:
            plt.plot(data)
        else:
            plt.plot(x, data)

    elif ndims == 2:
        # A contour plot
        
        if x == None:
            x = np.arange(size[1])
        if y == None:
            y = np.arange(size[0])
        
        if fill:
            #plt.contourf(data, colors=colors)
            cmap=None
            if mono: cmap = cm.gray
            plt.imshow(data, interpolation='bilinear', cmap=cmap)
        else:
            colors = None
            if mono: colors = 'k'
            
            plt.contour(x, y, data, colors=colors)

        # Add a color bar
        if colorbar:
            CB = plt.colorbar(shrink=0.8, extend='both')
        
    else:
        print "Sorry, can't handle %d-D variables" % ndims
        return
    
    if title != None:
        plt.title(title)
    if xtitle != None:
        plt.xlabel(xtitle)
    if ytitle != None:
        plt.ylabel(ytitle)
    
    if output != None:
        # Write to a file
        plt.savefig(output)
    else:
        # Plot to screen
        plt.show()

def test():
    """Test the plotdata routine."""
    # Generate and plot test data
    
    delta = 0.025
    x = np.arange(-3.0, 3.0, delta)
    y = np.arange(-2.0, 2.0, delta)
    X, Y = np.meshgrid(x, y)
    Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
    Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
    # difference of Gaussians
    Z = 10.0 * (Z2 - Z1)
    
    plotdata(Z, title="test data", fill=False, mono=False)
    plotdata(Z, title="Fill in mono", fill=True, mono=True)
