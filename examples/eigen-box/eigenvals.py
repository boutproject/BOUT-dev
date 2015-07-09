
from boutdata import collect

import matplotlib.pyplot as plt

from numpy import argmin, amax, amin

def plot_eigenvals(eigs):
  fig = plt.figure()
  ax = fig.add_subplot(111)

  eigs_r = eigs[:-1:2]
  eigs_i = eigs[1::2]
  
  range_r = amax(eigs_r) - amin(eigs_r)
  range_i = amax(eigs_i) - amin(eigs_i)
  

  ax.plot(eigs_r, eigs_i, 'x')

  def onclick(event):
    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
      event.button, event.x, event.y, event.xdata, event.ydata)
    
    # Find closest data point, but stretch axes so 
    # real and imaginary components are weighted equally
    
    dist = ((eigs_r - event.xdata)/range_r)**2 + ((eigs_i - event.ydata)/range_i)**2
    
    ind = argmin(dist)
    
    print("Eigenvalue number: %d (%e,%e)" % (ind, eigs_r[ind], eigs_i[ind]))
    
  
  cid = fig.canvas.mpl_connect('button_press_event', onclick)
  plt.show()
  


if __name__ == "__main__":
  path = "data"
  eigs = collect("t_array", path=path)
  plot_eigenvals(eigs)
