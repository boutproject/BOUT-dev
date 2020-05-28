# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as np
from boututils.bunch import Bunch
from .radial_grid import radial_grid
from .analyse_equil_2 import analyse_equil
from pylab import (cm, clabel, contour, draw, legend, plot, setp, show,
                   streamplot, subplot2grid, text, tick_params, figure)
from .ask import query_yes_no
from .read_geqdsk import read_geqdsk
from scipy import interpolate



def View2D(g, option=0):

    # plot and check the field
    fig=figure(num=2, figsize=(16, 6))
   # fig.suptitle('Efit Analysis', fontsize=20)

    ax = subplot2grid((3,3), (0,0), colspan=1, rowspan=3)
    nlev = 100
    minf = np.min(g.psi)
    maxf = np.max(g.psi)
    levels = np.arange(np.float(nlev))*(maxf-minf)/np.float(nlev-1) + minf
    ax.contour(g.r,g.z,g.psi, levels=levels)
#    ax.set_xlim([0,6])
    ax.set_xlabel('R')
    ax.set_ylabel('Z')
    ax.yaxis.label.set_rotation('horizontal')

    ax.set_aspect('equal')

   # fig.suptitle('Efit Analysis', fontsize=20)
   # title('Efit Analysis', fontsize=20)
    text(0.5, 1.08, 'Efit Analysis',
         horizontalalignment='center',
         fontsize=20,
         transform = ax.transAxes)


    draw()
    if option == 0 : show(block=False)

    plot(g.xlim,g.ylim,'g-')

    draw()

    csb=contour( g.r, g.z, g.psi,  levels=[g.sibdry])

    clabel(csb, [g.sibdry],  # label the level
            inline=1,
            fmt='%9.6f',
            fontsize=14)

    csb.collections[0].set_label('boundary')

  # pl1=plot(g.rbdry,g.zbdry,'b-',marker='x', label='$\psi=$'+ np.str(g.sibdry))
  #  legend(bbox_to_anchor=(0., 1.05, 1., .105), loc='upper left')



    draw()

   # fig.set_tight_layout(True)

  #  show(block=False)

# Function fpol and qpsi are given between simagx (psi on the axis) and sibdry (
# psi on limiter or separatrix). So the toroidal field (fpol/R) and the q profile are within these boundaries

    npsigrid=old_div(np.arange(np.size(g.pres)).astype(float),(np.size(g.pres)-1))

    fpsi = np.zeros((2, np.size(g.fpol)), np.float64)
    fpsi[0,:] = (g.simagx + npsigrid * ( g.sibdry -g.simagx ))
    fpsi[1,:] = g.fpol



    boundary = np.array([g.xlim, g.ylim])

    rz_grid = Bunch(nr=g.nx, nz=g.ny,   # Number of grid points
                   r=g.r[:,0], z=g.z[0,:], # R and Z as 1D arrays
                   simagx=g.simagx, sibdry=g.sibdry, # Range of psi
                   psi=g.psi, # Poloidal flux in Weber/rad on grid points
                   npsigrid=npsigrid, # Normalised psi grid for fpol, pres and qpsi
                   fpol=g.fpol, # Poloidal current function on uniform flux grid
                   pres=g.pres, # Plasma pressure in nt/m^2 on uniform flux grid
                   qpsi=g.qpsi, # q values on uniform flux grid
                   nlim=g.nlim, rlim=g.xlim, zlim=g.ylim) # Wall boundary




    critical = analyse_equil(g.psi,g.r[:,0],g.z[0,:])

    n_opoint = critical.n_opoint
    n_xpoint = critical.n_xpoint
    primary_opt = critical.primary_opt
    inner_sep   = critical.inner_sep
    opt_ri = critical.opt_ri
    opt_zi = critical.opt_zi
    opt_f  = critical.opt_f
    xpt_ri = critical.xpt_ri
    xpt_zi = critical.xpt_zi
    xpt_f  = critical.xpt_f


    psi_inner=0.6
    psi_outer=0.8,
    nrad=68
    npol=64
    rad_peaking=[0.0]
    pol_peaking=[0.0]
    parweight=0.0


    boundary = np.array([rz_grid.rlim, rz_grid.zlim])

  # Psi normalisation factors

    faxis = critical.opt_f[critical.primary_opt]

    fnorm = critical.xpt_f[critical.inner_sep] - critical.opt_f[critical.primary_opt]



  # From normalised psi, get range of f
    f_inner = faxis + np.min(psi_inner)*fnorm
    f_outer = faxis + np.max(psi_outer)*fnorm


    fvals = radial_grid(nrad, f_inner, f_outer, 1, 1, [xpt_f[inner_sep]], rad_peaking)


    ## Create a starting surface
    #sind = np.int(nrad / 2)
    #start_f = 0. #fvals[sind]


    # Find where we have rational surfaces
    # define an interpolation of psi(q)

    psiq = np.arange(np.float(g.qpsi.size))*(g.sibdry-g.simagx)/np.float(g.qpsi.size-1) + g.simagx
    fpsiq=interpolate.interp1d(g.qpsi, psiq)

    # Find how many rational surfaces we have within the boundary and locate x,y position of curves

    nmax=g.qpsi.max().astype(int)
    nmin=g.qpsi.min().astype(int)

    nr=np.arange(nmin+1,nmax+1)
    psi=fpsiq(nr)

    cs=contour( g.r, g.z, g.psi,  levels=psi)
    labels = ['$q='+np.str(x)+'$\n' for x in range(nr[0],nr[-1]+1)]

    for i in range(len(labels)):
        cs.collections[i].set_label(labels[i])

    style=['--', ':', '--', ':','-.' ]


#    gca().set_color_cycle(col)

#    proxy = [Rectangle((0,0),1,1, fc=col[:i+1])
#        for pc in cs.collections]
#
#    l2=legend(proxy, textstr[:i)
#
#    gca().add_artist(l1)


    x=[]
    y=[]
    for i in range(psi.size):
        xx,yy=surface(cs, i, psi[i],opt_ri[primary_opt], opt_zi[primary_opt], style, option)
        x.append(xx)
        y.append(yy)


    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #textstr = ['$q='+np.str(x)+'$\n' for x in range(psi.size)]
    #
    #ax.text(0.85, 1.15, ''.join(textstr), transform=ax.transAxes, fontsize=14,
    #verticalalignment='top', bbox=props)


    legend(bbox_to_anchor=(-.8, 1), loc='upper left', borderaxespad=0.)

    draw()

    #compute B - field

    Bp=np.gradient(g.psi)

    dr=old_div((np.max(g.r[:,0])-np.min(g.r[:,0])),np.size(g.r[:,0]))
    dz=old_div((np.max(g.z[0,:])-np.min(g.z[0,:])),np.size(g.z[0,:]))

    dpsidr=Bp[0]
    dpsidz=Bp[1]

    Br=-dpsidz/dz/g.r
    Bz=dpsidr/dr/g.r

    Bprz=np.sqrt(Br*Br+Bz*Bz)

    # plot Bp field
    if option == 0 :
        sm = query_yes_no("Overplot vector field")
        if sm :
            lw = 50*Bprz/Bprz.max()
            streamplot(g.r.T,g.z.T, Br.T,Bz.T, color=Bprz, linewidth=2, cmap=cm.bone)#density =[.5, 1], color='k')#, linewidth=lw)
            draw()


    # plot toroidal field

    ax = subplot2grid((3,3), (0,1), colspan=2, rowspan=1)
    ax.plot(psiq,fpsi[1,:])
    #ax.set_xlim([0,6])
    #ax.set_xlabel('$\psi$')
    ax.set_ylabel('$fpol$')
    ax.yaxis.label.set_size(20)
    #ax.xaxis.label.set_size(20)

    ax.set_xticks([])

    ax.yaxis.label.set_rotation('horizontal')
    #ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 20

    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #ax.text(0.85, 0.95, '$B_t$', transform=ax.transAxes, fontsize=14,
    #    verticalalignment='top', bbox=props)


    #
    draw()

        # plot pressure

    ax = subplot2grid((3,3), (1,1), colspan=2, rowspan=1)
    ax.plot(psiq,g.pres)
    #ax.set_xlim([0,6])
    #ax.set_xlabel('$\psi$')
    ax.set_ylabel('$P$')
    ax.yaxis.label.set_size(20)
    #ax.xaxis.label.set_size(20)

    ax.set_xticks([])

    ax.yaxis.label.set_rotation('horizontal')
    #ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 20

    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #ax.text(0.85, 0.95, '$B_t$', transform=ax.transAxes, fontsize=14,
    #    verticalalignment='top', bbox=props)


    #
    draw()


    # plot qpsi

    ax = subplot2grid((3,3), (2,1), colspan=2, rowspan=1)

    ax.plot(psiq,g.qpsi)
    ax.set_xlabel('$\psi$')
    ax.set_ylabel('$q$')
    ax.yaxis.label.set_rotation('horizontal')
    ax.yaxis.label.set_size(20)
    ax.xaxis.label.set_size(20)

    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 20

    tick_params(\
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on')


    draw()


    # Compute and draw Jpar
#
#    MU = 4.e-7*np.pi
#
#    jpar0 = - Bxy * fprime / MU - Rxy*Btxy * dpdpsi / Bxy
#

   # fig.set_tight_layout(True)


    fig.subplots_adjust(left=0.2, top=0.9, hspace=0.1, wspace=0.5)

    draw()

    if option == 0 : show(block=False)



    if option != 0:
        return Br,Bz, x, y, psi


## output to files
#np.savetxt('../data.in', np.reshape([g.nx, g.ny],(1,2)), fmt='%i, %i')
#f_handle = open('../data.in', 'a')
#np.savetxt(f_handle,np.reshape([g.rmagx, g.zmagx],(1,2)), fmt='%e, %e')
#np.savetxt(f_handle, np.reshape([np.max(g.r),np.max(g.z)],(1,2)))
#np.savetxt(f_handle, np.reshape([g.bcentr,g.rcentr],(1,2)))
#f_handle.close()
#
#
#f_handle = open('../input.field', 'w')
#np.savetxt(f_handle, np.reshape([dx, dy],(1,2)))
#np.savetxt(f_handle, (g.r[:,0],g.z[0,:]))
#np.savetxt(f_handle, Bx)  # fortran compatibility
#np.savetxt(f_handle, By)
#f_handle.close()
#
#np.savetxt('../tbound',[np.size(g.xlim)], fmt='%i')
#f_handle = open('../tbound', 'a')
#np.savetxt(f_handle, (g.xlim,g.ylim))
#f_handle.close()
#
#np.savetxt('../pbound',[np.size(g.rbdry)], fmt='%i')
#f_handle = open('../pbound', 'a')
#np.savetxt(f_handle, (g.rbdry,g.zbdry))
#f_handle.close()

def surface(cs, i, f, opt_ri, opt_zi, style, iplot=0):

  #  contour_lines( F, np.arange(nx).astype(float), np.arange(ny).astype(float), levels=[start_f])
#    cs=contour( g.r, g.z, g.psi,  levels=[f])
#    proxy = [Rectangle((0,0),1,1,fc = 'b')
#        for pc in cs.collections]
#
#    legend(proxy, ["q="+np.str(i)])



    p = cs.collections[i].get_paths()
 #
 #  You might get more than one contours for the same start_f. We need to keep the closed one
    vn=np.zeros(np.size(p))

 # find the closed contour

    for k in range(np.size(p)):
          v=p[k].vertices
          vx=v[:,0]
          vy=v[:,1]
          if [vx[0], vy[0]] == [vx[-1],vy[-1]] :
              xx=vx
              yy=vy


    x=xx
    y=yy
    #v = p[0].vertices
    #vn[0]=np.shape(v)[0]
    #xx=v[:,0]
    #yy=v[:,1]

    #if np.shape(vn)[0] > 1:
    #    for i in xrange(1,np.shape(vn)[0]):
    #        v = p[i].vertices
    #        vn[i]=np.shape(v)[0]
    #        xx = [xx,v[:,0]]
    #        yy = [yy,v[:,1]]

    #if np.shape(vn)[0] > 1 :
    ## Find the surface closest to the o-point
    #    ind = closest_line(np.size(xx), xx, yy, opt_ri, opt_zi)
    #    x=xx[ind]
    #    y=yy[ind]
    #else:
    #    ind = 0
    #    x=xx
    #    y=yy
    #
    if(iplot == 0):

    # plot the start_f line
        zc = cs.collections[i]
        setp(zc, linewidth=4, linestyle=style[i])

        clabel(cs, [f],  # label the level
            inline=1,
            fmt='%9.6f',
            fontsize=14)

    #    annotate('q= '+np.str(i+1),(x[0]+.1,y[0]+.1))

        draw()

        show(block=False)

    return x,y


if __name__ == '__main__':

    path='../../tokamak_grids/pyGridGen/'

    g=read_geqdsk(path+"g118898.03400")

    View2D(g, option=0)
    show()
