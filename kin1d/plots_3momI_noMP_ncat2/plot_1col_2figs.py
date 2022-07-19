#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, ticker,cm

#--------------------------------------------------------------------------------------------------
def plot_column(col):

    zet  = np.zeros((nt,nk), dtype='d')
    qit  = np.zeros((nt,nk), dtype='d')
    Frim = np.zeros((nt,nk), dtype='d')
    rhoi = np.zeros((nt,nk), dtype='d')
    Di   = np.zeros((nt,nk), dtype='d')
    Dhmax= np.zeros((nt,nk), dtype='d')
    muiB = np.zeros((nt,nk), dtype='d')   # column 13  (interporated from lookupTable)
    qc   = np.zeros((nt,nk), dtype='d')
    qr   = np.zeros((nt,nk), dtype='d')
    Wup  = np.zeros((nt,nk), dtype='d')
    Rliq = np.zeros((nt,nk), dtype='d')   #2D array for convenience only (only value at k=nk is relevant)
    Rsol = np.zeros((nt,nk), dtype='d')   #2D array for convenience only (only value at k=nk is relevant)
    Fliq = np.zeros((nt,nk), dtype='d')
    Nr = np.zeros((nt,nk), dtype='d')
    Nc = np.zeros((nt,nk), dtype='d')
    Nit = np.zeros((nt,nk), dtype='d')
    Dr   = np.zeros((nt,nk), dtype='d')

# extract from 'data'.  Last index is the column, correponding to order in 'out_p3.dat'
    for t in range(nt):
        Wup[t,:]  = data[col][(t)*nk:(t+1)*nk, 1]
        zet[t,:]  = data[col][(t)*nk:(t+1)*nk, 2]
        qc[t,:]   = data[col][(t)*nk:(t+1)*nk, 3] *1.e+3
        qr[t,:]   = data[col][(t)*nk:(t+1)*nk, 4] *1.e+3
        qit[t,:]  = data[col][(t)*nk:(t+1)*nk, 5] *1.e+3
        muiB[t,:] = data[col][(t)*nk:(t+1)*nk, 6]
        rhoi[t,:] = data[col][(t)*nk:(t+1)*nk, 7]
        Di[t,:]   = data[col][(t)*nk:(t+1)*nk, 8] *1.e+2   #convert to [cm]
        Dhmax[t,:]= data[col][(t)*nk:(t+1)*nk, 9] *1.e+3   #convert to [mm]
        Rliq[t,:] = data[col][(t)*nk:(t+1)*nk, 10]
        Rsol[t,:] = data[col][(t)*nk:(t+1)*nk, 11]
        Nc[t,:]   = data[col][(t)*nk:(t+1)*nk, 12] 
        Nr[t,:]   = data[col][(t)*nk:(t+1)*nk, 13] 
        Nit[t,:]  = data[col][(t)*nk:(t+1)*nk, 14]        
        Frim[t,:] = data[col][(t)*nk:(t+1)*nk, 15]
        Fliq[t,:] = data[col][(t)*nk:(t+1)*nk, 16]
        Dr[t,:]   = data[col][(t)*nk:(t+1)*nk, 17] *1.e+3 #convert to [mm]

    # apply mask:
    mask = qit > 0.001
    Wup  = np.where(Wup>0., Wup, -1.)
    Rliq = np.where(Rliq>0., Rliq, -1.)
    Rsol = np.where(Rsol>0., Rsol, -1.)
    muiB = np.where(muiB>0., muiB, 0.001)
    muiB = np.where(mask, muiB, -99.)
    #zet  = np.where(mask, zet,  -99.)
    Frim = np.where(mask, Frim, -99.)
    Di   = np.where(mask, Di,   -99.)
    Dhmax= np.where(mask, Dhmax,-99.)
    Fliq = np.where(mask, Fliq, -99.)

    cb_pad = 0.03
    lnwid  = 1.
    lncol  = 'black'

    im0 = ax[0].contourf(range(nt), z_agl, np.transpose(Wup), levels=[0.,1.,2.,3.,4.,5.], vmin=0., vmax=8., cmap='Greys')
    fig.colorbar(im0, ax=ax[0], shrink=0.8, pad=cb_pad)
    ax_b = ax[0].twinx()
    ax_b.plot(range(nt), Rliq, color='red')
    ax_b.plot(range(nt), Rsol, color='blue')
    ax_b.set_ylim([0.,80.])
    ax_b.tick_params(axis="y",direction="in", pad=-22)

    im1 = ax[1].contourf(range(nt), z_agl, np.transpose(qc), levels=levs_Qc,   cmap='Greens')
    fig.colorbar(im1, ax=ax[1], shrink=0.8, pad=cb_pad)
    ax[1].contour( range(nt), z_agl, np.transpose(qc), levels=[0.01,0.01001], linewidths=lnwid, linestyles='dotted', colors='Green')
    ax[1].contour(range(nt), z_agl, np.transpose(qr), levels=levs_Qc,   cmap='Reds')
    ax[1].contour( range(nt), z_agl, np.transpose(qr), levels=[0.01,0.01001], linewidths=lnwid, linestyles='dashed', colors='Red')

    im2 = ax[2].contourf(range(nt), z_agl, np.transpose(qit), levels=levs_Q,   cmap='Blues') #
    fig.colorbar(im2, ax=ax[2], shrink=0.8, pad=cb_pad)
    ax[2].contour( range(nt), z_agl, np.transpose(qit), levels=[0.001,0.001001], linewidths=lnwid, linestyles='dashed', colors=lncol)

    im3 = ax[3].contourf(range(nt), z_agl, np.transpose(Di), levels=levs_Di, colors=ccol_zet)
    fig.colorbar(im3, ax=ax[3], shrink=0.8, pad=cb_pad)
    ax[3].contour( range(nt), z_agl, np.transpose(qit), levels=[0.001,0.001001], linewidths=lnwid, linestyles='dashed', colors=lncol)

    im4 = ax[4].contourf(range(nt), z_agl, np.transpose(zet), levels=levs_Z, colors=ccol_zet)
    fig.colorbar(im4, ax=ax[4], shrink=0.8, pad=cb_pad)


    xloc = 93
    yloc = 12.5
    labsize = 12
    ax[0].text(xloc,yloc,'m s$^{-1}$', size=labsize)
    ax[1].text(xloc,yloc,'g kg$^{-1}$', size=labsize)
    ax[2].text(xloc,yloc,'g kg$^{-1}$', size=labsize)
    ax[3].text(xloc,yloc,'cm',     size=labsize)
    ax[4].text(xloc,yloc,'dBZ',    size=labsize)

    xloc = 5
    yloc = 10.5
    labsize = 18
    ax[0].text(xloc,yloc,'$w$', size=labsize)
    ax[0].text(67,yloc,'$R_{liq}$', size=labsize, color='red')
    ax[0].text(67,yloc-1.8,'$R_{sol}$', size=labsize, color='blue')
    ax[0].text(64,yloc-3.6,'(mm h$^{-1}$)', size=12)
    ax[1].text(xloc,yloc,'$Q_c,Q_r$', size=labsize)
    ax[2].text(xloc,yloc,'$Q_{i,tot}$', size=labsize)
    ax[3].text(xloc,yloc,'$D_m$',       size=labsize)
    ax[4].text(xloc,yloc,'$Z_e$',       size=labsize)


    xaxislabel = 'Time (min)'
    yaxislabel = 'Height (km)'
    if col==0:
        ax[0].yaxis.set_label_text(yaxislabel, size=14)
        ax[1].yaxis.set_label_text(yaxislabel, size=14)
        ax[2].yaxis.set_label_text(yaxislabel, size=14)
        ax[3].yaxis.set_label_text(yaxislabel, size=14)
        ax[4].yaxis.set_label_text(yaxislabel, size=14)

    bottom_ax = ax.shape[0]-1
    ax[bottom_ax].xaxis.set_label_text(xaxislabel, size=14)
    ax[0].set_title(plotTitle[col], size=16, fontweight='bold')

#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
def plot_column2(col):

    zet  = np.zeros((nt,nk), dtype='d')
    qit  = np.zeros((nt,nk), dtype='d')
    Frim = np.zeros((nt,nk), dtype='d')
    rhoi = np.zeros((nt,nk), dtype='d')
    Di   = np.zeros((nt,nk), dtype='d')
    Dhmax= np.zeros((nt,nk), dtype='d')
    muiB = np.zeros((nt,nk), dtype='d')   # column 13  (interporated from lookupTable)
    qc   = np.zeros((nt,nk), dtype='d')
    qr   = np.zeros((nt,nk), dtype='d')
    Wup  = np.zeros((nt,nk), dtype='d')
    Rliq = np.zeros((nt,nk), dtype='d')   #2D array for convenience only (only value at k=nk is relevant)
    Rsol = np.zeros((nt,nk), dtype='d')   #2D array for convenience only (only value at k=nk is relevant)
    Fliq = np.zeros((nt,nk), dtype='d')
    Nr = np.zeros((nt,nk), dtype='d')
    Nc = np.zeros((nt,nk), dtype='d')
    Nit = np.zeros((nt,nk), dtype='d')
    Dr   = np.zeros((nt,nk), dtype='d')
    TT   = np.zeros((nt,nk), dtype='d')

# extract from 'data'.  Last index is the column, correponding to order in 'out_p3.dat'
    for t in range(nt):
        Wup[t,:]  = data[col][(t)*nk:(t+1)*nk, 1]
        zet[t,:]  = data[col][(t)*nk:(t+1)*nk, 2]
        qc[t,:]   = data[col][(t)*nk:(t+1)*nk, 3] *1.e+3
        qr[t,:]   = data[col][(t)*nk:(t+1)*nk, 4] *1.e+3
        qit[t,:]  = data[col][(t)*nk:(t+1)*nk, 5] *1.e+3
        muiB[t,:] = data[col][(t)*nk:(t+1)*nk, 6]
        rhoi[t,:] = data[col][(t)*nk:(t+1)*nk, 7]
        Di[t,:]   = data[col][(t)*nk:(t+1)*nk, 8] *1.e+2   #convert to [cm]
        Dhmax[t,:]= data[col][(t)*nk:(t+1)*nk, 9] *1.e+3   #convert to [mm]
        Rliq[t,:] = data[col][(t)*nk:(t+1)*nk, 10]
        Rsol[t,:] = data[col][(t)*nk:(t+1)*nk, 11]
        Nc[t,:]   = data[col][(t)*nk:(t+1)*nk, 12] 
        Nr[t,:]   = data[col][(t)*nk:(t+1)*nk, 13] 
        Nit[t,:]  = data[col][(t)*nk:(t+1)*nk, 14]        
        Frim[t,:] = data[col][(t)*nk:(t+1)*nk, 15]
        Fliq[t,:] = data[col][(t)*nk:(t+1)*nk, 16]
        Dr[t,:]   = data[col][(t)*nk:(t+1)*nk, 17] *1.e+3 #convert to [mm]
        TT[t,:] = data[col][(t)*nk:(t+1)*nk, 18]

    # apply mask:
    mask = qit > 0.001
    Wup  = np.where(Wup>0., Wup, -1.)
    Rliq = np.where(Rliq>0., Rliq, -1.)
    Rsol = np.where(Rsol>0., Rsol, -1.)
    muiB = np.where(muiB>0., muiB, 0.001)
    muiB = np.where(mask, muiB, -99.)
    #zet  = np.where(mask, zet,  -99.)
    Frim = np.where(mask, Frim, -99.)
    Di   = np.where(mask, Di,   -99.)
    Dhmax= np.where(mask, Dhmax,-99.)
    Fliq = np.where(mask, Fliq, -99.)

    cb_pad = 0.03
    lnwid  = 1.
    lncol  = 'black'

    im0 = ax[0].contourf(range(nt), z_agl, np.transpose(Frim), levels=levs_Fr, colors=ccol_Fr)
    ax[0].contour(range(nt), z_agl, np.transpose(TT), levels=[0.],colors='k')
    fig.colorbar(im0, ax=ax[0], shrink=0.8, pad=cb_pad)

    im1 = ax[1].contourf(range(nt), z_agl, np.transpose(Nc), levels=levs_Nc,   cmap='Greens')
    fig.colorbar(im1, ax=ax[1], shrink=0.8, pad=cb_pad)
    ax[1].contour( range(nt), z_agl, np.transpose(Nc), levels=[1,1.001], linewidths=lnwid, linestyles='dotted', colors='Green')
    ax[1].contour(range(nt), z_agl, np.transpose(Nr), levels=levs_Nr,   cmap='Reds')
    ax[1].contour( range(nt), z_agl, np.transpose(Nr), levels=[1,1.001], linewidths=lnwid, linestyles='dashed', colors='Red')
    ax[1].contour(range(nt), z_agl, np.transpose(TT), levels=[0.],colors='k')

    im2 = ax[2].contourf(range(nt), z_agl, np.transpose(Nit), levels=levs_Ni,   cmap='Blues') #
    fig.colorbar(im2, ax=ax[2], shrink=0.8, pad=cb_pad)
    ax[2].contour( range(nt), z_agl, np.transpose(Nit), levels=[1,1.001], linewidths=lnwid, linestyles='dashed', colors=lncol)

    im3 = ax[3].contourf(range(nt), z_agl, np.transpose(Dr), levels=levs_Di, colors=ccol_zet)
    fig.colorbar(im3, ax=ax[3], shrink=0.8, pad=cb_pad)
    ax[3].contour( range(nt), z_agl, np.transpose(qr), levels=[0.001,0.001001], linewidths=lnwid, linestyles='dashed', colors=lncol)

    im4 = ax[4].contourf(range(nt), z_agl, np.transpose(Fliq), levels=levs_Fr, colors=ccol_Fr)
    ax[4].contour(range(nt), z_agl, np.transpose(TT), levels=[0.],colors='k')
    fig.colorbar(im4, ax=ax[4], shrink=0.8, pad=cb_pad)


    xloc = 93
    yloc = 12.5
    labsize = 12
    ax[0].text(xloc,yloc,'', size=labsize)
    ax[1].text(xloc,yloc,'# kg$^{-1}$', size=labsize)
    ax[2].text(xloc,yloc,'# kg$^{-1}$', size=labsize)
    ax[3].text(xloc,yloc,'mm',     size=labsize)
    ax[4].text(xloc,yloc,'',    size=labsize)

    xloc = 5
    yloc = 10.5
    labsize = 18
    ax[0].text(xloc,yloc,'$F_{i,rim}$', size=labsize)
    ax[1].text(xloc,yloc,'$N_c,N_r$', size=labsize)
    ax[2].text(xloc,yloc,'$N_{i,tot}$', size=labsize)
    ax[3].text(xloc,yloc,'$Dr_m$',       size=labsize)
    ax[4].text(xloc,yloc,'$F_{i,liq}$',       size=labsize)


    xaxislabel = 'Time (min)'
    yaxislabel = 'Height (km)'
    if col==0:
        ax[0].yaxis.set_label_text(yaxislabel, size=14)
        ax[1].yaxis.set_label_text(yaxislabel, size=14)
        ax[2].yaxis.set_label_text(yaxislabel, size=14)
        ax[3].yaxis.set_label_text(yaxislabel, size=14)
        ax[4].yaxis.set_label_text(yaxislabel, size=14)

    bottom_ax = ax.shape[0]-1
    ax[bottom_ax].xaxis.set_label_text(xaxislabel, size=14)
    ax[0].set_title(plotTitle[col], size=16, fontweight='bold')

#--------------------------------------------------------------------------------------------------


#-------- read in data from cld1d; store as individual field arrays:
outputDir = './'
path1 = './'

plotTitle0 = 'P3 v4'
data0 = np.loadtxt(path1 + 'out_p3.dat')
data = [data0]
plotTitle = [plotTitle0]


nk = 41                         # number of vertical levels
nt = int(data[0].shape[0]/nk)   # number of output times
#z_agl = data[0][0:nk,0]        # array of elevations AGL [m]  (column 0)
z_agl = data[0][0:nk,0]/1000.   # array of elevations AGL [km]  (column 0)


#--------  plotting intervals:
#levs_Z    = [-20., 0., 10., 20., 30., 40., 50., 60., 70., 80.]
levs_Z = [-30,-10,0,5,10,15,20,25,30,35,40,45,50,55,60,65,70]
#levs_Q    = [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1., 2., 3., 4., 5.,7.,10.]
levs_Q    = [0.1, 0.2, 0.3, 0.4, 0.5, 1., 2., 3., 4., 5.,6.,7.,8.,9.,10.]
levs_Qc   = [0.1, 0.2, 0.3, 0.4, 0.5, 1., 2., 3., 4., 5.]
levs_Fr   = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
levs_Vi   = [0.001, 0.2, 0.5, 1., 2., 3., 4., 5., 7., 10., 15., 20., 25.]
#levs_muB  = [0., 0.5, 1., 2., 3., 4., 5., 10., 15., 20., 25.]
levs_Di   = [0., 0.1,0.2,0.3, 0.4, 0.5, 1., 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0,6.0, 7.0, 8.0, 9.0]
levs_Dh   = levs_Di
levs_rhoi = [0.001, 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.]
levs_lami = [0.000, 25., 50., 75., 100., 200., 300., 400., 500.]
#levs_mui = [0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.]
levs_Nr    = [0,20000,40000,60000,80000,100000,120000,140000,160000,180000,200000]
levs_Nc    = [0,0.4e8,0.8e8,1.2e8,1.6e8,2.0e8,2.4e8,2.8e8,3.0e8,5e8]
levs_Ni    = [0,50000,75000,100000,250000,500000,750000,1000000,2500000,5000000,7500000]

ccol_zet  = ['lightgrey','darkgrey','grey','lightskyblue','dodgerblue','blue','limegreen','forestgreen',
	     'darkgreen','yellow','orange','salmon','red','darkred','darkviolet','Indigo','darkblue','black']
levs_mui = [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,99.]

ccol_Fr = ['lightgrey','lightskyblue','blue','limegreen','forestgreen',
	     'yellow','orange','salmon','red','darkred']
fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(9,14), sharex=True)

plot_column(0)

# adjust subplots
#plt.subplots_adjust(hspace=0.10)
#plt.subplots_adjust(wspace=0.01)

plt.tight_layout()

#--------------------------
plt.savefig('./fig1.png', format='png')
#plt.show()

fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(9,14), sharex=True)

plot_column2(0)

# adjust subplots
#plt.subplots_adjust(hspace=0.10)
#plt.subplots_adjust(wspace=0.01)

plt.tight_layout()

#--------------------------
plt.savefig('./fig2.png', format='png')
#plt.show()
