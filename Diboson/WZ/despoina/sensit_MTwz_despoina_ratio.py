#!/usr/bin/ python3

# from matplotlib import gridspec
# from  matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text',usetex=True)
# from pylab import *
import yoda
import numpy as np
# from decimal import Decimal
# # from termcolor import colorescolors[2]
# import math
import matplotlib as mpl
# import matplotlib.pyplot as plt
# import itertools
# import re
# import matplotlib.cbook as cbook
# # from ._color_data import BASE_COLORS, TABLEAU_COLORS, CSS4_COLORS, XKCD_COLORS
# import matplotlib.gridspec as gridspec
# import math
import os
# import sys

#parametrization of Wilson coeff
v = 0.246
cwil=10.

print("cW:")
file_1 = yoda.read('Inclusive/Inclusive_yodas_BSM1/cW.yoda')
hist_1 = file_1[ '/TESTDET_BSM1/d12-x01-y01']
edges_1 = hist_1.xEdges()
vals_1cw = hist_1.areas()
print("vals BSM1: ",vals_1cw, "Edges 1: ",edges_1 ,"XS (in fb)", np.sum(vals_1cw))
#
#
file_2 = yoda.read('Inclusive/Inclusive_yodas_BSM2/cW.yoda')
hist_2 = file_2[ '/TESTDET_BSM2/d12-x01-y01']
edges_2 = hist_2.xEdges()
vals_2cw = hist_2.areas()
print("vals BSM2: ",vals_2cw, "Edges 2: ",edges_2 ,"XS (in fb)", np.sum(vals_2cw))
#
#
file_3 = yoda.read('Inclusive/SM/cW.yoda')
hist_3 = file_3[ '/TESTDET/d12-x01-y01']
edges_3 = hist_3.xEdges()
vals_3cw = hist_3.areas()
print("vals SM: ",vals_3cw, "Edges 3: ",edges_3 ,"XS (in fb)", np.sum(vals_3cw))
#
#
file_4 = yoda.read('VBS/yoda_VBS/cW.yoda')
hist_4 = file_4[ '/TESTDET/d12-x01-y01']
edges_4 = hist_4.xEdges()
vals_4cw = hist_4.areas()
print("vals VBS: ",vals_4cw, "Edges 4: ",edges_4 ,"XS (in fb)", np.sum(vals_4cw))
#
#
print("SM:")
file_1 = yoda.read('Inclusive/Inclusive_yodas_BSM1/sm.yoda')
hist_1 = file_1[ '/TESTDET_BSM1/d12-x01-y01']
edges_1 = hist_1.xEdges()
vals_1sm = hist_1.areas()
print("vals BSM1: ",vals_1sm, "Edges 1: ",edges_1 ,"XS (in fb)", np.sum(vals_1sm))
#
#
file_2 = yoda.read('Inclusive/Inclusive_yodas_BSM2/sm.yoda')
hist_2 = file_2[ '/TESTDET_BSM2/d12-x01-y01']
edges_2 = hist_2.xEdges()
vals_2sm = hist_2.areas()
print("vals BSM2: ",vals_2sm, "Edges 2: ",edges_2 ,"XS (in fb)", np.sum(vals_2sm))
#
#
file_3 = yoda.read('Inclusive/SM/sm.yoda')
hist_3 = file_3[ '/TESTDET/d12-x01-y01']
edges_3 = hist_3.xEdges()
vals_3sm = hist_3.areas()
print("vals SM: ",vals_3sm, "Edges 3: ",edges_3 ,"XS (in fb)", np.sum(vals_3sm))
#
#
file_4 = yoda.read('VBS/yoda_VBS/sm.yoda')
hist_4 = file_4[ '/TESTDET/d12-x01-y01']
edges_4 = hist_4.xEdges()
vals_4sm = hist_4.areas()
print("vals VBS: ",vals_4sm, "Edges 4: ",edges_4 ,"XS (in fb)", np.sum(vals_4sm))
#
# Compare the different sets of cuts
#
#

vals_1 = vals_1cw/vals_1sm
vals_1[np.isnan(vals_1)] = 0
vals_2 = vals_2cw/vals_2sm
vals_2[np.isnan(vals_2)] = 0
vals_3 = vals_3cw/vals_3sm
vals_3[np.isnan(vals_3)] = 0
vals_4 = vals_4cw/vals_4sm
vals_4[np.isnan(vals_4)] = 0

print("Ratios: ", vals_1, vals_2, vals_3, vals_4)

#quit()

#################################################################################################
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
rescolors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
#################################################################################################
# Creating data sequence: middle of each bin
edges = edges_1
xBinning = edges
print("xBinning:" , xBinning)
xmin = edges[0]
xmax = edges[-1]
xData = []
binwidth = []
for i in range(len(edges)-1):
    xData = np.append(xData, edges[i]+(edges[i+1]-edges[i])/2)
    binwidth = np.append(binwidth, edges[i+1] - edges[i])
print("xdata", xData, "length", len(xData))
print("bin widths", binwidth, "length", len(binwidth))
#
# 
# 
## setup plot space ### this is just copied from madanalysis
fig_width=24
fig_aspect=1
fig = plt.figure(figsize=plt.figaspect(1))
width, height = fig.get_size_inches()
# 
gs = gridspec.GridSpec(2, 1, height_ratios=[2,1],  hspace=0, wspace=0, left=0.15, right=0.965, bottom=0.09, top=0.92)
title = 'Diboson  WZ ATLAS $c_W$ = 10,  $\sqrt{s}= 13 \, TeV$  '
plt.suptitle( title , fontsize=18)
padTop = plt.subplot(gs[0])
padBot = plt.subplot(gs[1])
inches_per_cm = 1./2.54
fig_height = fig_width*fig_aspect
fig.set_size_inches((fig_width*inches_per_cm, fig_height*inches_per_cm), forward=True)

h = padTop.hist(x=xData, bins=xBinning, weights= vals_1 ,\
        # density = True,\
        label= 'BSM1',\
        histtype="step", rwidth=1.0,\
        stacked=False,  hatch='**',
        color=rescolors[0], edgecolor=rescolors[0], linewidth=1.6, linestyle="solid",\
        bottom= None, cumulative=False, align="mid", orientation="vertical",fill=True,alpha=0.3)
#print("padtop = ", h)
#
padTop.hist(x=xData, bins=xBinning, weights= vals_2 ,\
        # density = True,\
        label= 'BSM2',\
        histtype="step", rwidth=1.0,\
        stacked=False,  hatch='//',
        color=rescolors[1], edgecolor=rescolors[1], linewidth=1.6, linestyle="solid",\
        bottom= None, cumulative=False, align="mid", orientation="vertical", fill=True,alpha=0.3 )

padTop.hist(x=xData, bins=xBinning, weights= vals_3 ,\
        # density = True,\
        label= 'SM',\
        histtype="step", rwidth=1.0, stacked=False, # hatch='//',
        color=rescolors[2], edgecolor=rescolors[2], linewidth=1.6, linestyle="solid",fill= False ,\
        bottom= None  , cumulative=False, align="mid", orientation="vertical")

padTop.hist(x=xData, bins=xBinning, weights= vals_4 ,\
        # density = True,\
        label= 'VBS',\
        histtype="step", rwidth=1.0, stacked=False, #hatch='**',
        color=rescolors[3], edgecolor=rescolors[3], linewidth=1.6, linestyle="solid",\
        #fill=True,alpha=0.3  ,\
        fill = False, 
        bottom= None  , cumulative=False, align="mid", orientation="vertical")


ylabel = r'$ \Delta \sigma~[\mathrm{fb}] $ '
xlabel = r' $  M^T_{WZ}~[\mathrm{GeV}]  $ '

padTop.set_ylabel(ylabel,fontsize=18,color="black", y = 0.75)
ymax = 6
ymin=  0.02 # log scale
print("ymin, ymax: ", ymin, ymax)
padTop.set_ylim(ymin, ymax)
padTop.set_xlim(xmin, xmax)
padTop.set_yscale('linear')
padTop.set_xscale('linear')
padTop.tick_params(axis='y', pad=1)
padTop.grid(False)
padTop.legend(loc=1,fontsize=15)
padTop.tick_params(which='both',direction='in',labelsize=15,right=True, left=True)
padTop.tick_params(which='major',length=3)
padTop.tick_params(which='minor',length=2)

# ratio plots #
padBot.sharex = padTop
padBot.sharey = padTop

hh = padBot.hist(x=xData, bins=xBinning, weights= vals_1/(np.sum(vals_1)) ,
            color=rescolors[0], edgecolor=rescolors[0],  histtype="step",\
                linestyle="solid", lw=2, hatch = '**' ,\
                bottom=None, cumulative=False, align="mid", orientation="vertical", fill=True,alpha=0.3 )
print("Padbot = ", hh)

padBot.hist(x=xData, bins=xBinning, weights= vals_2/(np.sum(vals_2)) ,
        color=rescolors[1], edgecolor=rescolors[1],  histtype="step",\
        linestyle="solid", lw=2,  hatch = '//' ,\
        bottom=None, cumulative=False, align="mid", orientation="vertical", fill=True,alpha=0.3 )

padBot.hist(x=xData, bins=xBinning, weights=  vals_3/(np.sum(vals_3) ) ,
        color=rescolors[2], edgecolor=rescolors[2],  histtype="step",\
        linestyle="solid", lw=2 ,\
        bottom=None, cumulative=False, align="mid", orientation="vertical", fill=False )

padBot.hist(x=xData, bins=xBinning, weights=  vals_4/(np.sum(vals_4) ) ,
        color=rescolors[3], edgecolor=rescolors[3],  histtype="step",\
        linestyle="solid", lw=2,  \
        bottom=None, cumulative=False, align="mid", orientation="vertical")

padBot.hlines(0, 0, edges[-1], colors="black", linestyles='solid',lw=2)

padBot.set_ylim(0.0, 1.1)
padBot.set_xlim(xmin,xmax)


padBot.set_xscale('linear') #test
padBot.set_xlabel( xlabel ,fontsize=20,color="black", x=1.3)
padBot.set_ylabel( 'normalised to 1' ,fontsize=18,color="black",ha='center', x=1)
padBot.tick_params(axis='both', pad=5)
padBot.grid(False)

plt.setp(padTop.get_xticklabels(), visible=False)
plt.setp(padBot.get_yticklabels(), family='serif')
plt.setp(padBot.get_xticklabels(), family='serif')

padBot.tick_params(which='both',direction='in',labelsize=15,right=True)
padBot.tick_params(which='major',length=3)
padBot.tick_params(which='minor',length=2)


plt.tight_layout(pad=1, w_pad=1, h_pad=1.0)

if not os.path.exists('sensitivityplots_desp_ratios/PNG'):
    os.makedirs('sensitivityplots_desp_ratios/PNG')
if not os.path.exists('sensitivityplots_desp_ratios/PDF'):
    os.makedirs('sensitivityplots_desp_ratios/PDF/')
plt.savefig('sensitivityplots_desp_ratios/PDF/WZ_atlas_MtWZ_ratios.pdf')
plt.savefig('sensitivityplots_desp_ratios/PNG/WZ_atlas_MtWZ_ratios.png')
# 
# 
# 
# # 
quit()
# 
