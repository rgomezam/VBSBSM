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

file_sm = yoda.read('SM_WZ_dib.yoda')
hist_sm = file_sm[ '/WZ_ATLAS/d10-x01-y01']
edges = hist_sm.xEdges()
print("xEdges: ",  edges )
vals_lo = hist_sm.areas()
print("vals lo", vals_lo, "XS (in fb)", np.sum(vals_lo))
print("Check lengths:", len(edges), len(vals_lo))
#
file_data = yoda.read('WZ_ATLAS.yoda')
hist_data = file_data [ '/REF/WZ_ATLAS/d10-x01-y01' ]
vals_data = hist_data.yVals()
print("vals DATA", vals_data, "XS (in fb)", np.sum(vals_data))
# 
file_eft = yoda.read('cW_WZ_dib_linear.yoda')
hist_eft = file_eft ['/WZ_ATLAS/d10-x01-y01']
# normalise by 10 
vals_eft = vals_lo + ((hist_eft.areas() - vals_lo)/10)
print("vals EFT (normalised to c+1)", vals_eft, "XS lin EFT (in fb)", np.sum(vals_eft))
# 
file_quad = yoda.read('cW_WZ_dib_quad.yoda')
hist_quad = file_quad ['/WZ_ATLAS/d10-x01-y01']
vals_quad = vals_eft + (hist_quad.areas()/100)
print("vals quad (normalised to c=1)", vals_quad, "XS quad (in fb)", np.sum(vals_quad))
# 
# 
def plot_lin():
    #################################################################################################
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    rescolors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
    #################################################################################################
    # Creating data sequence: middle of each bin
    xBinning = edges[0:-1]
    print("xBinning:" , xBinning)
    xmin = edges[0]
    xmax = edges[-2]
    xData = []
    binwidth = []
    for i in range(len(edges)-2):
        xData = np.append(xData, edges[i]+(edges[i+1]-edges[i])/2)
        binwidth = np.append(binwidth, edges[i+1] - edges[i])
    print("xdata", xData, "length", len(xData))
    print("bin widths", binwidth, "length", len(binwidth))
    #
    M_1_weights = vals_lo[0:-1]
    M_2_weights = vals_eft[0:-1]
    M_3_weights = vals_quad[0:-1]
    #
    M_1_of_weights = np.array(vals_lo[-1:])
    M_2_of_weights = vals_eft[-1:]
    M_3_of_weights = vals_quad[-1:]
    xData_of = (edges[-1] - edges[-2])/2
    xBinning_of = np.array([edges[-2], edges[-1]])
    # 
    ## setup plot space ### this is just copied from madanalysis
    fig_width=24
    fig_aspect=1
    fig = plt.figure(figsize=plt.figaspect(1))
    width, height = fig.get_size_inches()
    # 
    gs = gridspec.GridSpec(2, 2, height_ratios=[2,1],\
    width_ratios=[2,1], hspace=0, wspace=0, left=0.15, right=0.965, bottom=0.09, top=0.92)
    title = 'Diboson  WZ ATLAS,  $\sqrt{s}= 13 \, TeV$  '
    plt.suptitle( title , fontsize=18)
    padTop = plt.subplot(gs[0,0])
    padBot = plt.subplot(gs[1,0])
    padTR = plt.subplot(gs[0,1])
    padBR = plt.subplot(gs[1,1])

    inches_per_cm = 1./2.54
    fig_height = fig_width*fig_aspect
    fig.set_size_inches((fig_width*inches_per_cm, fig_height*inches_per_cm), forward=True)

    h = padTop.hist(x=xData, bins=xBinning, weights= M_1_weights ,\
            # density = True,\
            label= r'$\rm{SM} \quad \,  p p \to \ell \bar{\ell}  \ell \bar{\nu} $',\
            histtype="step", rwidth=1.0,\
            stacked=False,  #hatch='/',
            color=rescolors[0], edgecolor=rescolors[0], linewidth=1.6, linestyle="solid",\
            bottom= None, cumulative=False, align="mid", orientation="vertical",fill=True,alpha=0.3)
    print("padtop = ", h)
    
    padTop.hist(x=xData, bins=xBinning, weights= M_2_weights ,\
            # density = True,\
            label= r'$\rm{linear \, EFT} \,  c_W = 1 $',\
            histtype="step", rwidth=1.0,\
            stacked=False,  hatch='//',
            color=rescolors[0], edgecolor=rescolors[0], linewidth=2.0, linestyle="solid",\
            bottom= None, cumulative=False, align="mid", orientation="vertical",fill=False)

    padTop.hist(x=xData, bins=xBinning, weights= M_3_weights ,\
            # density = True,\
            label= r'$\rm{quad \, EFT} \, c_{W} = 1 $',\
            histtype="step", rwidth=1.0, stacked=True, 
            #  hatch='//',
            color=rescolors[1], edgecolor=rescolors[1], linewidth=2.0, linestyle="solid",fill=False,\
            bottom= None  , cumulative=False, align="mid", orientation="vertical")

    ylabel = r'$ \Delta \sigma~[\mathrm{fb}] $ '
    xlabel = r' $  P^T_{W}~[\mathrm{GeV}]  $ '

    padTop.set_ylabel(ylabel,fontsize=18,color="black", y = 0.75)
    ymax = 20
    ymin=  0.51 # log scale
    print("ymin, ymax: ", ymin, ymax)
    padTop.set_ylim(ymin, ymax)
    padTop.set_xlim(xmin, xmax)
    padTop.set_yscale('log')
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
    padBR.sharey = padBot
    # padBot.sharey = padTop

    hh = padBot.hist(x=xData, bins=xBinning, weights= -1 + (M_2_weights )/( M_1_weights ) ,
                color=rescolors[0], edgecolor=rescolors[0],  histtype="step",\
                  linestyle="solid", lw=2,\
                  bottom=None, cumulative=False, align="mid", orientation="vertical")
    print("Padbot = ", hh)
    
    padBot.hist(x=xData, bins=xBinning, weights= -1 + (M_3_weights )/(M_1_weights ) ,
            color=rescolors[1], edgecolor=rescolors[1],  histtype="step",\
            linestyle="solid", lw=2,  \
            bottom=None, cumulative=False, align="mid", orientation="vertical")

    padBot.hlines(0, 0, edges[-1], colors="black", linestyles='solid',lw=2)

    padBot.set_ylim(-0.3, 0.6)
    padBR.set_ylim(-4.5, 9)
    padBot.set_xlim(xmin,xmax)
    padTR.set_xlim(edges[-2], edges[-1])
    padTR.set_ylim(0.6,9)
    padTR.set_yscale('log')
    padBR.set_xlim(edges[-2], edges[-1])
    padBR.sharex = padTR
    # padBR.set_yscale("log") 

    padBot.set_xscale('linear') #test
    padBot.set_xlabel( xlabel ,fontsize=20,color="black", x=1.3)
    padBot.set_ylabel( r'$\rm{Ratio  \, EFT/SM}$' ,fontsize=18,color="black",ha='center', x=1)
    padBot.tick_params(axis='both', pad=5)
    padBot.grid(False)

    plt.setp(padTop.get_xticklabels(), visible=False)
    plt.setp(padBot.get_yticklabels(), family='serif')
    plt.setp(padBot.get_xticklabels(), family='serif')

    padBot.tick_params(which='both',direction='in',labelsize=15,right=True)
    padBot.tick_params(which='major',length=3)
    padBot.tick_params(which='minor',length=2)

    # Overflow Bins
    # print("test", xData_of, xBinning_of, M_1_of_weights )
    # print("test 2", type(xData_of), type(xBinning_of), type(M_1_of_weights) )
    # print("test 3", type(xData), type(xBinning), type(M_1_weights) )
    hhh = padTR.hist(x=xData_of, bins=xBinning_of, weights= M_1_of_weights ,\
            # density = True,\
            label= r'$\rm{SM} \quad \,  p p \to \ell \bar{\ell}  \ell \bar{\nu} $',\
            histtype="step", rwidth=2.0,\
            stacked=False,  #hatch='/',
            color=rescolors[0], edgecolor=rescolors[0], linewidth=1.6, linestyle="solid",\
            bottom= None, cumulative=False, align="mid", orientation="vertical",fill=True,alpha=0.3)
    print("padtop overflow = ", hhh)

    padTR.hist(x=xData_of, bins=xBinning_of, weights= M_2_of_weights ,\
            # density = True,\
            color=rescolors[0], edgecolor=rescolors[0],  histtype="step",\
            linestyle="solid", lw= 2,\
            stacked=False,  hatch='//',
            bottom= None, cumulative=False, align="mid", orientation="vertical",fill=False)
    
    padTR.hist(x=xData_of, bins=xBinning_of, weights= M_3_of_weights ,\
            # density = True,\
            histtype="step", rwidth=1.0, \
            stacked=False,  #hatch='//',
            color=rescolors[1], edgecolor=rescolors[1], linewidth=2, linestyle="solid",\
            bottom= None, cumulative=False, align="mid", orientation="vertical",fill=False)
    

    padBR.hist(x=xData_of, bins=xBinning_of, weights= -1 + (M_2_of_weights )/( M_1_of_weights ) ,
                color=rescolors[0], edgecolor=rescolors[0],  histtype="step",\
                 linestyle="solid", lw=2,\
                bottom=None, cumulative=False, align="mid", orientation="vertical")

    
    padBR.hist(x=xData_of, bins=xBinning_of, weights= -1 + (M_3_of_weights )/(M_1_of_weights ) ,
        color=rescolors[1], edgecolor=rescolors[1],  histtype="step",\
        linestyle="solid", lw=2,  \
        bottom=None, cumulative=False, align="mid", orientation="vertical")

    padBR.hlines(0, 0, edges[-1], colors="black", linestyles='solid',lw=1.5)


    plt.tight_layout(pad=1, w_pad=1, h_pad=1.0)

    if not os.path.exists('sensitivityplots/PNG'):
        os.makedirs('sensitivityplots/PNG')
    if not os.path.exists('sensitivityplots/PDF'):
        os.makedirs('sensitivityplots/PDF/')
    plt.savefig('sensitivityplots/PDF/WZ_atlas_ptW.pdf')
    plt.savefig('sensitivityplots/PNG/WZ_atlas_ptW.png')
# # end def


plot_lin()
print("exiting")
quit()



def plotSens(interf, op, factor, label):
    #################################################################################################
    rescolors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
    #################################################################################################
    # Creating data sequence: middle of each bin
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
    y3_M_1_weights = vals_lo  # tmp2 = Monte Carlo
    y3_M_3_weights = (vals_lo + factor*interf) # tmp3 = EFT interf c=1
    print("PLOT: SM weights:", y3_M_1_weights, "length", len(y3_M_1_weights))
    print(op +" EFT weights:", y3_M_3_weights,  "length", len(y3_M_3_weights),  \
    "tot", np.sum(y3_M_3_weights))
    print(op +" EFT ratios:", y3_M_3_weights/y3_M_1_weights )
    #
    ### setup plot space ### this is just copied from madanalysis
    fig_width=18
    fig_aspect=1
    fig = plt.figure(figsize=plt.figaspect(1))
    width, height = fig.get_size_inches()
    gs = gridspec.GridSpec(2, 1, height_ratios=[2,1], hspace=0, left=0.15, right=0.965, bottom=0.09, top=0.92)
    padTop = plt.subplot(gs[0])
    title = 'Diboson  WZ ATLAS,  $\sqrt{s}= 13 \, TeV$  '
    plt.suptitle( title , fontsize=18)
    padTop = plt.subplot(gs[0])
    padBotA = plt.subplot(gs[1])

    inches_per_cm = 1./2.54
    fig_height = fig_width*fig_aspect
    fig.set_size_inches((fig_width*inches_per_cm, fig_height*inches_per_cm), forward=True)


    padTop.hist(x=xData, bins=xBinning, weights= y3_M_1_weights ,\
               # density = True,\
               label= "SM ", histtype="step", rwidth=1.0,\
               stacked=False,  #hatch='/',
               color=rescolors[0], edgecolor=rescolors[0], linewidth=1.6, linestyle="solid",\
                bottom= None, cumulative=False, align="mid", orientation="vertical",fill=True,alpha=0.3)

    padTop.hist(x=xData, bins=xBinning, weights= y3_M_1_weights ,\
               # density = True,\
                histtype="step", rwidth=1.0,\
               stacked=False,  #hatch='/',
               color=rescolors[0], edgecolor=rescolors[0], linewidth=2.0, linestyle="solid",\
                bottom= None, cumulative=False, align="mid", orientation="vertical",fill=False)

    padTop.hist(x=xData, bins=xBinning, weights= y3_M_3_weights,\
               # density = True,\
               label= 'SM+EFT '+op+'='+ str(factor) +' TeV$^{-2})$', histtype="step", rwidth=1.0,\
                stacked=True, hatch='//',
                color=rescolors[1], edgecolor=rescolors[1], linewidth=2.0, linestyle="solid",fill=False,\
                bottom= None  , cumulative=False, align="mid", orientation="vertical")



    ymax = 100
    ymin=  0.001 # log scale
    print("ymin, ymax: ", ymin, ymax)

    totSM= np.sum(y3_M_1_weights)
    totEFT = np.sum(y3_M_3_weights)

    # ylabel = r'$ d \sigma / d P^T{Z}~[\mathrm{fb/GeV}] $ '
    ylabel = r'$ \Delta \sigma~[\mathrm{fb}] $ '
    xlabel = r' $  P^T_{Z}~[\mathrm{GeV}]  $ '

      # print ("ylabel is:", ylabel)
    padTop.set_ylabel(ylabel,fontsize=18,color="black")
    padTop.set_ylim(ymin, ymax)
    padTop.set_yscale('log')
    padTop.set_xscale('linear') #test
    padTop.tick_params(axis='y', pad=3)
    padTop.grid(False)
    padTop.legend(loc=1,fontsize=17)
    padTop.tick_params(which='both',direction='in',labelsize=15,right=True)
    padTop.tick_params(which='major',length=3)
    padTop.tick_params(which='minor',length=2)

    # ratio plots #
    padBotA.sharex = padTop
    padBotA.sharey = padTop

    padBotA.hist(x=xData, bins=xBinning, weights= -1 + (y3_M_3_weights )/( y3_M_1_weights ) ,
               color=rescolors[1], edgecolor=rescolors[1],  histtype="step",\
                 linestyle="solid", lw=2,\
                 bottom=None, cumulative=False, align="mid", orientation="vertical")

    padBotA.hlines(0, 0, 2000, colors=rescolors[0], linestyles='solid',lw=2)


    padBotA.set_ylim(-0.5, 0.5)
    padTop.set_xlim(xmin,xmax)
    padBotA.set_xlim(xmin,xmax)
    padBotA.set_xscale('linear') #test
    #padBotA.set_xlabel( xlabel ,fontsize=18,color="black",ha='center', x=1)
    padBotA.set_xlabel( xlabel ,fontsize=20,color="black")
    padBotA.set_ylabel( "$Ratio EFT/SM$" ,fontsize=18,color="black",ha='center', x=1)
    padBotA.tick_params(axis='both', pad=5)
    padBotA.grid(False)

    plt.setp(padTop.get_xticklabels(), visible=False)
    plt.setp(padBotA.get_yticklabels(), family='serif')
    plt.setp(padBotA.get_xticklabels(), family='serif')

    padBotA.tick_params(which='both',direction='in',labelsize=15,right=True)
    padBotA.tick_params(which='major',length=3)
    padBotA.tick_params(which='minor',length=2)

    plt.tight_layout(pad=1, w_pad=1, h_pad=1.0)

    if not os.path.exists('sensitivityplots2/PNG'):
        os.makedirs('sensitivityplots2/PNG')
    if not os.path.exists('sensitivityplots2/PDF'):
        os.makedirs('sensitivityplots2/PDF/')
    plt.savefig('sensitivityplots2/PDF/signal_only_'+op+label+'.pdf')
    plt.savefig('sensitivityplots2/PNG/signal_only_'+op+label+'.png')
# end def

############################################################
#### PLOT EACH OP
# line = "run summary: \n"
# operators_short = ['cW','cHW']
# for op in operators_short:
# #       #pick up the relative eft interference histo
#         eft_name = 'yodafiles/' +  op + '.yoda'
#         print("Doing op:", op, "file", eft_name)
#         temp = yoda.read(eft_name)[histoname]
#         vals_eft = (temp.areas())*2
#         # inter(vals_eft)
#         #line += op+' ratio to SM: '+str(inter(vals_eft)/vals_lo)+'\n'
#         # if np.max(interf/vals_lo) > 2. or np.min(interf/vals_lo) < - 2.:
#         #     print(op, np.max(interf/vals_lo) , np.min(interf/vals_lo) )
#         #     plotSens(interf, op, 0.1 , '0p1')
#         # elif np.max(interf/vals_lo) > 1.3 or np.min(interf/vals_lo) < - 1.3:
#         plotSens(inter(vals_eft), op, 0.5 , '0p5')
#         # elif np.max(interf/vals_lo) < 0.1 and np.min(interf/vals_lo) > -0.1:
#         #     print(op, np.max(interf/vals_lo) , np.min(interf/vals_lo) )
#         #     plotSens(interf, op, 5 , '5')
#         # else:
#         #     print(op, np.max(interf/vals_lo) , np.min(interf/vals_lo) )
#         #     plotSens(interf, op, 1. , '1')

#### PLOT 3 OPERATORS
# operators_short = ['cW','cHW','cHB']
plot3('cW','cHW','cHB',0.5,'0p5')
plot3('cW','cHW','cHB',0.1,'0p1')
#
# for op in operators_short:
# #       #pick up the relative eft interference histo
#         eft_name = 'yodafiles/' +  op + '.yoda'
#         temp = yoda.read(eft_name)[histoname]
#         vals_eft = (temp.areas())/6.176*1E03
#         interf = np.append(interf , ( vals_eft - vals_lo ) / (v*v)  / cwil)
#
# print("interf 1,2,3", operators_short , interf[0], interf[1], interf[2])
# print("length interf", len(interf))
# # def plot3(op1, interf1 , op2 , interf2 , op3 , interf3 , factor1, factor2, factor3, label)
# plot3(operators_short[0], interf[0] , operators_short[1] , interf[1] , operators_short[2] , \
#  interf[2] , 0.5 , 0.1 , 0.1,  '0p5')
# # plot2(operators_short[1] , interf[1] , operators_short[2] , interf[2] , 0.5 , '0p5')
# # #
# # #
# writefile(line, 'output.log')
