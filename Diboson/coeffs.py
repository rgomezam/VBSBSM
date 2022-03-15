import yoda
import matplotlib as mpl 
import numpy as np
import math
#from termcolor import colored

# Best fit values:
# cpB = -0.22  cpW = 0.04  cpWB = 0.12 cpd = 0.85  cpD = -0.26  c3W = 0.21 
# Note the sign of cpd needs to be inverted to match our conventions 
#
#  Despoina

dirs = ['WZ/despoina/Inclusive/SM/', 'WZ/despoina/Inclusive/Inclusive_yodas_BSM1/', 'WZ/despoina/Inclusive/Inclusive_yodas_BSM2/']
analysis = ['/TESTDET/', '/TESTDET_BSM1/', '/TESTDET_BSM2/']


operators = [  'sm', 'cW' , 'cHW',  'cHB', 'cHDD', 'cHWB'] #, 'cHWtil', 'cHBtil', 'cHWBtil','cWtil']
best_fit_vals = [1 ,  0.21,  0.04 ,  -0.22 , -0.26 , 0.12 ] #,      1,        1,        1,        1  ]
#
files = ['d12-x01-y01', 'd08-x01-y01','d10-x01-y01','d14-x01-y01', 'd16-x01-y01' ]
distribs = ['MT_WZ', 'PT_Z', 'PT_W','Delta Phi WZ' , 'PT_nu'] 
# the distribs array is hard coded, comparing the hepdata entried with the
# labeling in the yoda files or directly comparing with the rivet analysis.cc file (See  )


# missing distribs 'd05-x01-y01' -> 'fid_XS_ratio'


#define sensitivity, to circumbvent many 0/0 results
def sensit(x,y): # x and y are arrays
    sensit_array = np.empty(len(x))
    for i in range(len(x)):        
        # if (x[i] == 0 and y[i] != 0) or  (x[i] != 0 and y[i] == 0):
        #     print("weird: EFT=", x[i], ' SM=', y[i])    
        if  x[i] == 0 or y[i]==0:
            sensit_array[i]= 00.00 # testmath.pi
        else:    
            sensit_array[i] = (x[i]/y[i]).round(decimals=3)   
    return sensit_array


# First print total XSEC, then loop over the other distributions
print('###########    Total XSEC:    ##########################' )
print('\n ### we first print the SM as a check ### ')
# now loop over operators:
for op,bf in zip(operators,best_fit_vals):
    print('Operator: ' + op + ", best fit val: "+ str(bf))
    for dir,an in zip(dirs,analysis):
        hist_sm = yoda.read(dir+'sm.yoda')[an + files[0]] #any histo is fine, we take 0 for example
        vals_sm = hist_sm.areas()
        filename = yoda.read(dir +  op  + '.yoda')
        hist = filename[ an + files[0]] #any histo is fine, we take 0 for example
        vals_lo = (hist.areas())
        print( an + " Total XS (in fb)", np.sum(vals_lo).round(decimals=3) ,  " SM XS (in fb)", np.sum(vals_sm).round(decimals=3) ,
            ' ratio to SM (in %, for c=10): ',  (1- np.sum(vals_lo)/np.sum(vals_sm) ).round(decimals=2) *100  , 
              "best fit ratio: (in %, for c="+ str(bf) +")" ,  ((1- np.sum(vals_lo)/np.sum(vals_sm) ) *(bf/10) ).round(decimals=3) *100       )   
    print('\n \n')



#loop over distributions
for i in range(len(files)):
    print('###########     Distribution:  ' + distribs[i] + '  ##########################' )
    print('\n ### we first print the SM as a check ### \n ')
    
    for op,bf in zip(operators,best_fit_vals) :
        print('Operator: ' + op , "best fit val: ", bf)
        for dir,an in zip(dirs,analysis):
            hist_sm = yoda.read(dir+'sm.yoda')[an + files[i] ] 
            vals_sm = hist_sm.areas()
            filename = yoda.read(dir +  op  + '.yoda')
            hist = filename[ an + files[i]] 
            vals_lo = hist.areas()
            print(an + "vals linear EFT", vals_lo.round(decimals=3))
            print(an + "Sensitivities (in %, for c=10): ", np.absolute(1-sensit(vals_lo,vals_sm))*100 )
            print(an + "Best fit Sensitivities (in %, for c="+ str(bf) +")",   np.absolute( (bf/10)*(1-sensit(vals_lo,vals_sm)))*100 )
            print("\n")
        print( '\n \n')



quit()
    # now loop over operators (also over SM for completeness)
#     for op in operators:
#         filename = yoda.read(dir1 + op +'.yoda')
#         hist = filename[ '/TESTDET/'+ files[i]]
#         edges = hist.xEdges()
#         #print("xEdges: ",  edges )
#         vals_lo = (hist.areas()).round(decimals=3)
#         #
#         print('Operator: ' + op )
#         print("vals linear EFT", vals_lo)
#         print("Sensitivities", (vals_lo/vals_sm).round(decimals=3) , '\n \n')
        

# quit()

# # Now we do the same for the BSM1/BSM2 selections
