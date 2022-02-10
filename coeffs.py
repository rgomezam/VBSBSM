import yoda
import matplotlib as mpl 
import numpy as np


# Best fit values 

# cw = .....
# 
# 
# 
# 
#  Despoina

file_sm = yoda.read('SM_WZ_dib.yoda')
hist_sm = file_sm[ '/WZ_ATLAS/d12-x01-y01']
edges = hist_sm.xEdges()
print("xEdges: ",  edges )
vals_lo = hist_sm.areas()
# print("vals lo", vals_lo, "XS (in fb)", np.sum(vals_lo))
print("Check lengths:", len(edges), len(vals_lo))
#
file_data = yoda.read('WZ_ATLAS.yoda')
hist_data = file_data [ '/REF/WZ_ATLAS/d12-x01-y01' ]
vals_data = hist_data.yVals()
print("vals DATA", vals_data, "XS (in fb)", np.sum(vals_data))
# 
print("vals lo", vals_lo, "XS (in fb)", np.sum(vals_lo))
# 
file_eft = yoda.read('cW_WZ_dib_linear.yoda')
hist_eft = file_eft ['/WZ_ATLAS/d12-x01-y01']
# normalise by 10 
vals_eft = vals_lo + ((hist_eft.areas() - vals_lo)/10)
print("vals EFT (normalised to c+1)", vals_eft, "XS lin EFT (in fb)", np.sum(vals_eft))
# 
file_quad = yoda.read('cW_WZ_dib_quad.yoda')
hist_quad = file_quad ['/WZ_ATLAS/d12-x01-y01']
vals_quad = vals_eft + (hist_quad.areas()/100)
print("vals quad (normalised to c=1)", vals_quad, "XS quad (in fb)", np.sum(vals_quad))
# 