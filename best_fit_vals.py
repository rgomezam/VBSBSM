import yoda
import matplotlib as mpl 
import numpy as np


# Best fit values from SMEFiT

cpB_ind = [-0.005,0.002, 'cpB']
cpB_mar = [-0.739,0.289, 'cpB']

cpW_ind = [-0.018,0.007, 'cpW']
cpW_mar =[-0.592,0.677, 'cpW']

cpWB_ind = [-2.905,0.490, 'cpWB']
cpWB_mar = [-0.462,0.694, 'cpWB']

cpd_ind = [-0.428,1.214, 'cpd']
cpd_mar = [-2.002,3.693, 'cpd']

cpD_ind = [-4.066,2.657, 'cpD']
cpD_mar = [-1.498,0.974, 'cpD']

cWWW_ind = [-1.057,1.318, 'c3W']
cWWW_mar = [-1.049,1.459, 'c3W']
 

# For the best fit values we take the marginalized
def best_fit(array):
    bf = np.round((array[0] + array[1])/2,2)
    print( array[2], array[0:2] , bf)


ops = [cpB_mar, cpW_mar, cpWB_mar, cpd_mar, cpD_mar, cWWW_mar]


print("Marginalised, best fit")
for i in range(len(ops)):
    best_fit(ops[i])


print("Central values")
for i in range(len(ops)):
    tmp = np.round((ops[i][0] + ops[i][1])/2,2)
    central_value = ops[i][2] + " = " + str(tmp) + " "
    print(central_value)