import numpy as np
import scipy as sc
import os
import matplotlib.pyplot as plt








## Prompt user for input
print("I'll need to know some things first")
R0 = input ("Enter major radius (float, i.e. 1.0)")
a0 = input ("Enter vertical offset (float, i.e. 1.0)")
R0 = float(R0)
a0 = float(a0)











### First part of code import data into handy-dandy format
print('Importing data...')

# Get current directory from OS
top = os.getcwd()

# Make array containing all the file names
files = []
# r=root, d=directories, f = files
for r, d, f in os.walk(top):
    for file in f:
        if ('out') in file:
            if ('.out') not in file:
                if ('-') not in file:
                    files.append(file)



# Create list containing numbers
numbers= []
for idx, f in enumerate(files):
    f = f.replace('out', '')
    numbers.append(f)
numbers = np.asarray(numbers)
numbers = numbers.astype(int)



# Check how many files there are
n_files = len(numbers)


# Make list containing all R Z coordinates. arrays[x] corresponds to outx
arrays = [0]*n_files
for idx, i in enumerate(numbers):
    vals      = np.loadtxt(files[idx],usecols = (1,2))
    arrays[i] = vals


    
    
    
## Convert data to r,theta coordinates
print("Performing coordinate transform...")
for i in range(0,n_files):
    RZ_coor    = arrays[i]
    R_coor     = RZ_coor[:,0]
    Z_coor     = RZ_coor[:,1]
    r_coor     = np.sqrt((R_coor - R0)**2 + (Z_coor - a0)**2)
    theta_coor = np.arctan2((Z_coor-a0),(R_coor-R0))
    new_coors  = np.transpose(np.array([r_coor,theta_coor]))
    new_coors  = new_coors[np.argsort(new_coors[:, 1])]
    arrays[i]  = new_coors



maxarr = [0]*n_files
# Calculate island width
print("Calculating island width...")

# We check for the maximum excursion in r for an theta-ordered array
for i in range(0,n_files):
    coors         = arrays[i]
    R_coors       = coors[:,0]
    R_coors_shift = np.roll(R_coors,1)
    diff          = abs(R_coors - R_coors_shift)
    maxdr         = np.amax(diff[1:])
    maxarr[i]     = [int(i),maxdr]

## We order the array based on the largest excursion - the number corresponding to the largest is expected to be an island seperatrix
maxarr     = np.array(maxarr)
maxarr     = maxarr[np.argsort(maxarr[:, 1])]
island_sep = int(maxarr[-1,0])
island_w   = maxarr[-1,1]


print(island_w)
print(island_sep)


plotexam = arrays[island_sep]
plt.scatter(plotexam[:,1],plotexam[:,0])
plt.show()
