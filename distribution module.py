# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 11:21:42 2023
@author: BWCoo
"""

import numpy as np
import matplotlib.pyplot as plt

# Define required constants
A = 4e3
PMEa = 4.5e4 # J mol^-1
R = 8.314
VmMolar = 1.65e-5 # m^3 mol^-1, molar volume of monomer
rmMolar = (3/4*np.pi) * pow(VmMolar, 1/3) # effective radius of monomer (from the molar volume)
rm = 1.87e-10
Vm = (4/3)*np.pi*pow(rm,3)
D = 1e-10
NA = 6.022e23
Vtot = 0.0252e-3
MInf = 0.1 #mol m^-3 #
u = 0.42 # coagulation parameter
gamma = 0.4 # J m^-2
kB = 1.38e-23
kr = 1#1.6e-9 # -------------------------------------------
density = 5.8e3 # kg m^-3
MW = 0.096 # kg mol^-1
Q = (4*np.pi*density) / 3 * MW * Vtot * MInf

# Define temp variable & heating rate
Temp = 453
Tf = 567
HR = 0.025 # K s^-1 (1.5 K min^-1)

# Define supersaturation population variable
SS = 5e3


sigma = np.sqrt(kB*Temp/2)


rCrit = (2*gamma*VmMolar) / (R*Temp*np.log(SS))
print("rCrit " + str(rCrit))
r = 0

"""
# Calculate g & its components & print (for inspection of values)
gnormaliseconst = (1/(sigma*np.sqrt((2*np.pi))))
gtest = (-k*np.power((r-rCrit), 2)) / (2*np.power(sigma, 2))
g = (1/(sigma*np.sqrt((2*np.pi)))) * np.exp((-k*np.power((r-rCrit), 2)) / (2*np.power(sigma, 2))) 
print("g " + str(g))
print("gnormaliseconstant " + str(gnormaliseconst))
print("g exponent value " + str(gtest))
"""

# setting the radius bins
rdiff =  3e-11 #0.3 Angstroms
rmax = 8.1e-9
rBins = np.linspace(0,rmax,int((rmax+rdiff)/rdiff)*1)

gArray = []
# Calculate the distribution across the r values in rBins
for r in rBins:
    giteration = (1 / (1e-10 * np.sqrt((2*np.pi)))) * np.exp((-np.power((r - 1e-9), 2)) / (2*(1e-10)**2)) # distribution with mean radius 1*10^-9 standard deviation of 10%

    #giteration = (1/(sigma*np.sqrt((2*np.pi)))) * np.exp((-k*np.power((r-rCrit), 2)) / (2*np.power(sigma, 2))) 
    gArray.append(giteration)

#g = np.random.normal(1e-9, 1e-10, size = 270)
#g = g.tolist()
#print(g)
# Plot the distribution
#plt.scatter(rBins,g)
plt.plot(rBins,gArray)
plt.show()
