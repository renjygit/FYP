# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:16:23 2023

@author: BWCoo
"""
import numpy as np
import matplotlib.pyplot as plt

# Define constants
R = 8.314
VmMolar = 3.29e-5 # m^3 mol^-1, molar volume of monomer
rmMolar = (3/4*np.pi) * pow(VmMolar, 1/3) # effective radius of monomer (from the molar volume)
rm = 1.87e-10
Vm = (4/3)*np.pi*pow(rm,3)
D = 1e-11 #--
NA = 6.022e23
Vtot = 0.0252e-3
MInf = 0.1 #mol m^-3 #
u = 0.42 # coagulation parameter
gamma = 0.4 # J m^-2
kB = 1.38e-23
k = 1e12 # J mol^-1 m^-1 distribution exponent constant --------------------------
kr = 1.6e-9 # -------------------------------------------

# Set temperature value
Temp = 500


# Produce 3 plots of instantaneous growth rate against radius, for varying Damkohler numbers, supersaturations and gamma values.

# Produce first plot, for varying supersaturations.

# Setting the radius bins
rdiff =  3e-14 #0.3 Angstroms
rmax = 8.1e-9
rBins = np.linspace(1e-9, rmax,int(rmax/rdiff))

# Create empty array to hold growth rate values
growthRateListList = []

# Create supersaturation bins to iterate through
SSBinsList = np.linspace(2, 200, 4)

# Set gamma & Damkohler number values
gamma = 0.2
damkohler = 1

SSCount = -1
#
for SS in SSBinsList:
    SSCount += 1
    growthRateListList.append([])
    for r in rBins:
        growthRate = ((D * VmMolar * SS * (1-np.exp((2*gamma*VmMolar)/((r)*R*Temp)))) / ((r) + damkohler))
        growthRateListList[SSCount].append(growthRate)
    
# Plot growth rate against radius
for n in np.linspace(0,SSCount,(SSCount+1)):
    #fig = plt.figure()
    #plt.subplot(2, 1, 1)
    plt.plot(rBins, growthRateListList[int(n)], label = "Supersaturation = " + str(SSBinsList[int(n)]))
    plt.xlabel("Radius")
    plt.ylabel("Instantaneous Growth Rate")
    plt.axhline(y=0, color='black')
    plt.legend()
    plt.yscale('symlog')
    plt.show()





# Produce second plot, for varying gamma values.

# Create empty array to hold growth rate values
growthRateListList = []


# Create gamma bins to iterate through
gammaBinsList = np.linspace(0.075, 1, 6)

# Set supersaturation and Damkohler number values
SS = 100
damkohler = 1

gammaCount = -1
#
for gamma in gammaBinsList:
    gammaCount += 1
    growthRateListList.append([])
    for r in rBins:
        growthRate = (D * VmMolar * SS * (1-np.exp((2*gamma*VmMolar)/((r)*R*Temp)))) / ((r) + damkohler)
        growthRateListList[gammaCount].append(growthRate)

# Plot growth rate against radius
for n in np.linspace(0,gammaCount,(gammaCount+1)):
    plt.plot(rBins, growthRateListList[int(n)], label = "gamma = " + str(gammaBinsList[int(n)]))
    plt.xlabel("Radius")
    plt.ylabel("Instantaneous Growth Rate")
    plt.yscale('symlog')
    plt.axhline(y=0, color='black')
    plt.legend()
    plt.show()





# Produce third plot, for varying Damkohler numbers.

# Create empty array to hold growth rate values
growthRateListList = []

# Create Damkohler bins to iterate through
damkohlerBinsList = [0.001, 0.01,0.1,1,10,100,1000]

# Set supersaturation and gamma values
SS = 100
gamma = 0.2

damkohlerCount = -1
#
for damkohler in damkohlerBinsList:
    damkohlerCount += 1
    growthRateListList.append([])
    for r in rBins:
        #print(np.exp((2*gamma*VmMolar)/((r)*R*Temp)))
        growthRate = (D * VmMolar * SS * (1-np.exp((2*gamma*VmMolar)/(r*R*Temp)))) / ((r) + damkohler)
        print("r " + str(r))
        print("growthRatepos exponent " + str((2*gamma*VmMolar)/(r*R*Temp)))
        print("growthRatepos e " + str(np.exp((2*gamma*VmMolar)/(r*R*Temp))))
        print("growthRatepos 1-e " + str((1-np.exp((2*gamma*VmMolar)/(r*R*Temp)))))
        growthRateListList[damkohlerCount].append(growthRate)

# Plot the growth rate against radius
for n in np.linspace(0,damkohlerCount,(damkohlerCount+1)):
    plt.plot(rBins, growthRateListList[int(n)], label = "damkohler = " + str(damkohlerBinsList[int(n)]))
    plt.yscale('symlog')
    plt.xlabel("Radius")
    plt.ylabel("Instantaneous Growth Rate")
    plt.axhline(y=0, color='black')
    plt.legend()
    plt.show()