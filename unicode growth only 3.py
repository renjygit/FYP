# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 13:30:13 2023

@author: BWCoo
"""

import numpy as np
import matplotlib.pyplot as plt
import statistics as stats
import sys

# Define the Van Leer flux limiter function
def vanLeerFunc(X):
    #print((X + np.abs(X)) / (1 + X))
    return (X + np.abs(X)) / (1 + X)


# Define required constants
A = 4e3
PMEa = 4.5e4 # J mol^-1
R = 8.314
VmMolar = 1.65e-5 # m^3 mol^-1, molar volume of monomer
rmMolar = (3/4*np.pi) * pow(VmMolar, 1/3) # effective radius of monomer (from the molar volume)
rm = 1.87e-10
Vm = (4/3)*np.pi*pow(rm,3)
D = 1e-12
NA = 6.022e23
Vtot = 0.0252e-3
MInf = 1e-2 #mol m^-3 #
u = 0.42 # coagulation parameter
gamma = 0.4 # J m^-2
kB = 1.38e-23
k = 0.1 # 1e12 # J mol^-1 m^-1 distribution exponent constant --------------------------------------------
kr = 1.6e-9
Damkohler = 1e-3
density = 5.8e3 # kg m^-3
MW = 0.096 # kg mol^-1
Q = (4*np.pi*density) / 3 * MW * Vtot * MInf

# Define time bins
tmax = 5e-6 # 3600
tdiff = 1e-6 #1e-3
timeArray = np.linspace(tdiff, tmax, int(tmax/tdiff))
timeCount = 0


# Define radius bins
rdiff =  3e-11 #0.3 Angstroms
rmax = 8.1e-9
rBins = np.linspace(0.5e-9,rmax,int((rmax+rdiff)/rdiff))

# Define temp variable & heating rate
Temp = 503 #Temp = 453
Tf = 567
HR = 0.025 # K s^-1 (1.5 K min^-1)

# Define precursor population variable & array
Ppop = 1e24 # change
PPopArray = [Ppop]

# Define supersaturation population variable & array
SS = 5e3
SSArray = [SS]

# Define nuclei arrays
NArraysArray = [[]]
NConc = 60e-6 # moles of nanocrystals initially
NNum = NConc * NA # number of nanocrystals initially
print("NNum = " + str(NNum))
# Create an initial distribution of nanoparticles (to analyse growth stage)
sum = 0
for r in rBins:
    g = (1 / (1e-10 * np.sqrt((2*np.pi)))) * np.exp((-np.power((r - 1e-9), 2)) / (2*(1e-10)**2)) # distribution with mean radius 1*10^-9 standard deviation of 10%
    sum+=g*rdiff # check probabilities sum to NNum (approx)
    #print(g)
    NArraysArray[0].append(g) # g going to 0 for all r
print("sum " + str(sum))
# Define avg array
#NArraysArray = [[0] * len(rBins)]
NArrayAvgR = [stats.fmean(NArraysArray[0])]
    
#print("initial N array " + str(NArraysArray[0]))


# Main loop steps forward in time
for time in timeArray: #start at timestep tdiff not 0
    timeCount +=1 # Keep track of number of time iterations
    #print(timeCount)
    
    plt.plot(rBins, NArraysArray[0])
    plt.title("Nanoparticle radius distribution t = " + str(time-tdiff))
    plt.show()
    
    # Increase the temperature by the heating rate up to the maximum
    if (Temp < Tf):
        Temp += tdiff*HR
    #print("temp " +str(Temp))
    
    
    """
    # Calculate the change to the precursor population and add the result to an array
    Ppop -= A*np.exp(-PMEa/(R*Temp)) * Ppop * tdiff
    PPopArray.append(Ppop)
    """

    """
    # De-dimensionalised values
    phi = (R * Temp) / (2 * gamma *Vm)
    psi = phi**2 * D * Vm * MInf
    Damkohler = (D * phi) / kr
    tau = time * psi
    """  
    """
    # Calculate critical radius, 'p' exponent and thus nucleation rate.
    rCrit = (2*gamma*VmMolar) / (R*Temp*np.log(SS))
    print("rCrit " + str(rCrit))
    p = np.power((rCrit/rm), 3) # Need actual value of rm not molar value
    # Use logarithms to reduce computation (exponents very large)
    lnRnuc = np.log(8*np.pi*rm*D*NA) + (p*u+1) * np.log(SS) * 2*np.log(MInf)  - ((4*np.pi*np.power(rCrit, 2)*gamma)/(3*kB*Temp)) #Need actual value of rm not molar value
    
    #print("Rnuc exponent " +str(np.exp((-4*np.pi*np.power(rCrit, 2)*gamma)/(3*kB*Temp))))
    #print("Rnuc " + str(Rnuc))
    """
    
    # Calculate standard deviation for this temperature.
    sigma = np.sqrt(kB*Temp/2)
    #print("sigma " + str(sigma))
    
    
    # Add a new empty array to N array to hold newly calculated values.
    NArraysArray.append([])
    # Create empty array for SS integral elements (resets each iteration)
    SSSumsArray = []
    NHalfPosList = []
    NHalfNegList = []
    
    rCount = -1 # Keep track of number of radius iterations
    for r in rBins:
        rCount += 1 # Keep track of number of radius iterations
        
        #beta = r * phi # de-dimensionalised radius
        
        Nprime = NArraysArray[0][rCount] + 0 # 0 for diffusion case
        #print("Nprime " + str(Nprime))
        
        # Set N to 0 for end cases since deltapos/neg need (i+1) or (i-1) values of radius, which don't exist at the boundaries.
        if r == rBins[0] or r == rBins[len(rBins) - 1]:
            NHalfPos = NArraysArray[0][rCount]
            NHalfNeg = NArraysArray[0][rCount]
            NHalfPosList.append(NHalfPos)
            NHalfNegList.append(NHalfNeg)
        else:
            """
            # Calculate value of distribution for this temperature and radius. Again use logarithms.
            lng = -np.log(sigma*np.sqrt((2*np.pi))) - ((-k*np.power((r-rCrit), 2)) / (2*np.power(sigma, 2))) 
            #print("g " + str(g))
            
            # Calculate intermediate value of N using N value at previous timestep
            #print("lnRnuc + lng " + str(lnRnuc + lng))

            """   
            

            

            
            # Calculate deltapos and deltaneg
            deltapos = NArraysArray[0][rCount+1] - NArraysArray[0][rCount]
            deltaneg = NArraysArray[0][rCount] - NArraysArray[0][rCount-1]
            #print("deltapos " + str(deltapos))
            
            # If either delta = 0 one or both of their ratios will diverge, so need a condition to avoid this.
            if deltapos == 0 or deltaneg == 0:
                NHalfPos = NArraysArray[0][rCount]
                NHalfNeg = NArraysArray[0][rCount]
                NHalfPosList.append(NHalfPos)
                NHalfNegList.append(NHalfNeg)
            else:
                # Calculate N at radii sizes (i+1/2) and (i-1/2) using van Leer
                NHalfPos = NArraysArray[0][rCount] + 0.5*vanLeerFunc(deltapos/deltaneg)*deltaneg
                NHalfPosList.append(NHalfPos)
                NHalfNeg = NArraysArray[0][rCount] - 0.5*vanLeerFunc(deltaneg/deltapos)*deltapos
                NHalfNegList.append(NHalfNeg)
                #print("NHalfPos  " + str(NHalfPos))
                
        
        # Calculate the growth rates at radii (i+1/2) and (i-1/2)
        growthRatepos = (D * VmMolar * SS * (1-np.exp((2*gamma*VmMolar)/((r+0.5*rdiff)*R*Temp)))) / ((r+0.5*rdiff) + Damkohler) # Find growth rate for i+1/2, i.e. current radius plus half a step
        growthRateneg = (D * VmMolar * SS * (1-np.exp((2*gamma*VmMolar)/((r-0.5*rdiff)*R*Temp)))) / ((r-0.5*rdiff) + Damkohler) # Find growth rate for i-1/2, i.e. current radius minus half a step
        #print("r " + str(r))
        #print("growthRatepos exponent " + str((2*gamma*VmMolar)/((r+0.5*rdiff)*R*Temp)))
        #print("growthRatepos e " + str(np.exp((2*gamma*VmMolar)/((r+0.5*rdiff)*R*Temp))))
        #print("growthRatepos 1-e " + str(1-np.exp((2*gamma*VmMolar)/((r+0.5*rdiff)*R*Temp))))
        #print("growthRatepos  " + str(growthRatepos))
            
        # Check Courant condition - just prints warning at present
        if (growthRatepos * tdiff/rdiff) > 1 or (growthRateneg * tdiff/rdiff) >1:
            print("Courant condition not satisfied.")
            print("Courant value for growthratepos " + str(growthRatepos * tdiff/rdiff))
            print("Courant value for growthrateneg " + str(growthRateneg * tdiff/rdiff))
            sys.exit()        
        
        # Calculate new N & add to the array
        NNew = Nprime - (tdiff/rdiff) * ((growthRatepos*NHalfPos) - (growthRateneg*NHalfNeg)) # Calculate new N, i.e. N for time (k+1)
        NArraysArray[1].append(NNew)
        #print("Narraysarray second element " + str(NArraysArray[1][1]))
        #print("   ")

        #Calculate element of SS integral for current r & add to array
        SSIntegralElement = pow(r,3) * (1/Vtot) * (NArraysArray[1][rCount] - NArraysArray[0][rCount])
        SSSumsArray.append(SSIntegralElement)

    plt.plot(rBins, NHalfPosList)
    plt.title("NHalfPos at t = " + str(time))
    plt.show()
    plt.plot(rBins, NHalfNegList)
    plt.title("NHalfNeg at t = " + str(time))
    plt.show()
    
    # Add average radius of N to array
    NArrayAvgR.append([stats.fmean(NArraysArray[1])])
    #print("Avg N " + str(NArrayAvgR[timeCount]))
    # Delete N array for previous timestep (to conserve memory since it is no longer needed)
    del NArraysArray[0]

    """
    # Calculate molar precursor concentration
    PConc = Ppop / (Vtot * NA)
    """
    
    # Calculate SS value for this timestep & add to array
    SS += - (Q * rdiff*np.sum(SSSumsArray)) # Sum all integral element values to approximate the integral, no precursor part for growth (just have an initial pop of monomer)
    SSArray.append(SS)
    print("SS " + str(SS))
    #print("   ")
    #print(NArraysArray[0])
    


timeArrayFull = np.linspace(0, tmax, int((tmax/tdiff) + 1)) # Create new time array to include t = 0

"""
# Plot precusor population against time
plt.plot(timeArrayFull, PPopArray, label="Heating rate = x K/min")
plt.title("Precursor population against time")
plt.show()
"""

plt.plot(rBins, NArraysArray[0])
plt.show()

# Plot supersaturation against time
plt.semilogx(timeArrayFull, SSArray, label="Heating rate = x K/min")
plt.title("Supersaturation against time")
plt.show()

