import numpy as np
import matplotlib.pyplot as plt
import sys

# Define the van Leer flux limiter function
def vanLeerFunc(X):
    return (X + np.abs(X)) / (1 + X)

# Define required constants
A = 4e3
R = 8.314
NA = 6.022e23
VmMolar = 1.65e-5 # m^3 mol^-1, molar volume of monomer
rmMolar = (3/4*np.pi) * pow(VmMolar, 1/3) # effective radius of monomer (from the molar volume)
Vm = VmMolar / NA #(4/3)*np.pi*pow(rm,3)
rm =  1.87e-10 # pow((3*Vm) / (4*np.pi), (1/3))
D = 5e-11
Vtot = 1e-4 #0.0252-3 # m^3
MInf = 4e-3 #mol m^-3 #
u = 0.48 # coagulation parameter
gamma = 1.1 #0.4 # J m^-2
kB = 1.38e-23
#k = 1e12 # J mol^-1 m^-1 distributio n exponent constant
kr = 1.6e-9
density = 5.816e3 # kg m^-3
MW = 0.09586 # kg mol^-1
Q = (4*np.pi*density) / (3 * MW * Vtot * MInf)
PMEnergy = 1

# Define timesteps to iterate through & create an array
tmax = 1 # s
tdiff = 5e-6
tmin = tdiff
timeArray = np.linspace(tdiff, tmax, int((tmax-tmin+tdiff)/tdiff))
timeCount = 0 # Used to track number of time iterations


# Define radius bins
rdiff =  3e-11 # m
rmax = 8.1e-9
rmin = 5e-10
rBins = np.linspace(rmin,rmax,int((rmax-rmin+rdiff)/rdiff))

# Define temperature variable & heating rate
Temp = 573 # 453
#Tf = 567
#HR = 0.025 # K s^-1 (1.5 K min^-1)

# Define precursor population variable & array
PPop = 1e24
PConc = PPop / (NA * Vtot)
PArray = [PConc]

# Define supersaturation variable & array
SS = 8e3
SSArray = [SS]

# Monomers
MConc = 32 # mol m^-3

NArraysArray = [np.zeros(len(rBins))]

NNumArray = [np.sum(NArraysArray[0])]

"""
# Create an array of arrays to hold number of nuclei in each radius bin for different timesteps
NArraysArray = [[]]
# Give number of nanopartices initially
NNumMole = 60e-6 # mol
NNum = NNumMole * NA # number
print("Initial number of nanoparticles = " + str(NNum))

# Create the input distribution of nanoparticles
distributionArray = []
sum = 0
for r in rBins:
    g = (1 / (1e-10 * np.sqrt((2*np.pi)))) * np.exp(-0.5*(np.power((r - 1e-9)/(1e-10), 2))) # Normal distribution with mean radius 1e-9 m & standard deviation 1e-10 m
    gN = NNum * g * rdiff
    
    sum += gN # Calculated to check this sums to NNum (approximately)
    
    distributionArray.append(g)
    NArraysArray[0].append(gN)
print("Sum of nanoparticles in initial distribution = " + str(sum))

# Plot the initial probability density function
plt.plot(rBins, distributionArray)
plt.title("Nanoparticle radius distribution t = 0")
plt.ylim(0, 4.5e9)
plt.show()
"""

# Plot the initial distribution of nanoparticles
plt.plot(rBins, NArraysArray[0])
plt.title("Nanoparticles number against radius at t = 0")
plt.ylim(0, 4.5e18)
plt.show()


# Create an array to hold the mean radius values
NArrayAvgR = [0] #[np.average(rBins, weights = NArraysArray[0])]
print("Initial avg N = " + str(NArrayAvgR[0]))

# Main loop steps forward in time
for time in timeArray:
    timeCount += 1

    # Calculate the change to the precursor population and add the result to an array
    PConc -= (A*np.exp(-PMEnergy/Temp) * PConc * tdiff)
    PArray.append(PConc)

    #print("temp " +str(Temp))
    
    # Calculate the critical radius
    rCrit = (2*gamma*VmMolar) / (R*Temp*np.log(SS))
    
    # Plot the nanoparticle distribution every few hundred iterations
    if timeCount % 400 == 0:
        plt.plot(rBins, NArraysArray[0])
        plt.title("Nanoparticles number against radius at t = " + f'{(time-tdiff):.3g}') # Take away tdiff since we are actually plotting the graph from the previous here
        #plt.ylim(0, 4.5e18)
        #plt.plot(rCrit,1e9,'ro')
        plt.show()
    
    # Calculate the de-dimensionalised variables
    phi = (R * Temp) / (2 * gamma *VmMolar)
    psi = phi**2 * D * VmMolar * MInf
    zeta = 1e5 # (D * phi) / kr
    tau = time * psi
    
    # Calculate p value
    p = np.power((rCrit/rm), 3) # Need actual value of rm not molar value
    # Calculate nucleation rate
    RNuc = 8*np.pi*rm*D*NA * pow(SS, (p*u+1)) * pow(MInf, 2) * np.exp((-4*np.pi*pow(rCrit, 2)*gamma)/(3*kB*Temp)) #Need actual value of rm not molar value

    # Add a new empty array to the N array of arrays to hold the newly calculated N values
    NArraysArray.append([])
    # Create empty arrays for the SS integral elements and the N values at radii values of (i+1/2) and (i-1/2) (reset each time iteration)
    SSSumsArray = []
    NHalfPosList = []
    NHalfNegList = []
    
    rCount = -1 # Keep track of the number of radius iterations (resets each time iteration)
    for r in rBins:
        rCount += 1
        
        # Calculate standard deviation
        sigma = np.sqrt((kB*Temp)/2) # (kB*Temp/2)
        # Calculate value of distribution for current r value
        g = (1 / (sigma * np.sqrt((2*np.pi)))) * np.exp((-np.power((r - rCrit), 2)) / (2*np.power(sigma, 2)))
        #print("g = " + str(g))
        #print("radius = " + str(r))
        
        # Calculate the intermediate value of N, i.e. the nucleated N (for this radius value)
        Nprime = NArraysArray[0][rCount] + (tdiff * RNuc * g)
        
        # Deltapos & deltaneg use number of nanoparticles at radii (i+1) & (i-1), which don't exist at the boundaries, so set NHalfPos and NHalfNeg to previous N for these edge cases
        if r == rBins[0] or r == rBins[len(rBins) - 1]:
            NHalfPos = NArraysArray[0][rCount]
            NHalfNeg = NArraysArray[0][rCount]
            NHalfPosList.append(NHalfPos)
            NHalfNegList.append(NHalfNeg)
        else:
            # Calculate deltapos and deltaneg
            deltapos = NArraysArray[0][rCount+1] - NArraysArray[0][rCount]
            deltaneg = NArraysArray[0][rCount] - NArraysArray[0][rCount-1]
            
            # If one or both of the deltas = 0, one or both of their ratios will diverge, so include a condition to avoid this
            if deltapos == 0 or deltaneg == 0:
                NHalfPos = NArraysArray[0][rCount] # Set NHalfPos and NHalfNeg values equal to the previous N value
                NHalfNeg = NArraysArray[0][rCount]
                NHalfPosList.append(NHalfPos)
                NHalfNegList.append(NHalfNeg)
            else:
                # Calculate NHalfPos and NHalfNeg using the van Leer function
                NHalfPos = NArraysArray[0][rCount] + 0.5*vanLeerFunc(deltapos/deltaneg)*deltaneg
                NHalfPosList.append(NHalfPos)
                NHalfNeg = NArraysArray[0][rCount] - 0.5*vanLeerFunc(deltaneg/deltapos)*deltapos
                NHalfNegList.append(NHalfNeg)
                
        # Find beta values for radii values (i+1/2) and (i-1/2)
        betaPos = (r+0.5*rdiff) * phi # de-dimensionalised radius
        betaNeg = (r-0.5*rdiff) * phi # de-dimensionalised radius

        # Calculate the growth rates using these beta values & then redimensionalise them
        growthRatePosDD = (SS - np.exp(1/betaPos))/(betaPos + zeta) # De-dimensionalised
        growthRateNegDD = (SS - np.exp(1/betaNeg))/(betaNeg + zeta) # De-dimensionalised
        growthRatePos = growthRatePosDD * (psi/phi) # Re-dimensionalised
        growthRateNeg = growthRateNegDD * (psi/phi) # Re-dimensionalised
            
        # Check Courant condition & stop the program if not satisfied
        if (np.abs(growthRatePos) * tdiff/rdiff) > 1 or (np.abs(growthRateNeg) * tdiff/rdiff) >1:
            print("Courant condition not satisfied.")
            print("Courant value for growthratepos " + str(np.abs(growthRatePos) * tdiff/rdiff))
            print("Courant value for growthrateneg " + str(np.abs(growthRateNeg) * tdiff/rdiff))
            sys.exit()        
        
        # Calculate new N & add to the array
        NNew = Nprime - (tdiff/rdiff) * ((growthRatePos*NHalfPos) - (growthRateNeg*NHalfNeg)) # Calculate new N, i.e. N for current timestep
        NArraysArray[1].append(NNew)

        # Calculate element of SS integral for current r & add to array
        SSIntegralElement = pow(r,3) * (NArraysArray[1][rCount] - NArraysArray[0][rCount])
        SSSumsArray.append(SSIntegralElement)
        #print("Q * r^3 = " + str(Q*pow(r,3)))
        #print("SS Element = " + str(SSIntegralElement))
    
    # Add mean radius of N to the array
    NArrayAvgR.append(np.average(rBins, weights = NArraysArray[1]))
    #print("Avg N radius = " + str(NArrayAvgR[timeCount]))
     
    # Delete N array for previous timestep (to conserve memory, since it is no longer needed)
    del NArraysArray[0]
    
    # Calculate SS value for this timestep & add to the array
    SS = SS - (Q * np.sum(SSSumsArray)) + ((A*PConc*tdiff/MInf)*np.exp(-PMEnergy/(R*Temp))) # Sum all integral element values and multiply by Q to approximate the integral (no precursor part to consider for growth only case)
    SSArray.append(SS)
    #print("SS " + str(SS))
    #print("   ")
    
    NNumArray.append(np.sum(NArraysArray[0]))

    
    """
    # Increase the temperature by the heating rate up to the maximum
    if (Temp < Tf):
        Temp += tdiff*HR
    """
    
# Create new time array to include t = 0
timeArrayFull = np.linspace(0, tmax, int(((tmax-tmin+tdiff)/tdiff) + 1))

# Plot the final distribution of nanoparticles
plt.plot(rBins, NArraysArray[0])
plt.title("Final nanoparticle size distribution")
#plt.ylim(0, 4.5e18)
plt.plot(rCrit,1e9,'ro')
plt.show()


# Plot the change in precursor population over time
plt.plot(timeArrayFull, PArray)
plt.title('Precursor concentration against time')
plt.show()

print(NArrayAvgR)
# Plot mean radius against time
plt.plot(timeArrayFull, NArrayAvgR)
plt.title("Mean radius against time")
plt.xscale('log')
plt.show()

# Plot supersaturation against time
plt.plot(timeArrayFull, SSArray)
plt.title("Supersaturation against time")
plt.xscale('log')
plt.yscale('log')
plt.show()

# Plot the total number of nanoparticles over time
plt.plot(timeArrayFull, NNumArray)
plt.title("Number of nanoparticles over time")
plt.xscale('log')
plt.show()
