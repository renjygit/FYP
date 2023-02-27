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
rdiff =  3e-11 #0.3 Angstroms
rmax = 8.1e-9
rmin = 0.5e-9
rBins = np.linspace(rmin, rmax,int((rmax-rmin+rdiff)/rdiff))

# Create empty array to hold growth rate values
growthRateListList = []

# Create supersaturation bins to iterate through
SSBinsList = [2., 10., 50., 100., 200.]

# Set gamma & Damkohler number (zeta) values
gamma = 0.2
zeta = 1

##########################################################################
##########################################################################

#right, let's sort out dimensionless units first for simple computations#

##########################################################################
##########################################################################

phi = (R*Temp)/(2*gamma*VmMolar)
psi = (phi**2)*D*VmMolar*MInf
#beta = r*phi
#tau = t*psi
print("phi = " + str(phi))
print("psi = " + str(psi))

SSCount = -1
#
for SS in SSBinsList:
    SSCount += 1
    growthRateListList.append([])
    for r in rBins:

        #sort out dimensionless rad
        beta = r*phi
        #print("beta = " + str(beta))

        #now sort out growth rate
        growthRate = (SS - np.exp(1/beta))/(beta + zeta)
        dimGrowthRate = growthRate*(psi/phi)

        growthRateListList[SSCount].append(dimGrowthRate)


# Plot growth rate against radius
for n in np.linspace(0,SSCount,(SSCount+1)):
    #fig = plt.figure()
    #plt.subplot(2, 1, 1)
    plt.plot(rBins, growthRateListList[int(n)], label = "Supersaturation = " + str(SSBinsList[int(n)]))

plt.xlabel("Radius")
plt.ylabel("Instantaneous Growth Rate")
plt.axhline(y=0, color='black', ls='--')
plt.legend()
plt.yscale('symlog')
plt.ylim(-1*10**(-6), 2*10**(-6))
plt.show()

# Produce second plot, for varying gamma values.
rmax = 8.1e-9
rmin = 0.2e-9
rBins = np.linspace(rmin, rmax,int((rmax-rmin+rdiff)/rdiff))
# Create empty array to hold growth rate values
growthRateListList = []

# Create gamma bins to iterate through
gammaBinsList = [0.075,0.1,0.15,0.2,0.5,1.]
# Set supersaturation and Damkohler number values
SS = 100
zeta = 1
gammaCount = -1
#
for gamma in gammaBinsList:
    gammaCount += 1
    
    phi = (R*Temp)/(2*gamma*VmMolar)
    psi = (phi**2)*D*VmMolar*MInf
    
    growthRateListList.append([])
    for r in rBins:
        
        #sort out dimensionless rad
        beta = r*phi

        #now sort out growth rate
        growthRate = (SS - np.exp(1/beta))/(beta + zeta)
        dimGrowthRate = growthRate*(psi/phi)
                
        growthRateListList[gammaCount].append(dimGrowthRate)

# Plot growth rate against radius
for n in np.linspace(0,gammaCount,(gammaCount+1)):
    plt.plot(rBins, growthRateListList[int(n)], label = "gamma = " + str(gammaBinsList[int(n)]))

plt.xlabel("Radius")
plt.ylabel("Instantaneous Growth Rate")
plt.axhline(y=0, color='black', ls='--')
plt.legend()
plt.yscale('symlog')
plt.ylim(-1*10**(-6), 2*10**(-6))
plt.show()
    
    
    

# Produce third plot, for varying Damkohler numbers.

# Set the radius bins
rmax = 8.1e-9
rmin = 0.5e-9
rBins = np.linspace(rmin, rmax,int((rmax-rmin+rdiff)/rdiff))

# Create empty array to hold growth rate values
growthRateListList = []

# Create Damkohler bins to iterate through
zetaBinsList = [0.001, 0.01,0.1,1,10,100,1000]
# Set supersaturation and gamma values
SS = 100
gamma = 0.2

phi = (R*Temp)/(2*gamma*VmMolar)
psi = (phi**2)*D*VmMolar*MInf

zetaCount = -1
#
for zeta in zetaBinsList:
    zetaCount += 1
    
    growthRateListList.append([])
    for r in rBins:

        #sort out dimensionless rad
        beta = r*phi

        #now sort out growth rate
        growthRate = (SS - np.exp(1/beta))/(beta + zeta)
        dimGrowthRate = growthRate*(psi/phi)
                
        growthRateListList[zetaCount].append(dimGrowthRate)


# Plot the growth rate against radius
for n in np.linspace(0,zetaCount,(zetaCount+1)):
    plt.plot(rBins, growthRateListList[int(n)], label = "zeta = " + str(zetaBinsList[int(n)]))

plt.xlabel("Radius")
plt.ylabel("Instantaneous Growth Rate")
plt.axhline(y=0, color='black', ls='--')
plt.legend()
plt.yscale('symlog')
plt.ylim(1*10**(-12), 3*10**(-6))
plt.show()
