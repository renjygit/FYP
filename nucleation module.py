import numpy as np
import matplotlib.pyplot as plt
import statistics as stats

# Define required constants
R = 8.314
VmMolar = 1.65e-5 # m^3 mol^-1, molar volume of monomer
rmMolar = (3/4*np.pi) * pow(VmMolar, 1/3) # effective radius of monomer (from the molar volume)
rm = 1.87e-10
Vm = (4/3)*np.pi*pow(rm,3)
D = 1e-10
NA = 6.022e23
MInf = 1e-3 #mol m^-3 #
u = 0.46 # coagulation parameter
gamma = 0.6 # J m^-2
kB = 1.38e-23


"""
# Define time bins
tmax = 100 # 3600
tdiff = 0.01 #1e-3
timeArray = np.linspace(tdiff, tmax, int(tmax/tdiff))
timeCount = 0


# Define radius bins
rdiff =  3e-10 #3 Angstroms
rmax = 8.1e-9
rBins = np.linspace(0,rmax,int((rmax+rdiff)/rdiff))

# Define temp variable & heating rate
Temp = 453
Tf = 567
HR = 0.025 # K s^-1 (1.5 K min^-1)

# Define supersaturation
SS = 5e3

# Define nuclei arrays
NArraysArray = [[0] * len(rBins)]
NArrayAvgR = [stats.fmean(NArraysArray[0])]

# Main loop steps forward in time
for time in timeArray: #start at timestep tdiff not 0
    timeCount +=1 # Keep track of number of time iterations
    
    #print(timeCount)
    
    # Increase the temperature by the heating rate up to the maximum
    if (Temp < Tf):
        Temp += tdiff*HR
    #print("temp " +str(Temp))
    
    # Calculate critical radius, 'p' exponent and thus nucleation rate.
    rCrit = (2*gamma*VmMolar) / (R*Temp*np.log(SS))
    #print("rCrit " + str(rCrit))
    
    p = np.power((rCrit/rm), 3) # Need actual value of rm not molar value
    
    Rnuc = 8*np.pi*rm*D*NA * pow(SS, (p*u+1)) * pow(MInf, 2) * np.exp((-4*np.pi*pow(rCrit, 2)*gamma)/(3*kB*Temp)) #Need actual value of rm not molar value
    #lnRnuc = np.log(8*np.pi*rm*D*NA) + (p*u+1) * np.log(SS) * 2*np.log(MInf)  - ((4*np.pi*np.power(rCrit, 2)*gamma)/(3*kB*Temp)) #Need actual value of rm not molar value    


"""




# Produce 2 plots of nucleation rate against supersaturation, for varying temperature and gamma values.

# Create supersaturation bins to iterate through
SSBinsList = np.linspace(1, 20000, 20000)
    
# Produce first plot, for varying temperatures:

# Create temperature bins to iterate through    
tempBinsList = [373,423,473,523,573]

# Create empty arrays to hold values to plot later
RNucListList = []
pListList = []

# Choose a gamma value
gamma = 0.6 # J m^-2

tempCount=-1
# Produce nucleation rates for different supersaturations at several temperatures (given by tempBinsList) 
for Temp in tempBinsList:
    tempCount +=1
    RNucListList.append([])
    pListList.append([])
    for SS in SSBinsList: #start at timestep tdiff not 0
        
        # Calculate critical radius, 'p' exponent and thus nucleation rate.
        rCrit = (2*gamma*VmMolar) / (R*Temp*np.log(SS))
        
        p = np.power((rCrit/rm), 3) # Need actual value of rm not molar value
        pListList[tempCount].append(p)
        
        Rnuc = 8*np.pi*rm*D*NA * pow(SS, (p*u+1)) * pow(MInf, 2) * np.exp((-4*np.pi*pow(rCrit, 2)*gamma)/(3*kB*Temp)) #Need actual value of rm not molar value
        
        #print("SS^(pu+1) " + str(pow(SS, (p*u+1))))
        #print("exp((-4*np.pi*pow(rCrit, 2)*gamma)/(3*kB*Temp)) " + str(np.exp((-4*np.pi*pow(rCrit, 2)*gamma)/(3*kB*Temp))))
        
        RNucListList[tempCount].append(Rnuc*NA)
        
        #lnRnuc = np.log(8*np.pi*rm*D*NA) + (p*u+1) * np.log(SS) * 2*np.log(MInf)  - ((4*np.pi*np.power(rCrit, 2)*gamma)/(3*kB*Temp)) #Need actual value of rm not molar value

# Plot nucleation rate against supersaturation for each temperature in tempBinsList
# Also plot p against supersaturation for each temperature in tempBinsList
for n in np.linspace(0,tempCount,(tempCount+1)):
    plt.plot(SSBinsList, RNucListList[int(n)], label="Temperature = " + str(tempBinsList[int(n)]))
    
    #plt.plot(SSBinsList, pListList[int(n)], label="Temperature = " + str(tempBinsList[int(n)]))
plt.xlabel("Supersaturation")
plt.ylabel("Nucleation rate or p")
plt.yscale('log')
plt.legend()
plt.show()





# Produce second plot, for varying gamma values.

# Create gamma bins to iterate through
gammaBinsList = [0.2,0.4,0.6,0.8,1.0]

# Create empty arrays to hold values to plot later
RNucListList = []
pListList = []

# Choose a temperature value
Temp = 573

gammaCount=-1
# Produce nucleation rates for different supersaturations at several gamma values (given by gammaBinsList) 
for gamma in gammaBinsList:
    gammaCount +=1
    RNucListList.append([])
    pListList.append([])
    for SS in SSBinsList: #start at timestep tdiff not 0
        
        # Calculate critical radius, 'p' exponent and thus nucleation rate.
        rCrit = (2*gamma*VmMolar) / (R*Temp*np.log(SS))
        
        p = np.power((rCrit/rm), 3) # Need actual value of rm not molar value
        pListList[gammaCount].append(p)
        
        Rnuc = 8*np.pi*rm*D*NA * pow(SS, (p*u+1)) * pow(MInf, 2) * np.exp((-4*np.pi*pow(rCrit, 2)*gamma)/(3*kB*Temp)) #Need actual value of rm not molar value
        
        #print("SS^(pu+1) " + str(pow(SS, (p*u+1))))
        #print("exp((-4*np.pi*pow(rCrit, 2)*gamma)/(3*kB*Temp)) " + str(np.exp((-4*np.pi*pow(rCrit, 2)*gamma)/(3*kB*Temp))))
        
        RNucListList[gammaCount].append(Rnuc*NA)
        
        #lnRnuc = np.log(8*np.pi*rm*D*NA) + (p*u+1) * np.log(SS) * 2*np.log(MInf)  - ((4*np.pi*np.power(rCrit, 2)*gamma)/(3*kB*Temp)) #Need actual value of rm not molar value

# Plot nucleation rate against supersaturation for each gamma value in tempBinsList
# Also plot p against supersaturation for each gamma value in tempBinsList
for n in np.linspace(0,gammaCount,(gammaCount+1)):
    plt.plot(SSBinsList, RNucListList[int(n)], label="gamma = " + str(gammaBinsList[int(n)]))
    
    #plt.plot(SSBinsList, pListList[int(n)], label="gamma " + str(gammaBinsList[int(n)]))
plt.xlabel("Supersaturation")
plt.ylabel("Nucleation rate or p")
plt.yscale('log')
plt.legend()

plt.show()
