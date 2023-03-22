import numpy as np
import matplotlib.pyplot as plt
import sys
import csv
import multiprocessing as mp


# Define the van Leer flux limiter function
def vanLeerFunc(X):
    return (X + np.abs(X)) / (1 + X)

def f(x):
    return x*x

def NucAndGrow(radius, rCount, NVal, NValTake1, NValPlus1, phi, zeta, psi, SS, tdiff, rdiff, radiusmin, radiusmax): # deal with rcount
    #print("NucAndGrow Start")
        
    # Calculate the intermediate value of N, i.e. the nucleated N (for this radius value)
    Nprime = NVal + 0 # This adds 0 for the growth only case, since there is no nucleation
    
    # Deltapos & deltaneg use number of nanoparticles at radii (i+1) & (i-1), which don't exist at the boundaries, so set NHalfPos and NHalfNeg to previous N for these edge cases
    if radius == radiusmin or radius == radiusmax:
        NHalfPos = NVal
        NHalfNeg = NVal
    else:
        # Calculate deltapos and deltaneg
        deltapos = NValPlus1 - NVal
        deltaneg = NVal - NValTake1
        
        # If one or both of the deltas = 0, one or both of their ratios will diverge, so include a condition to avoid this
        if deltapos == 0 or deltaneg == 0:
            NHalfPos = NVal # Set NHalfPos and NHalfNeg values equal to the previous N value
            NHalfNeg = NVal
        else:
            # Calculate NHalfPos and NHalfNeg using the van Leer function
            NHalfPos = NVal + 0.5*vanLeerFunc(deltapos/deltaneg)*deltaneg
            NHalfNeg = NVal - 0.5*vanLeerFunc(deltaneg/deltapos)*deltapos
                
    # Find beta values for radii values (i+1/2) and (i-1/2)
    betaPos = (radius+0.5*rdiff) * phi # de-dimensionalised radius
    betaNeg = (radius-0.5*rdiff) * phi # de-dimensionalised radius
    
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

    # Calculate element of SS integral for current r & add to array
    SSIntegralElement = pow(radius,3) * (NNew - NVal)
    #print("Q * r^3 = " + str(Q*pow(r,3)))
    #print("SS Element = " + str(SSIntegralElement))
    
    #print("NucAndGrow End")
    
    return NNew, SSIntegralElement


if __name__ == '__main__':
    
    # Define required constants
    R = 8.314
    kB = 1.38e-23
    NA = 6.022e23
    A = 4e3
    VmMolar = 3.2996e-5 # m^3 mol^-1 - molar volume of monomer
    rmMolar = (3/4*np.pi) * pow(VmMolar, 1/3)
    Vm = VmMolar / NA # (4/3)*np.pi*pow(rm,3)
    rm = pow((3*Vm) / (4*np.pi), (1/3)) # 1.87e-10
    D = 1e-12
    Vtot = 1e-4 #0.0252-3 # m^3
    MInf = 1e-2 #mol m^-3 #
    u = 0.42 # coagulation parameter
    gamma = 0.06 #0.4 # J m^-2 - surface energy
    k = 1e12 # J mol^-1 m^-1 - constant in the exponent of the distribution function
    kr = 1.6e-9
    density = 5.816e3 # kg m^-3
    MW = 0.09586 # kg mol^-1 - monomeric molecular weight
    Q = (4*np.pi*density) / (3 * MW * Vtot * MInf)
    

    # Define radius bins
    rdiff =  3e-11 # m
    rmax = 8.1e-9
    rmin = 4e-10
    rBins = np.linspace(rmin,rmax,int((rmax-rmin+rdiff)/rdiff))

    # Define timesteps to iterate through & create an array
    tmax = 1 #s
    tdiff = 0.6e-5
    tmin = tdiff
    timeArray = np.linspace(tdiff, tmax, int((tmax-tmin+tdiff)/tdiff))
    timeCount = 0 # Used to track number of time iterations
    

    
    # Define temperature variable & heating rate
    Temp = 503 # 453
    #Tf = 567
    #HR = 0.025 # K s^-1 (1.5 K min^-1)
    
    # Define supersaturation variable & array
    SS = 5e2
    SSArray = [SS]
    
    # Create an array of arrays to hold number of nuclei in each radius bin for different timesteps
    NArraysArray = [[]]
    
    # Give number of nanopartices initially
    NNumMole = 60e-6 # mol
    NNum = NNumMole * NA # number
    print("Initial number of nanoparticles = " + str(NNum))
    
    # Create the input distribution of nanoparticles
    sum = 0
    for r in rBins:
        g = (1 / (1e-10 * np.sqrt((2*np.pi)))) * np.exp(-0.5*(np.power((r - 1e-9)/(1e-10), 2))) # Normal distribution with mean radius 1e-9 m & standard deviation 1e-10 m
        gN = NNum * g * rdiff
        
        sum += gN # Calculated to check this sums to NNum (approximately)
        
        NArraysArray[0].append(gN)
    print("Sum of nanoparticles in initial distribution = " + str(sum))
    
    NNumArray = [np.sum(NArraysArray[0])]
    
    # Plot the initial distribution of nanoparticles
    plt.plot(rBins, NArraysArray[0])
    plt.title("Nanoparticles number against radius at t = 0")
    #plt.ylim(0, 4.5e18)
    plt.show()
    
    
    # Create an array to hold the mean radius values
    NArrayAvgR = [np.average(rBins, weights = NArraysArray[0])]
    print("Initial avg N = " + str(NArrayAvgR[0]))
    
    
    
    
    

    """
    list1 = [1, 2, 3, 4]
    print("before")
    with mp.Pool(5) as p:
        print(p.map(f, list1))
        print("after")
    """

    # Main loop steps forward in time
    for time in timeArray:
        timeCount += 1
        
        
        if timeCount % 100 == 0:
            print("time = " + str(time))
    
        #print("temp " +str(Temp))
        
        # Calculate the critical radius
        rCrit = (2*gamma*VmMolar) / (R*Temp*np.log(SS))
        
        # Calculate the de-dimensionalised variables
        phi = (R * Temp) / (2 * gamma *VmMolar)
        psi = phi**2 * D * VmMolar * MInf
        zeta = 1e-3 # (D * phi) / kr
        tau = time * psi
    
        # Add a new empty array to the N array of arrays to hold the newly calculated N values
        NArraysArray.append([0]* len(NArraysArray[0]))
        # Create empty arrays for the SS integral elements and the N values at radii values of (i+1/2) and (i-1/2) (reset each time iteration)
        SSSumsArray = [0] * len(NArraysArray[0])

        
          
        
        rCountArray = range(0,len(rBins)-1)
        
        NArraysArrayPlus1 = []
        NArraysArrayTake1 = []
        aCount = -1
        for a in NArraysArray[0]:
            aCount += 1
            if aCount != (len(NArraysArray[0])-1):
                NArraysArrayPlus1.append(NArraysArray[0][aCount+1])
            else:
                NArraysArrayPlus1.append(0)
            if aCount !=0:
                NArraysArrayTake1.append(NArraysArray[0][aCount-1])
            else:
                NArraysArrayTake1.append(0) # just holds the place of the element
            
        
        
        phiArray = [phi] * len(rBins)
        zetaArray = [zeta] * len(rBins)
        psiArray = [psi] * len(rBins)
        SSSingleArray = [SS] * len(rBins)
        tdiffArray = [tdiff] * len(rBins)
        rdiffArray = [rdiff] * len(rBins)
        rminArray = [rBins[0]] * len(rBins)
        rmaxArray = [rBins[len(rBins) - 1]] * len(rBins)
        
        with mp.Pool(40) as p:
            #print("before starmap")
            mapOut = p.starmap(NucAndGrow, zip(rBins,rCountArray, NArraysArray[0], NArraysArrayPlus1, NArraysArrayTake1, phiArray, zetaArray, psiArray, SSSingleArray, tdiffArray, rdiffArray, rminArray, rmaxArray))
            #print("after starmap")

            
        nCount = -1
        for n in mapOut:
            nCount += 1
            NArraysArray[1][nCount] = mapOut[nCount][0]
            SSSumsArray[nCount] = mapOut[nCount][1]

    
        
        # Add mean radius of N to the array
        NArrayAvgR.append(np.average(rBins, weights = NArraysArray[1]))
        #print("Avg N radius = " + str(NArrayAvgR[timeCount]))
         
        # Delete N array for previous timestep (to conserve memory, since it is no longer needed)
        del NArraysArray[0]
        
        # Calculate SS value for this timestep & add to the array
        SS -= (Q * np.sum(SSSumsArray)) # Sum all integral element values and multiply by Q to approximate the integral (no precursor part to consider for growth only case)
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
    plt.title("Final nanoparticle size distribution gamma = " + str(gamma))
    #plt.ylim(0, 4.)
    plt.savefig("Final_distrib_gamma_" + str(gamma))
    plt.show()
    
    with open("Final_distrib_gamma_" + str(gamma) + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(NArraysArray[0])
    
    print(NArrayAvgR)
    # Plot mean radius against time
    plt.plot(timeArrayFull, NArrayAvgR)
    plt.title("Mean radius against time gamma = " + str(gamma))
    plt.xscale('log')
    plt.savefig("NAvg_gamma_" + str(gamma))
    plt.show()
    
    with open("NAvg_gamma_" + str(gamma) + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(NArraysArray[0])
    
    # Plot supersaturation against time
    plt.plot(timeArrayFull, SSArray)
    plt.title("Supersaturation against time gamma = " + str(gamma))
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("Supersaturation_" + str(gamma))
    plt.show()
    
    with open("Supersaturation_" + str(gamma) + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(NArraysArray[0])
    
    # Plot the total number of nanoparticles over time
    plt.plot(timeArrayFull, NNumArray)
    plt.title("Number of nanoparticles over time gamma = " + str(gamma))
    plt.xscale('log')
    plt.savefig("NNum_gamma_" + str(gamma))
    plt.show()
    
    with open("NNum_gamma_" + str(gamma) + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(NArraysArray[0])
