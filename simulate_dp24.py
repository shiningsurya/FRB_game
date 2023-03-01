import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
import os
##########################

def data(num, seed):

    n=num #Number of sources/data points to be generated
    s=seed  #seed value for random
    ##################################

    #Generating various parameters
    #print("Simulaiting data...")
    np.random.seed(s)
    RA=np.random.uniform(0, 359.999, n)
    Dec=np.random.uniform(-90, 89.999, n)
    period=np.random.uniform(0.003472, 1, n) #period in days
    d_cycle=np.random.uniform(0.01,0.499, n) #duty cycle  
    s_phase=np.random.uniform(0, 0.999, n) #start phase
    name=[("sample_%d_%d"%(s,i)) for i in range(n)]
    ##############################

    #Coverting into DataFrame and saving as CSV file
    #print("Saving data...")
    val = list(zip(name,RA, Dec,period,d_cycle,s_phase)) 
    df = pd.DataFrame(val, columns=['Name','RA', 'Dec','Period','Duty_cycle','phase'])
    outdir = 'simulated_data_CHIME_DP24'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    #df.to_csv("simulated_frb(%d,%d).csv"%(n,s))
    df.to_csv(f"{outdir}/simulated_frbs(%d,%d).csv"%(n,s))
    print("Data file saved successfully")
    return(df)
###############################

def plot(data):
    #plotting frbs ins sky
    resp= input("Do you want to visualize distribution of simulated FRB sources in sky? [Y/N]: ")
    resp=resp.upper()
    if resp=='Y':

        RA=data.RA
        Dec=data.Dec
        
        plt.scatter(RA,Dec,s=0.5)
        plt.title("FRBs in sky")
        plt.xlabel('RA')
        plt.xticks(range(0, 360,60))
        plt.ylabel('Dec')
        plt.yticks(range(-90, 90,20))
        plt.grid()
        plt.show()

        RA_rad=[]
        Dec_rad=[]
        RA[RA > 180] = RA[RA > 180] - 360
        
        #degress into radians, required for projection
        for i in range(len(RA)):
            RA_rad.append(np.deg2rad(RA[i]))
            Dec_rad.append(np.deg2rad(Dec[i]))

        fig = plt.figure(figsize = (6,6))
        ax = fig.add_subplot(111, projection='mollweide')
        ax.scatter(RA_rad, Dec_rad, s=1, marker='o', color='b')
        ax.grid(True) 
        plt.title('FRBs in sky')
        plt.xlabel('RA')
        plt.ylabel('Dec')
        plt.show()
        
    else:
        print("Invalid input")
