import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astroplan import Observer
import warnings
import time
import gc

#user defined file
import simulate_dp24

n= 50 #number of repeating FRBs
ndatasets= 100 #Number of Datasets to be generated, saved, and analysed
ndays= 2*365 #Number of observation days for a given dataset

R= Time('2000-01-01T12:00:00') #reference epoch
chime= Observer.at_site("CHIME") #CHIME location
lat= chime.location.lat.rad #CHIME Latitude in radian

#Beam ellipse, x and y are points to check
a= 2.5/2 #semi-major axis
b= 120/2 #Semi-minor axis
h= 180 #center point az
k= 90 #center point alt
h2= 0
k2= 90
h3= 360
k3= 90

for nd in range(190, 190+ndatasets): #range of seed values 

    st = time.time()
    print("Dataset: %d"%(nd+1))
    data=simulate_dp24.data(n, nd+1) #n= number of FRBs, i+1 = Seed for random numbers

    data.drop(data[data['Dec'] <= -11].index, inplace=True) #below declination −11 is outside the CHIME/FRB field of view
    data= data.reset_index(drop=True) #reindexing

    RA= data.RA
    Dec= data.Dec
    phs= data.phase #start phase
    P= data.Period
    D= data.Duty_cycle

    obs_frbs= [] #Observed FRBs
    ent= []
    ext= []
    new_phs=[]

    for i in range(ndays*288): #Every 5 mins
        
        #print("Min:",i*5)
        T= Time('2018-01-01T12:00:00') #time T
        T= T+i*(5*u.min) 
        
        P= data.Period
        
        new_phs= ((T.mjd-R.mjd)%P)/P
        
        x1= new_phs<=0.5+0.5*D
        x2= new_phs>=0.5-0.5*D
        x= x1 & x2
        Act_frbs= data[x] #Active FRBs
        Act_frbs= Act_frbs.reset_index(drop=True)
           
        start= T
        end= T+1*(5*u.min)
        twin= start + (end - start) * np.linspace (0, 1, 5) # Observing every minute
        
        ra= Act_frbs.RA
        dec= np.deg2rad(Act_frbs.Dec) #in radians
        lst= chime.local_sidereal_time(twin.value).deg
        my_saa= []

        #Converting Equatorial (RA-Dec) into horizontal (Alt-Az) coordinates for CHIME
        for j in range(len(ra)): #loop to calculate for each source

            ha= np.deg2rad(lst-ra[j])

            #if HA value is in negative, add 360 to it
            ha[np.where(ha<=0)]+= 2*np.pi

            alt= np.arcsin(np.sin(dec[j])*np.sin(lat) + np.cos(dec[j])*np.cos(lat)*np.cos(ha))
            az= np.arccos((np.sin(dec[j]) - np.sin(alt)*np.sin(lat))/(np.cos(alt)*np.cos(lat)))   
        
            alt= np.rad2deg(alt)
            az= np.rad2deg(az)
        
            #if sin(HA) is positve, subtract az from 360
            m= np.sin(ha)>0
            az[m]= 360-az[m]

            my_saa.append([alt, az])    

        #Determining active FRBs passing through beam region, and when (entry and exit time) 
        for m in range(len(my_saa)):
        
            x= my_saa[m][1]
            y= my_saa[m][0]
            r= (pow((x - h), 2) / pow(a, 2)) + (pow((y - k), 2) / pow(b, 2)); 
            #r<1 inside, r=1 on the ellipse, r>1 outside
            r= r<=1 
            r2 = (pow((x - h2), 2) / pow(a, 2)) + (pow((y - k2), 2) / pow(b, 2)); 
            r2= r2<=1
            r3 = (pow((x - h3), 2) / pow(a, 2)) + (pow((y - k3), 2) / pow(b, 2)); 
            r3= r3<=1
            r2+= r3
            r+= r2
            if r.sum()>=1: #checking if at any point passes through beam region
                dr= np.diff(r, prepend=0)
                sp= np.where(dr!=0)[0]
                if len(sp)%2!=0:
                    sp= np.append(sp, sp[-1]+1)
                    sp= sp.reshape((-1, 2))
                    [ent.append( start+ ((end - start) * (l/len(twin))) ) for l in sp[:,0]]
                    [ext.append( start+ ((end - start) * (l/len(twin))) ) for l in sp[:,1]]
                    for l in range(len(sp)): obs_frbs.append(Act_frbs.iloc[m])
                else:
                    sp= sp.reshape((-1, 2))
                    [ent.append( start+ ((end - start) * (l/len(twin))) ) for l in sp[:,0]]
                    [ext.append( start+ ((end - start) * (l/len(twin))) ) for l in sp[:,1]]
                    for l in range(len(sp)): obs_frbs.append(Act_frbs.iloc[m])#indexing is same for Act_frbs

    obs_frbs= pd.DataFrame(obs_frbs)
    obs_frbs= obs_frbs.reset_index(drop=True)
    #Saving all observed frbs in file
    #outdir = 'E:/FRB Project Data/Observed_FRBs_Vanilla_CHIME_DP24'
    outdir = "Observed_FRBs_Vanilla_CHIME_DP24"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    st2 = time.time()
    #Adding entry and exit times into dataframe obs_frbs
    obs_frbs= pd.concat([obs_frbs, pd.DataFrame(ent, columns=["Entry_time"]), pd.DataFrame(ext, columns=["Exit_time"])], axis = 1)
    
    #df.to_csv(f"{outdir}/observed_frbs(%s).csv"%(nd+1)) #saving as csv file
    np.savetxt(f"{outdir}/observed_frbs(%s).csv"%(nd+1), obs_frbs.values, fmt='%s,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%s', header=','.join(obs_frbs.columns), comments='')
    print("Observed data file saved successfully")

    print(f'For saving data: {time.time() - st2} seconds')
    print(f'For One 2 years run: {time.time() - st} seconds')

    del obs_frbs
    del my_saa
    gc.collect()
##################################################################################################
#Piece of code to convert string values of entry and exit time to Astropy time
# from astropy.time import Time
# Time(read_data.Exit_time.tolist(), format='isot', scale='utc')
    
