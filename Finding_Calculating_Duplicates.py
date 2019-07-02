
# coding: utf-8

# In[152]:


import numpy as np
import astropy
import regions
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits 
from astropy import units as u 
from regions import read_ds9, write_ds9
from astropy.coordinates import SkyCoord
import glob, os
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion
from astropy.coordinates import SkyOffsetFrame
from astropy import cosmology 
import math
from astropy.table import Table, Column, unique
from matplotlib import pyplot
import random
import scipy
import collections 
from collections import Counter


# In[193]:


file = 'C://Users/Janel/Desktop/master_BCG_Cluster_Data2.txt' #master cluster data
file2 = 'C://Users/Janel/Documents/Duplicates.txt' #names of duplicates that I found in previous code 
file3 = 'C://Users/Janel/Desktop/all_BCG_coords.txt'#original BCG coordinates 


outfil = 'C://Users/Janel/Desktop/Multiple_BCGs.txt'#writing info of multiple BCGs in this 

data =ascii.read(file,format = 'basic') #opening master file for reading 
data2 = ascii.read(file3) #opening BCG coordinaates for reading 
dup = open(file2, 'r') #opening duplicates file 



newdata = Table(names=('Name','SZ_RA','SZ_Dec','Xpeak_RA','Xpeak_Dec', 'BCG_RA', 'BCG_Dec'), dtype=('U17','f8','f8','f8','f8','f8','f8'))


 
cnames = data['Name']
szra = data['SZ_RA']
szdec = data['SZ_Dec']
xra = data['Xpeak_RA']
xdec = data['Xpeak_Dec']
bra = data ['BCG_RA']
bdec = data ['BCG_Dec']


doubles = []
for i in range(len(data)): 
    doubles = Counter(cnames).most_common() #identifying and counting the duplicate data 

for lines in dup:
    dup_names1 = lines.split()
    dup_names = '/n'.join(dup_names1) #getting the names for the duplicates 
    for i in range(len(data)): #for the length og data 
        if cnames[i] == dup_names: #if cnames matches dup_name 
            newdata.add_row((cnames[i], szra[i], szdec[i], xra[i], xdec[i],bra[i],bdec[i])) #write the data into the new file 

#print([19:21])
#newdata.write(outfil,format='ascii',overwrite=True)

#cluster names of doubles were copied and pasted into a new text document called, "duplicates"


# In[201]:


file4 = 'C://Users/Janel/Documents/Midpoint_Coordinates.txt'
file5 = 'C://Users/Janel/Desktop/urop_spt_clusters.txt'
file6 = 'C://Ussers/Janel/Documents/Average_Separations'
data2 = ascii.read(file4)
data3 = ascii.read(file5)

#Midpoint data:
RA = data2['RA_mp']
Dec = data2['Dec_mp']
cra = data3['col2']            
cdec = data3['col3'] 


c1 = SkyCoord(RA, Dec, unit='deg', frame = 'fk5')
c2 = SkyCoord(cra,cdec, unit='deg', frame = 'fk5')
#print(c2)
#sep1 = c1.separation(c2)  #SZ and XRAY CENTER 
#nsep1 = sep1.rad          #converting to radians 
#asep1 = sep1.arcsec       #converting to arcseconds 
#sin1 = math.sin(nsep1)    #taking the sine using the math module in order to start calculations 
#distance1 = np.multiply(sin1,adj) #multiplying the value above by the defined parameter , "adj"
#ndistance1 = abs(np.multiply(distance1,1000))
#print(ndistance1)


# In[ ]:


#write into a list 
#list of chosen BCGs -> make a new list 
#new file where every cluster has 1 bcg 
#measure the seperation of each bcg from SZ and xray
    #average seperations -> 50/50
#midpoint of the two coordinates of the BCGs and that will be the new location of the BCG 
    #average the coordinates 
    #if the distance is 0 it tells us that the BCG is seperated from the center 
#make a seperate file for each 
#should have a similar format for all 

