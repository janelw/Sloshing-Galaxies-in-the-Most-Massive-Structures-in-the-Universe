
# coding: utf-8

# In[1]:


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
from scipy import stats


# In[2]:


#dfile = 'C://Users/Janel/Desktop/all_BCG_coords.fits'
#ftable = fits.open(dfile)
#data = Table(ftable[1].data)
#tfile = 'C://Users/Janel/Desktop/master_BCG_Cluster_Data.txt' 
tfile2 = 'C://Users/Janel/Desktop/master_BCG_Cluster_Data2.txt' #pulls master text file from directory 
#data = ascii.read(tfile,format = 'basic') # read that texr file into an ascii table 
data =ascii.read(tfile2,format = 'basic')  
ncl = len(unique(data, keys = 'Name'))



bcg_sz = data['BCG_Xp_AngDist'] #bcg-sz seperation 
bcg_peak = data['BCG_Xp_AngDist'] # bcg - xray peak separation 
bcg_xcen = data['BCG_Xc_AngDist'] #bcg - xray center separation 
sz_xcen = data['SZ_Xray_AngDist'] #sz - xcen separation 
spt_size = data['SPT_size'] #spt size 
spt_sn = data ['SPT_SN'] #spt - sn 
z_cls = data['Redshift'] #redshift of the galaxy clusters 
mass = data['m500']
masse = data['m500err']
beam = 0.0002*(180/math.pi)*3600 #beam value in radians (lamda/D -> radians (arcsec)) this is one value 

a = data['SZ_RA'] #here and below are indexed values within the table 
b = data['SZ_Dec']
c = data['Xcen_RA']
d = data['Xcen_Dec']
e = data['Xpeak_RA']
f = data['Xpeak_Dec']
g = data['BCG_RA']
h = data['BCG_Dec']

SZc = SkyCoord(a,b, unit = 'deg', frame ='fk5') #SZ center
XRc = SkyCoord(c,d, unit = 'deg', frame = 'fk5') #xray center
XRp = SkyCoord(e,f, unit = 'deg', frame ='fk5')#xray peak 
BCG = SkyCoord(g,h, unit = 'deg',frame = 'fk5')#bcg 

sep6 = SZc.separation(XRp) #SZ and XRAY Peak 
nsep6 = sep6.rad
asep6 = sep6.arcsec 


sep5 = SZc.separation(XRc)  #SZ and XRAY CENTER 
nsep5 = sep5.rad
asep5 = sep5.arcsec

sep7 = BCG.separation(XRp) #BCG and XRAY Peak
nsep7 = sep7.rad
asep7 = sep7.arcsec 



#tfile = 'C://Users/Janel/Desktop/all_BCG_coords.txt'
#bdata = ascii.read(tfile)


#bdata = Table(names=('Name', 'Redshift', 'Redshift_Err','N_Members','vd_biw','vd_err','SZ_RA','SZ_Dec','Xpeak_RA',
#'Xpeak_Dec','Xcen_RA','Xcen_Dec','BCG_RA','BCG_Dec','BCG_SZ_AngDist','BCG_Xp_AngDist','BCG_Xc_AngDist','SZ_Xray_AngDist','SZ_Xray_Arcsec')

cos = astropy.cosmology.FlatLambdaCDM(H0 = 70, Om0 = 0.3, Tcmb0 = 2.725) #defining cosmological values 

simu = [] #indexing simiuation values 
simuk = [] #simulation of 
random.seed()
for i in range(len(data)): #for the length of cdata
    szerr = math.sqrt(beam**2 + spt_size[i]**2 )/ spt_sn[i] #this is the expected error for the SZ cluster, it is delta theta in the equation
    #and sigma in application #
    rval = np.random.rand(1000) 
    x = np.sqrt(-2*np.log(rval)*szerr**2) #this is the simulated measurements. now i have to put this into an array 
    simu = np.append(simu,x)
    adj = cos.angular_diameter_distance(z_cls[i]).value 
    fac = (3600*(180/math.pi))
    xr = np.divide(x,fac)
    xmpc = np.multiply(xr,adj)
    xkpc = np.multiply(xmpc,1000)
    simuk = np.append(simuk,xkpc)
    
        
#print(np.median(simu))
#print(np.std(simu))        
#print(max(simu))
#print(min(simu))

simux = []
simux_k = []
#random.seed()
for i in range(len(data)): #for the length of cdata
    xerr = 1 #this is the expected error for the SZ cluster, it is delta theta in the equation
    #and sigma in application #
    rval1 = np.random.rand(1000)
    x = np.sqrt(-2*np.log(rval1)*xerr**2) #this is the simulated measurements. now i have to put this into an array 
    simux = np.append(simux,x)
  
    adj = cos.angular_diameter_distance(z_cls[i]).value 
    fac = (3600*(180/math.pi))
    xr = np.divide(x,fac)
    xmpc = np.multiply(xr,adj)
    xkpc = np.multiply(xmpc,1000)
    simux_k = np.append(simux_k,xkpc)
    


# In[3]:


#BCG and SZ in kpc
sep = np.array(bcg_sz)
sep1 = np.array(bcg_peak)

prange = [0, 400]
bins = 10
 
# plotting a histogram (_hist)
#plt.hist(sep, bins, prange, color = 'purple', histtype = 'bar', rwidth = 1,density = True)

#plotting pdf histogram (_pdf)
#plt.hist(sep, bins, prange, color = 'purple', histtype = 'bar', rwidth = 1, density = True)

#plotting the cumulative distribution (_cdf)
plt.hist(sep, bins, prange, color = 'cyan', histtype = 'step', rwidth = 2, cumulative = True, density = True, label = 'BCG - SZ')

print(np.median(sep))
print(np.median(sep1))
print(np.median(simuk))
print(np.percentile(sep1, .1, ))
print(np.percentile(sep, .1, ))
print(np.percentile(simuk, .1, ))
# x-axis label
plt.xlabel('Seperation Value (kpc)')
# frequency labe
plt.ylabel('Cumulitive Value (%)')
# plot title



prange = [0, 400]
bins = 100
#plt.hist(simuk, bins, prange, color = 'orange', histtype = 'step', rwidth = 1, density = True)
plt.hist(simuk, bins, prange, color = 'black', histtype = 'step', rwidth = 2, cumulative = True, density = True, label = 'Simulated from SZ Measurement Errors')
plt.hist(sep1, bins, prange, color = 'blue', histtype = 'step', rwidth = 2, cumulative = True, density = True, label = 'BCG - Xray')

plt.legend(loc='upper left', frameon=False )
  
plt.ylim(0, 1.5)    
plt.xlim(0, 400)
   


#plt.savefig('C://Users/Janel/Desktop/MultiGraph2.png')  
     #   values = lines.split(',')
      #  print(values)


# In[4]:


#BCG and XRAY PEAK in kpc 

sep1 = np.array(bcg_peak)

print(max(bcg_peak))
print(min(bcg_peak))
prange = [0, 400]
bins = 10
 
# plotting a histogram (_hist)
plt.hist(sep1, bins, prange, color = 'blue', histtype = 'bar', rwidth = 1, density = True)

#plotting pdf histogram (_pdf)
#plt.hist(sep1, bins, prange, color = 'blue', histtype = 'bar', rwidth = 1, density = True)

#plotting the cumulative distribution (_cdf)
#plt.hist(sep1, bins, prange, color = 'blue', histtype = 'bar', rwidth = 1, cumulative = True, density = True)


# x-axis label
plt.xlabel('Seperation Value (kpc)')
# frequency labe
plt.ylabel('N')
# plot title
plt.title('BCG and Xray Peak Offset')
#plt.show()
     #   values = lines.split(',')
      #  print(values)
        
#prange = [0, 20]
#bins = 100
#plt.hist(simux_k, bins, prange, color = 'orange', histtype = 'step', rwidth = 1, density = True)


#plt.hist(simux_k, bins, prange, color = 'orange', histtype = 'step', rwidth = 1, density = True)

#plt.savefig('C://Users/Janel/Desktop/BCG_Xrayp_Graph_kpc_hist.png')    


# In[7]:


#BCG and XRAY CENTER in kpc
sep2 = np.array(bcg_xcen)


prange = [0, 300]
bins = 8
 
# plotting a histogram
#plt.hist(sep2, bins, prange, color = 'green', histtype = 'bar', rwidth = 1)

#plotting pdf histogram (_pdf)
plt.hist(sep2, bins, prange, color = 'green', histtype = 'bar', rwidth = 1, density = True)

#plotting the cumulative distribution (_cdf)
#plt.hist(sep2, bins, prange, color = 'green', histtype = 'bar', rwidth = 1, cumulative = True, density = True)

# x-axis label
plt.xlabel('Seperation Value (kpc)')
# frequency labe
plt.ylabel('N')
# plot title
plt.title('BCG and Xray Center Offset')

        
prange = [0, 350]
bins = 100
plt.hist(simuk, bins, prange, color = 'orange', histtype = 'step', rwidth = 1, density = True)
  
#plt.savefig('C://Users/Janel/Desktop/BCG_Xcen_graph_kpc_hist.png')    


# In[8]:


#SZ and Xray center seperation in kpc 
sep3 = np.array(sz_xcen)

prange = [0, 350]
bins = 6
 
# plotting a histogram (_hist)
#plt.hist(sep3, bins, prange, color = 'pink', histtype = 'bar', rwidth = 1)

#plotting pdf histogram (_pdf)
plt.hist(sep3, bins, prange, color = 'pink', histtype = 'bar', rwidth = 1, density = True)

#plotting the cumulative distribution (_cdf)
#plt.hist(sep3, bins, prange, color = 'pink', histtype = 'bar', rwidth = 1, cumulative = True, density = True) 

# x-axis label
plt.xlabel('Seperation Value (kpc)')
# frequency labe
plt.ylabel('N')
# plot title
plt.title('SZ and Xray Center Offset')
     
prange = [0, 350]
bins = 100
plt.hist(simuk, bins, prange, color = 'orange', histtype = 'step', rwidth = 1, density = True) 

#plt.hist(simux, bins, prange, color = 'orange', histtype = 'step', rwidth = 1, density = True)
#plt.savefig('C://Users/Janel/Desktop/SZ_Xcen_Graph_kpc_hist.png')    


# In[9]:


#SZ and Xray Center in arcseconds 
prange = [0, 50]
bins = 8

#plotting a histogram (_hist)
#plt.hist(asep5, bins, prange, color = 'orange', histtype = 'bar', rwidth = 1)
 
#plotting pdf histogram (_pdf)
plt.hist(asep5, bins, prange, color = 'orange', histtype = 'bar', rwidth = 1, density = True)

#plotting the cumulative distribution (_cdf)
#plt.hist(asep5, bins, prange, color = 'Orange', histtype = 'bar', rwidth = 1, cumulative = True)    
    
    
# x-axis label
plt.xlabel('Seperation Value (Arcseconds)')
# frequency labe
plt.ylabel('N')
# plot title
plt.title('SZ and Xray Center Offset')



prange = [0, 45]
bins = 100
plt.hist(simu, bins, prange, color = 'black', histtype = 'step', rwidth = 1, density = True)
#plt.hist(simux_k, bins, prange, color = 'blue', histtype = 'step', rwidth = 1, density = True)

#plt.savefig('C://Users/Janel/Desktop/SZ_X_Arc_Graph_arc_hist.png')        


# In[10]:


#SZ and Xray Peak seperation in arcseconds 
prange = [0, 50]
bins = 7

#plotting a histogram (_hist)
#plt.hist(asep6, bins, prange, color = 'cyan', histtype = 'bar', rwidth = 1)

#plotting pdf histogram (_pdf)
plt.hist(asep6, bins, prange, color = 'cyan', histtype = 'bar', rwidth = 1, density = True)

#plotting the cumulative distribution (_cdf)
#plt.hist(asep6, bins, prange, color = 'cyan', histtype = 'bar', rwidth = 1, cumulative = True, density = True)


# x-axis label
plt.xlabel('Seperation Value (Arcseconds)')
# frequency labe
plt.ylabel('N')
# plot title
plt.title('SZ AND Xray Peak Offsets')



prange = [0, 50]
bins = 85
plt.hist(simu, bins, prange, color = 'orange', histtype = 'step', rwidth = 1, density = True)
 
#plt.savefig('C://Users/Janel/Desktop/SZ_Xp_Arc_Graph_arc_hist.png') 


# In[11]:


sep1 = np.array(bcg_peak)

print(max(asep7))
print(min(asep7))
prange = [0, 90]
bins = 50
 
# plotting a histogram (_hist)

plt.hist(asep7, bins, prange, color = 'blue', histtype = 'step', rwidth = 1, cumulative = True, density = True)
plt.hist(sep, bins, prange, color = 'purple', histtype = 'step', rwidth = 1, cumulative = True, density = True)
#plotting pdf histogram (_pdf)
#plt.hist(sep1, bins, prange, color = 'blue', histtype = 'bar', rwidth = 1, density = True)

#plotting the cumulative distribution (_cdf)
#plt.hist(sep1, bins, prange, color = 'blue', histtype = 'bar', rwidth = 1, cumulative = True, density = True)


# x-axis label
plt.xlabel('Seperation Value (arcseconds)')
# frequency labe
plt.ylabel('Number of Galaxies')
# plot title
plt.title('BCG and Xray peak offset')



bins = 10
plt.hist(simux, bins, prange, color = 'orange', histtype = 'step', rwidth = 1, cumulative = True, density = True)

#plt.savefig('C://Users/Janel/Desktop/BCG_Xrayp_Graph_arc_hist.png')  


# In[14]:


prange = [0, 1.3]
bins = 10

# plotting a histogram (_hist)
plt.hist(z_cls, bins, prange, color = 'blue', histtype = 'step', rwidth = 1, density = True)


# x-axis label
plt.xlabel('Redshift')
# frequency labe
plt.ylabel('Number of Clusters')
# plot title
print(np.average(z_cls))
plt.show()


mass = data['m500']
masse = data['m500err']


#plt.scatter(z_cls,mass,color='G')
#plt.scatter(x,y, z_cls, color = 'b')
#plt.errorbar(z_cls, mass, yerr=masse, fmt = 'bo')
#plt.xlabel('Redshift')
#plt.ylabel('Mass')


plt.savefig('C://Users/Janel/Desktop/Redshift_1.png')
#plt.show()
     #   values = lines.split(',')
      #  print(values)


# In[12]:


mass = data['m500']
masse = data['m500err']

print(mass)
plt.scatter(z_cls,mass,color='G')
#plt.scatter(x,y, z_cls, color = 'b')
plt.errorbar(z_cls, mass, yerr=masse, fmt = 'bo')
plt.xlabel('Redshift')
plt.ylabel('Mass')

plt.savefig('C://Users/Janel/Desktop/Mass1_Graph.png')
#spt is well defined in mass, high mass very rare 


# In[13]:


# old code:
# with open('C:\\Users\Janel\Desktop\BCG-Xray_peak_separation.txt') as f: #opens seperations file 
#    sep = f.readlines() #reads into seperations file
#    sep1 = np.array(sep)


#put titles on the graphs *****
#SZ is large so 
#purple in arcsec 
#correlations between orange and purple graph signify what the uncertainty
#account for the noise of the SZ which is accounting for the 
#may want to make a cumuitive distribution for seperatoins (like in paper) and a theoretical plot like in the paper*****
#page 18 - reread section of the paper on figures 
#most noises are gaussian 
#delta theta is sigma for our curve 
#delta theta will often be 20 arcsec 
#need to take the integral to find the probability (true probability)
#probability density functions: divide by the total number of clusters then you get the fraction of cluster in each bin *****

#new versions of graphs that are probability 


# In[14]:


#make plots of the mass, mass error, and make legends, decide on how to handle ones with 2 bcgs vs redshift (do this first)

#for 2 bcgs: 1). choose one thing, 2) average the positions(the coordinates) 3) average the seperations(average the offsetts) <-- talk to mike
#the avg offset will be smaller
#average the coordinates, and then average the seperation values 
#make new version of the plots based off of those 
#15 clusters have 2 BCGs 

