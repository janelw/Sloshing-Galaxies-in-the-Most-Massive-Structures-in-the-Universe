
# coding: utf-8

# In[2]:


#This code was the primary part of my research at MIT Kavli Institute of Astrophysics and Space Research 
#In this code, I look into the files and I coded the data inside to be readable by Python 
#Text with their own comment line of code are descriptions of the sections of the following section of data. this also includes 
#subtext of data that's indented immediately under some lines of code, such as when I define cosmological parameters 
#code with comment lines next to them are descriptions of the code within that line 


# In[2]:


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
from astropy.table import Table, Column
import collections


# In[3]:


#in order to use the data in the files you must first define them into variables.

#defining files in the directory to be used in the code
filename = 'C://Users/Janel/Desktop/urop_spt_clusters.txt'
filename1 = 'C://Users/Janel/Desktop/mbb_win_compat_2'
filename2 = 'C://Users/Janel/Desktop/BCG_xray_mod.txt'
#defining files to be written into 
ofile = 'C://Users/Janel/Desktop/All_separation_outputs.txt'   #master file with all of the values we calulated and the numbers they are derived from (old) 
outf1 = 'C://Users/Janel/Desktop/BCG-Xray_peak_separation.txt' #just the bcg-xray peak sep numbers 
outf2 = 'C://Users/Janel/Desktop/BCG-Xray_cent_sed_clusters.txt' #matched files for writing
boutfil = 'C://Users/Janel/Desktop/all_BCG_coords.txt'#master file with all the values we calculated and the numbers they are derived from (new)
coutfits= 'C://Users/Janel/Desktop/SZ_Xray_matched_clusters.fitparation.txt' #just the BCG-xray center sep numbers - although we wont be using this file 
#defining more files to be written intno 
coutfil = 'C://Users/Janel/Desktop/SZ_Xray_matches'#fits version of the above file (although we did not end up using it, here for back up in case the other text file crashes)
boutfits = 'C://Users/Janel/Desktop/all_BCG_coords.fits'#fits version of the above file ''''''''''


#my code, for some reason, could not find the file under, "filename1" in my directory, so I created a path below 
os.chdir(filename1) #creates a path for the file to be read into windows format??


#below is the data in, "filename" being read into an ascii file 
cdata = ascii.read(filename)   #this reads the context of filename into a "Table" object which is apparently a thing in python
                               #Look up Data Tables (astropy.table) for how they work, see: # docs.astropy.org/en/stable/table/index.html

#the following are data from the chart that has been indexed into varaibles
#this allows us to use each type of variable independently 
#using brackets lets you index columns by column  

cnames1 = cdata['col1']        #cluster names 
szra = cdata['col2']           #SZ right ascension 
szdec = cdata['col3']          #SZ declination 
Nmem = cdata['col4']           #Number of cluster members 
z_cls = cdata['col5']          #redshift 
z_cl_err = cdata['col6']       #redshift error
vd_biw = cdata['col7']         #veolcity?
vd_err = cdata['col8']         #velocity error 

#below are the beginnings of a table being made. Here, I'm just creating the header names for my data 
mdata = Table(names=('Name', 'Redshift', 'Redshift_Err','N_Members','vd_biw','vd_err','SZ_RA','SZ_Dec','Xpeak_RA','Xpeak_Dec','Xcen_RA','Xcen_Dec'), dtype=('U17','f8','f8','i8','f8','f8','f8','f8','f8','f8','f8','f8'))          
                               #dtype: formatting based on amount of memory needed to process the numbers. 
                               #32 vs 64 bit = smaller vs larger numbers 
                               #cof1.write(cnames1[i],z_cls[i],z_cl_err[i],Nmem[i],vd_biw[i],vd_err[i],szra[i],szdec[i],c4.ra.deg[0],c4.dec.deg[0],c3.ra.deg[0],c3.dec.deg[0])
bdata = Table(names=('Name', 'Redshift', 'Redshift_Err','N_Members','vd_biw','vd_err','SZ_RA','SZ_Dec','Xpeak_RA','Xpeak_Dec','Xcen_RA','Xcen_Dec','BCG_RA','BCG_Dec','BCG_SZ_AngDist','BCG_Xp_AngDist','BCG_Xc_AngDist','SZ_Xray_AngDist'), dtype=('U17','f8','f8','i8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
                               #above are the headers for the masterfile 
                               #using brackets lets you index row by row
                               #using brackets lets you index row by row

#below is when I actually start to read the files 
f = open(filename,'r')         #opens urop_spt_clusters.txt file and 'r' reads it 
f2 = open(filename2,'r')       #opens up BCG_xray_mod.txt file and 'r' reads it 

#this opens that files 
ofil = open(ofile,'w')         #opens up ofile for writing 
outf = open(outf1, 'a')        #opens up outf1 for writing and appending values 
bof1 = open(boutfil,'a')       #opens up boutfil for writing 
                               # print a header portion of the file with info about the columns:
hdrstring = '# Cluster_Name Cluster_z BCG_cnt  SZ-Xcen            BCG-Xcen            BCG-SZ         BCG-Xpeak'+'\n' #headers for seperations masterfile 
hdrstring2 = '# BCG-Xpeak'+'\n'#header for bcg_xraypeak seperations 
ofil.write(hdrstring)          #writes headers into file 
outf.write(hdrstring2)         #writes headers into file 

#here, I'm defining my cosomology 
cos = astropy.cosmology.FlatLambdaCDM(H0 = 70, Om0 = 0.3, Tcmb0 = 2.725)  #defining cosmology                                                               
                               #flat = no geometric cuvature, CDM = cold dark matter
                               #lambda = cosmological constant/normalization of the energy vacuum content of the universe
                               #H0 = hubble constant, Om0 = omega matter, normalized energy content of the universe
                               #Lambda and Om0 have to add to 1
                               #Tcmb0 = tempurature of the CMB
                               # We could now, for example, compute and print out the age of the universe at a redshift of 1 or 10 with:
                               #old cluster sky coordinate c1 = SkyCoord(cra, cdec, unit='deg', frame = 'fk5')
adj = cos.angular_diameter_distance(z_cls).value #distance that relates angular size to physical size in Mpc


#here is the beginning of my loop that I use to actually calculate my separations 
xraycount = 0
for line in f2:
    xraycount = xraycount+1        #counts the number of 
    columns1 = line.split()        #splits file into columns 
    xray_name1 = columns1[0:1]     #indexing row by row to get values 
    xray_peak_ra1 = columns1[1:2]  #xray peak right ascension 
    xray_peak_dec1 = columns1[2:3] #xray peal declination 
    xray_centr_ra1 = columns1[3:4] 
    xray_centr_dec1 = columns1[4:5] 
    #print(xray_name1)
            
                                # we want to drop the outer loop that reads in, line by line, the urop_clusters_whatever.txt file, and instead 
                                # read that file in as a Table (w/ ascii.read). We then need to match the cluster name in the X-ray file to the 
                                # cluster names in the urop_clusters_whatever.txt file, and where there is a match grab the appropriate redshift 
                                # and SPT SZ coordinates 
    xray_name = '/n'.join(xray_name1)
    print(xray_name)
    c4 = SkyCoord(xray_peak_ra1,xray_peak_dec1, unit = 'deg', frame ='fk5')
    c3 = SkyCoord(xray_centr_ra1,xray_centr_dec1, unit = 'deg',frame = 'fk5')
                                # Here is where we want to match xray_name with the cname array from the urop_spt_clusters.txt 
                                # file, and whatever index has that match, we use it to grab the relevant SPT SZ cluster data, 
                                # such as the redshift and the SZ coordinates (RA and Dec)
#Old loop    
#    for line in cdata:
#        if xray_name in line:
#            cnames1 = line['col1']
#            szra = line['col2']
#            szdec = line['col3']
#            Nmem = line['col4']
#            z_cls = line['col5']
#            z_cl_err = line['col6']
#            vd_biw = line['col7']
#            vd_err = line['col8']
    
#the loop below reads the length of cdata and matacheds the cluster names to the xray names. 
#Relevant data matches print out under "if"    
    for i in range(len(cdata)): #for the length of cdata
        if cnames1[i] == xray_name: #find matches within each of the files 
            mdata.add_row((cnames1[i],z_cls[i],z_cl_err[i],Nmem[i],vd_biw[i],vd_err[i],szra[i],szdec[i],c4.ra.deg[0],c4.dec.deg[0],c3.ra.deg[0],c3.dec.deg[0]))
            #create rows for matching clusters with data 
            adj = cos.angular_diameter_distance(z_cls[i]).value         #Prints angular distance using matching redshift values in Mpc
            c1 = SkyCoord(szra[i],szdec[i],unit = 'deg', frame = 'fk5') #SZ right ascension and declination 

            for file in glob.glob(cnames1[i] +'*_bcgs*'): #reads bcg file in order according to the cnames1 file file 
                bcgreg = read_ds9(file)
                count3 = 0
              
                for region in bcgreg: #this reads the regions in the data within bcgreg
                    cen = region.center  #defining the centers
                    bra = cen.ra.deg     #right asecention from the center coordinates of the bcgs being converted into degrees 
                    bdec = cen.dec.deg   #bcg declinations from the center coordinates of the bcgs being converted into degrees 
                    c2 = SkyCoord(bra, bdec, unit='deg', frame = 'fk5') #setting the coordinates to a variable  
                    
                    sep1 = c1.separation(c3)  #SZ and XRAY CENTER 
                    nsep1 = sep1.rad          #converting to radians 
                    asep1 = sep1.arcsec       #converting to arcseconds 
                    sin1 = math.sin(nsep1)    #taking the sine using the math module in order to start calculations 
                    distance1 = np.multiply(sin1,adj) #multiplying the value above by the defined parameter , "adj"
                    ndistance1 = abs(np.multiply(distance1,1000)) #taking the absolute value of the product the the previous value and 1000
            
            
                    sep2 = c2.separation(c3)  #BCG AND XRAY CENTER 
                    nsep2 = sep2.rad #converting the above calculation to radians
                    asep2 = sep2.arcsec #converting to arc seconds 
                    sin2 = math.sin(nsep2) #this is based off of the relationship between the angle and side length of the seperation 
                    distance2 = np.multiply(sin2,adj) 
                    ndistance2 = abs(np.multiply(distance2,1000)) #converting from meters to kilometers 
                    
                    #repeating the top steps for the BCG and SZ seperations
                    sep3 = c2.separation(c1)  #BCG AND SZ
                    nsep3 = sep3.rad
                    asep3 = sep3.arcsec
                    sin3 = math.sin(nsep3)
                    distance3 = np.multiply(sin3,adj)
                    ndistance3 = abs(np.multiply(distance3,1000))

                    #repeating steps for the BCG and Xray calculations 
                    sep4 = c2.separation(c4)  #BCG AND XRAY Peak
                    nsep4 = sep4.rad
                    asep4 = sep4.arcsec
                    sin4 = math.sin(nsep4)
                    distance4 = np.multiply(sin4,adj)
                    ndistance4 = abs(np.multiply(distance4,1000))

                    #the bdata is the data table that is created after calulations have been done. These are the values that will be used in the graphing code
                    #bdata.add_row((cnames1[i],z_cls[i],z_cl_err[i],Nmem[i],vd_biw[i],vd_err[i],szra[i],szdec[i],c4.ra.deg[0],c4.dec.deg[0],c3.ra.deg[0],c3.dec.deg[0],bra[0],bdec[0],ndistance3,ndistance4,ndistance2,ndistance1))
                    #this row reformats a table with the new data
                    #outstring =xray_name+'  '+str(z_cls[i])+'  '+str(count3)+'  '+str(ndistance1)+' '+str(ndistance2)+' '+str(ndistance3)+' '+str(ndistance4)+' \n'
                    #this outputs the actual table 
                    #ofil.write(outstring)
                    #outstring_peak  = ''+str(ndistance4)+'\n'
                    #outf.write(outstring_peak)
                    
                 
                   


#mdata.write(coutfil,format='ascii',overwrite=True) #this writes to the defined file 
#mdata.write(coutfits,format='fits',overwrite=True)
print(bdata)
#bdata.write(boutfil,format='ascii',overwrite=True)
#bdata.write(boutfits,format='fits',overwrite=True)

