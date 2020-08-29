import os
import numpy as np
import xml.etree.ElementTree as ET
from deg2str import *
from astropy import coordinates as coord
from astropy import units as u
from simbad import *
from datetime import *
from PIL import Image
from query_wsa_fits import *
from query_pso_fits import *
from jdcal import *
import pdb
import glob
import re
import astropy.coordinates as coord
import astropy.units as u

def psocolor(source_name,keepfiles=False,allcolor='#FFFF00',rejcolor='b',tm_color='r',plot=False,savepdf=True,secondary='',skipdownloads=False,circle_radius=0.0015,size=2.0,override_directory=None,primarypos_label=None,title=None,filename=None,buffer=False):
    # Set $FINDER_PATH in your bash_profile if you would like to control where the finder charts are output
    # size: arcmin
    # allwise: overplot AllWISE catalog positions
    # rejallwise = overplot AllWISE reject positions
    # tmass = overplot 2MASS psc positions
    # keepfiles = keep fits, tbl, and xml files
    # allcolor = color of allwise symbols
    # rejcolor = color of allwise reject symbols
    # tm_color = color of tmass symbols
    # plot = show plot (otherwise, finder is just made)
    # savepdf = save a pdf of the finder    
    
    circle_width = 2.5
    circle_alpha = 0.8
    fig_xsize = 11
    fig_ysize = 8.5
    
    #Convert sexagecimal
    if source_name.find(' ') == -1:
    	symbol = '+'
    	if source_name.find('+') == -1:
    		symbol = '-'
    	radec = re.split('[\+ | -]',source_name)
    	ras = radec[0][0:2]+' '+radec[0][2:4]+' '+radec[0][4:]
    	decs = radec[1][0:2]+' '+radec[1][2:4]+' '+radec[1][4:]
    	ra = coord.Angle(ras, unit=u.hour)
    	dec = coord.Angle(symbol+decs,unit=u.degree)
    	source_name = str(ra.degree)+' '+str(dec.degree)
    
    #Use buffer if needed
    if buffer:
        import matplotlib
        matplotlib.use('Agg')
        plot=False
    import pylab,pyfits,aplpy
    
    #Verify whether a working directory is set in the bash profile
    main_dir = None
    if override_directory:
        main_dir = override_directory
    else:
        proc = subprocess.Popen(["echo $FINDER_PATH"], stdout=subprocess.PIPE,shell=True)
        (out, err) = proc.communicate()
        if out != '\n':
            main_dir = out.split('\n')[0]
    if main_dir:
        initial_dir = os.getcwd()
        if not os.path.exists(main_dir):
            os.makedirs(main_dir)
        os.chdir(main_dir)
    
    #List of colors
    color_blue = '#377eb8'#RGB=[55,126,184]
    color_red = '#e41a1c'#RGB=[228,26,28]
    color_purple = '#b27bba'#RGB=[178,123,186]
    color_green = '#4daf4a'#RGB = [77,175,74]
    
    t1 = datetime.now()    
    
    ra,de = simbad(source_name)
    ra2 = None
    de2 = None
    if secondary:
        ra2,de2 = simbad(secondary)
    if filename is None:
        filename = source_name
    
    if skipdownloads is not True:
        print("Downloading Pan-Starrs data")
        #Remove previous data
        os.system("rm *_PSO_TMP.fits*")
        query_pso_fits(ra,de,size=size,output_file='PSO_TMP.fits')
    
    #If no PSO data could be downloaded, return
    if len(glob.glob('*_PSO_TMP.fits*')) == 0:
        print("No PSO data !")
        return
    
    if plot:
        pylab.ion()
    fig = pylab.figure(figsize=(fig_xsize,fig_ysize))
    pylab.rcParams['font.family'] = 'serif'
    
    #Create and plot RGB PSO image
    files = ['y_PSO_TMP.fits','i_PSO_TMP.fits','g_PSO_TMP.fits']
    aplpy.make_rgb_cube(files,'PSO_rgb.fits')
    impsoc = aplpy.FITSFigure('PSO_rgb_2d.fits',figure=fig)
    impsoc.add_label(0.05,0.9,'PSO $y$/$i$/$g$',relative=True,size='medium',color='k',bbox=dict(facecolor='white', alpha=0.5),horizontalalignment='left')

    meds = []
    mads = []
    devs = []
    for filei in files:
        datai = pyfits.getdata(filei)
        medi = np.nanmedian(datai)
        madi = np.nanmedian(abs(datai - np.nanmedian(datai)))
        devi = np.nanpercentile(datai,95) - np.nanpercentile(datai,5)
        meds.append(medi)
        mads.append(madi)
        devs.append(devi)
    
    mins = []
    maxs = []
    for i in range(0,len(files)):
        mini = (meds[i] - 2.0*mads[i])
        #maxi = (meds[i] + 10.0*mads[i])
        maxi = meds[i] + 2.0*devs[i]
        mins.append(mini)
        maxs.append(maxi)
    
    aplpy.make_rgb_image('PSO_rgb.fits','PSO_rgb.png',vmin_r=mins[0],vmin_g=mins[1],vmin_b=mins[2],vmax_r=maxs[0],vmax_g=maxs[1],vmax_b=maxs[2])
    
    impsoc.show_rgb('PSO_rgb.png')
    impsoc.hide_tick_labels()
    impsoc.ticks.hide()
    impsoc.hide_xaxis_label()
    impsoc.hide_yaxis_label()
    impsoc.recenter(ra,de,width=(size/60.0),height=(size/60.0))
    impsoc.show_circles(ra,de,edgecolor=color_red,linewidth=circle_width,facecolor='none',radius=circle_radius,alpha=circle_alpha)
    if secondary:
        impsoc.show_circles(ra2,de2,edgecolor=color_blue,linewidth=circle_width,facecolor='none',    radius=circle_radius,alpha=circle_alpha)
    
    # Remove files (or not)
    if keepfiles:
        pass
    else:
        print("Removing files...")
        cmdrm1 = "rm PSO_rgb.png PSO_rgb*.fits"
        os.system(cmdrm1)
    if savepdf:
        pylab.savefig(filename+'_pso_color.pdf')
    
    #Return to initial directory
    if main_dir:
        os.chdir(initial_dir)
    
    t2 = datetime.now()
    tdiff = (t2 - t1)
    print("PSO image creation took %s seconds" % (round(tdiff.total_seconds(),0)))