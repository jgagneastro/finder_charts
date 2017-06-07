import pylab,aplpy,pyfits,os
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

def finder(source_name,allwise=False,rejallwise=False,tmass=False,rejtmass=False,PSO=True,UKIDSS=True,VHS=True,keepfiles=False,allcolor='#FFFF00',rejcolor='b',tm_color='r',plot=False,savepdf=True,secondary='',addtext='',addtext2='',skipdownloads=False,circle_radius=0.0025,size=3.0,override_directory=None,primarypos_label=None,secondarypos_label=None):
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
    
    nxplot = 5
    nyplot = 3
    fig_xsize = 11
    fig_ysize = 8.5

    t1 = datetime.now()    

    ra,de = simbad(source_name)
    ra2 = None
    de2 = None
    if secondary:
        ra2,de2 = simbad(secondary)
    
    #Download xml file from IRSA
    xmlfile = "source.xml"
    if skipdownloads is not True:
        print "Getting xml file..."
        cmd = "wget -O "+xmlfile+" 'http://irsa.ipac.caltech.edu/applications/finderchart/servlet/api?locstr="+str(ra)+"+"+str(de)+"&subsetsize="+str(size)+"' "
        os.system(cmd)
    
    # parse xml file
    print "Parsing xml file..."
    
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    
    images = []
    for image in root.iter('image'):
        for child in image.getchildren():
            if child.tag == 'surveyname':
                surveyname = child.text
            if child.tag == 'band':
                band = child.text
            if child.tag == 'obsdate':
                obsdate = child.text
            if child.tag == 'fitsurl':
                fitsurl = child.text
        images.append([surveyname,band,obsdate,fitsurl])
    
    if skipdownloads is not True:
        print "Downloading DSS data..."
        for i in range(len(images)):
            if images[i][1] == 'DSS1 Blue':
                cmd1 = "wget -O DSS1_Blue.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'DSS1 Red':
                cmd1 = "wget -O DSS1_Red.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'DSS2 Blue':
                cmd1 = "wget -O DSS2_Blue.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'DSS2 Red':
                cmd1 = "wget -O DSS2_Red.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'DSS2 IR':
                cmd1 = "wget -O DSS2_IR.fits '"+images[i][3]+"'"
                os.system(cmd1)
        
        print "Downloading 2MASS data..."
        for i in range(len(images)):
            if images[i][1] == 'J':
                cmd1 = "wget -O 2MASS_J.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'H':
                cmd1 = "wget -O 2MASS_H.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'K':
                cmd1 = "wget -O 2MASS_K.fits '"+images[i][3]+"'"
                os.system(cmd1)
        
        print "Downloading WISE data"
        for i in range(len(images)):
            if images[i][1] == 'w1':
                cmd1 = "wget -O AllWISE_w1.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'w2':
                cmd1 = "wget -O AllWISE_w2.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'w3':
                cmd1 = "wget -O AllWISE_w3.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'w4':
                cmd1 = "wget -O AllWISE_w4.fits '"+images[i][3]+"'"
                os.system(cmd1)
        
        if UKIDSS:
            print "Downloading UKIDSS data"
            #Remove previous data
            os.system("rm *_UKIDSS_TMP.fits.gz")
            query_wsa_fits(ra,de,size=size,output_file='UKIDSS_TMP.fits.gz',filter='all',catalog='UKIDSS')
        
        if VHS:
            print "Downloading VHS data"
            #Remove previous data
            os.system("rm *_VHS_TMP.fits.gz")
            query_wsa_fits(ra,de,size=size,output_file='VHS_TMP.fits.gz',filter='all',catalog='VHS')
        
        if PSO:
            print "Downloading Pan-Starrs data"
            #Remove previous data
            os.system("rm *_PSO_TMP.fits*")
            query_pso_fits(ra,de,size=size,output_file='PSO_TMP.fits')
    
    #If no UKIDSS data could be downloaded, turn off the UKIDSS option
    if glob.glob('*_UKIDSS_TMP.fits*') is not True:
        UKIDSS = None
    
    #If no VHS data could be downloaded, turn off the VHS option
    if glob.glob('*_VHS_TMP.fits*') is not True:
        VHS = None
    
    #If no PSO data could be downloaded, turn off the PSO option
    pdb.set_trace()
    if glob.glob('*_PSO_TMP.fits*') is not True:
        PSO = None
    
    #Determine the amount of additional rows needed
    vertical_spacing = 0
    ukidss_spacing = 0
    vhs_spacing = 0
    tmass_spacing = 0
    allwise_spacing = 1
    if PSO:
        vertical_spacing += 1
        ukidss_spacing += 1
        vhs_spacing += 1
        tmass_spacing += 1
        allwise_spacing += 1
    if UKIDSS:
        vertical_spacing += 1
        vhs_spacing += 1
        tmass_spacing += 1
        allwise_spacing += 1
    if VHS:
        vertical_spacing += 1
        tmass_spacing += 1
        allwise_spacing += 1
    
    #Adapt window size
    nyplot += vertical_spacing
    fig_ysize *= np.sqrt(13.0/8.5)**(vertical_spacing-1.0)
    
    # Download the catalog data
    allwise_ra = None
    allwise_de = None
    if allwise:
        if skipdownloads is not True:
            cmd1 = 'wget -O allwise.tbl "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=Box&size=120.0&radunits=arcsec&objstr='+str(ra)+','+str(de)+'&catalog=allwise_p3as_psd&selcols=ra,dec&outfmt=1"'
            os.system(cmd1)
        try:
            awise = np.loadtxt('allwise.tbl',skiprows=27,unpack=True,usecols=(0,1))
            allwise_ra,allwise_de = awise
        except:
            print "No AllWISE sources found!"
            allwise=False
    
    rejallwise_ra = None
    rejallwise_de = None
    if rejallwise:
        if skipdownloads is not True:
            cmd2 = 'wget -O rejallwise.tbl "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=Box&size=120.0&radunits=arcsec&objstr='+str(ra)+','+str(de)+'&catalog=allwise_p3as_psr&selcols=ra,dec&outfmt=1"'
            os.system(cmd2)
        try:
            rejawise = np.loadtxt('rejallwise.tbl',skiprows=27,unpack=True,usecols=(0,1))
            rejallwise_ra,rejallwise_de = rejawise
        except:
            print "No AllWISE Reject sources found!"
            rejallwise=False
    
    tmass_ra = None
    tmass_de = None
    if tmass or rejtmass:
        if skipdownloads is not True:
            cmd3 = 'wget -O tmass.tbl "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=Box&size=120.0&radunits=arcsec&objstr='+str(ra)+','+str(de)+'&catalog=fp_psc&selcols=ra,dec&outfmt=1"'
            os.system(cmd3)
        try:
            tmass_psc = np.loadtxt('tmass.tbl',skiprows=37,unpack=True,usecols=(0,1))
            tmass_ra,tmass_de = tmass_psc
        except:
            print "No 2MASS sources found!"
            tmass=False
    
    rejtmass_ra = None
    rejtmass_de = None
    if rejtmass:
        if skipdownloads is not True:
            cmd4 = 'curl -o rejtmass.tbl "http://irsa.ipac.caltech.edu/TAP/sync?FORMAT=IPAC_TABLE&QUERY=SELECT+ra,dec,rel+FROM+pt_src_rej+WHERE+CONTAINS(POINT(\'J2000\',ra,dec),CIRCLE(\'J2000\','+str(ra)+','+str(de)+','+str(size/60.0)+'))=1"'
            #This is a temporary work-out, but only fetches rel='A' entries
            #cmd4 = 'wget -O rejtmass.tbl "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=Box&size='+str(size*60.0)+'&radunits=arcsec&objstr='+str(ra)+','+str(de)+'&catalog=pt_src_rej&selcols=ra,dec&outfmt=1"'
            os.system(cmd4)
        try:
            rejtmass_psc = np.loadtxt('rejtmass.tbl',skiprows=14,unpack=True,usecols=(0,1))
            rejtmass_ra,rejtmass_de = rejtmass_psc
        except:
            print "No 2MASS-Reject sources found!"
            rejtmass=False

    if plot:
        pylab.ion()
    fig = pylab.figure(figsize=(fig_xsize,fig_ysize))
    pylab.rcParams['font.family'] = 'serif'
    
    for i in range(len(images)):
        if images[i][1] == 'DSS1 Blue':
            min1, max1, im = oplotfits(fig,'DSS1_Blue.fits',nyplot,nxplot,1,ra,de,'DSS1 B',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)

        if images[i][1] == 'DSS1 Red':
            oplotfits(fig,'DSS1_Red.fits',nyplot,nxplot,2,ra,de,'DSS1 R',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
    
        if images[i][1] == 'DSS2 Blue':
            oplotfits(fig,'DSS2_Blue.fits',nyplot,nxplot,3,ra,de,'DSS2 B',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
    
        if images[i][1] == 'DSS2 Red':
            oplotfits(fig,'DSS2_Red.fits',nyplot,nxplot,4,ra,de,'DSS2 R',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        
        if images[i][1] == 'DSS2 IR':
            oplotfits(fig,'DSS2_IR.fits',nyplot,nxplot,5,ra,de,'DSS2 IR',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        
        if images[i][1] == 'J':
            oplotfits(fig,'2MASS_J.fits',nyplot,nxplot,6+tmass_spacing*nxplot,ra,de,'2MASS $J$',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,rejtmass_ra=rejtmass_ra,rejtmass_de=rejtmass_de,circle_radius=circle_radius)
        
        if images[i][1] == 'H':
            oplotfits(fig,'2MASS_H.fits',nyplot,nxplot,7+tmass_spacing*nxplot,ra,de,'2MASS $H$',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,rejtmass_ra=rejtmass_ra,rejtmass_de=rejtmass_de)
        
        if images[i][1] == 'K':
            oplotfits(fig,'2MASS_K.fits',nyplot,nxplot,8+tmass_spacing*nxplot,ra,de,'2MASS $K_S$',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,rejtmass_ra=rejtmass_ra,rejtmass_de=rejtmass_de)
        
        if images[i][1] == 'w1':
            wmin1, wmax1, im = oplotfits(fig,'AllWISE_w1.fits',nyplot,nxplot,6+allwise_spacing*nxplot,ra,de,'W1',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
            im.axis_labels.set_xtext('RA (arcmin)')
            im.axis_labels.set_ytext('Dec (arcmin)')
            
        if images[i][1] == 'w2':
            wmin2, wmax2, void = oplotfits(fig,'AllWISE_w2.fits',nyplot,nxplot,7+allwise_spacing*nxplot,ra,de,'W2',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        
        if images[i][1] == 'w3':
            wmin3, wmax3, void = oplotfits(fig,'AllWISE_w3.fits',nyplot,nxplot,8+allwise_spacing*nxplot,ra,de,'W3',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        
        if images[i][1] == 'w4':
            wmin4, wmax4, void = oplotfits(fig,'AllWISE_w4.fits',nyplot,nxplot,9+allwise_spacing*nxplot,ra,de,'W4',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
    
    # Plot PSO images
    if PSO:
        try:#g band
            #Get date and then make plot
            fitsfile = 'g_PSO_TMP.fits'
            hdulist = pyfits.open(fitsfile)
            jd = hdulist[0].header['MJD-OBS']
            hdulist.close()
            gregdate = jd2gcal(2400000.5, jd)
            year = '{0:.0f}'.format(gregdate[0]+gregdate[1]/12.0)
            wminpsog, wmaxpsog, void = oplotfits(fig,fitsfile,nyplot,nxplot,6,ra,de,'PSO $g$',year=year,ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
        try:#r band
            #Get date and then make plot
            fitsfile = 'r_PSO_TMP.fits'
            hdulist = pyfits.open(fitsfile)
            jd = hdulist[0].header['MJD-OBS']
            hdulist.close()
            gregdate = jd2gcal(2400000.5, jd)
            year = '{0:.0f}'.format(gregdate[0]+gregdate[1]/12.0)
            wminpsor, wmaxpsor, void = oplotfits(fig,fitsfile,nyplot,nxplot,7,ra,de,'PSO $r$',year=year,ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
        try:#i band
            #Get date and then make plot
            fitsfile = 'i_PSO_TMP.fits'
            hdulist = pyfits.open(fitsfile)
            jd = hdulist[0].header['MJD-OBS']
            hdulist.close()
            gregdate = jd2gcal(2400000.5, jd)
            year = '{0:.0f}'.format(gregdate[0]+gregdate[1]/12.0)
            wminpsoi, wmaxpsoi, void = oplotfits(fig,fitsfile,nyplot,nxplot,8,ra,de,'PSO $i$',year=year,ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
        try:#z band
            #Get date and then make plot
            fitsfile = 'z_PSO_TMP.fits'
            hdulist = pyfits.open(fitsfile)
            jd = hdulist[0].header['MJD-OBS']
            hdulist.close()
            gregdate = jd2gcal(2400000.5, jd)
            year = '{0:.0f}'.format(gregdate[0]+gregdate[1]/12.0)
            wminpsoz, wmaxpsoz, void = oplotfits(fig,fitsfile,nyplot,nxplot,9,ra,de,'PSO $z$',year=year,ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
        try:#y band
            #Get date and then make plot
            fitsfile = 'y_PSO_TMP.fits'
            hdulist = pyfits.open(fitsfile)
            jd = hdulist[0].header['MJD-OBS']
            hdulist.close()
            gregdate = jd2gcal(2400000.5, jd)
            year = '{0:.0f}'.format(gregdate[0]+gregdate[1]/12.0)
            wminpsoy, wmaxpsoy, void = oplotfits(fig,fitsfile,nyplot,nxplot,10,ra,de,'PSO $y$',year=year,ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
    
    # Plot UKIDSS images
    if UKIDSS:
        try:#Y band
            #Get date and then make plot
            fitsfile = 'Y_UKIDSS_TMP.fits.gz'
            hdulist = pyfits.open(fitsfile)
            year = hdulist[0].header['UTDATE'][0:4]
            hdulist.close()
            oplotfits(fig,fitsfile,nyplot,nxplot,6+ukidss_spacing*nxplot,ra,de,'UKIDSS Y',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
        try:#J band
            #Get date and then make plot
            fitsfile = 'J_UKIDSS_TMP.fits.gz'
            hdulist = pyfits.open(fitsfile)
            year = hdulist[0].header['UTDATE'][0:4]
            hdulist.close()
            oplotfits(fig,fitsfile,nyplot,nxplot,7+ukidss_spacing*nxplot,ra,de,'UKIDSS J',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
        try:#H band
            #Get date and then make plot
            fitsfile = 'H_UKIDSS_TMP.fits.gz'
            hdulist = pyfits.open(fitsfile)
            year = hdulist[0].header['UTDATE'][0:4]
            hdulist.close()
            oplotfits(fig,fitsfile,nyplot,nxplot,8+ukidss_spacing*nxplot,ra,de,'UKIDSS H',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
        try:#K band
            #Get date and then make plot
            fitsfile = 'K_UKIDSS_TMP.fits.gz'
            hdulist = pyfits.open(fitsfile)
            year = hdulist[0].header['UTDATE'][0:4]
            hdulist.close()
            oplotfits(fig,fitsfile,nyplot,nxplot,9+ukidss_spacing*nxplot,ra,de,'UKIDSS K',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
    
    # Plot VHS images
    if VHS:
        try:#Y band
            #Get date and then make plot
            fitsfile = 'Y_VHS_TMP.fits.gz'
            hdulist = pyfits.open(fitsfile)
            year = hdulist[1].header['DATE'][0:4]
            hdulist.close()
            oplotfits(fig,fitsfile,nyplot,nxplot,6+vhs_spacing*nxplot,ra,de,'VHS $Y$',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
        try:#J band
            #Get date and then make plot
            fitsfile = 'J_VHS_TMP.fits.gz'
            hdulist = pyfits.open(fitsfile)
            year = hdulist[1].header['DATE'][0:4]
            hdulist.close()
            oplotfits(fig,fitsfile,nyplot,nxplot,7+vhs_spacing*nxplot,ra,de,'VHS $J$',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
        try:#H band
            #Get date and then make plot
            fitsfile = 'H_VHS_TMP.fits.gz'
            hdulist = pyfits.open(fitsfile)
            year = hdulist[1].header['DATE'][0:4]
            hdulist.close()
            oplotfits(fig,fitsfile,nyplot,nxplot,8+vhs_spacing*nxplot,ra,de,'VHS $H$',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
        try:#Ks band
            #Get date and then make plot
            fitsfile = 'Ks_VHS_TMP.fits.gz'
            hdulist = pyfits.open(fitsfile)
            year = hdulist[1].header['DATE'][0:4]
            hdulist.close()
            oplotfits(fig,fitsfile,nyplot,nxplot,9+vhs_spacing*nxplot,ra,de,'VHS $K_S$',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius)
        except:
            pass
    
    # Create and plot RGB AllWISE image
    files = ['AllWISE_w3.fits','AllWISE_w2.fits','AllWISE_w1.fits']
    aplpy.make_rgb_cube(files,'AllWISE_rgb.fits')
    im20 = aplpy.FITSFigure('AllWISE_rgb_2d.fits',figure=fig,subplot=(nyplot,nxplot,10+allwise_spacing*nxplot))
    meds = []
    mads = []
    devs = []
    for filei in files:
        datai = pyfits.getdata(filei)
        medi = np.nanmedian(datai)
        madi = np.nanmedian(abs(datai - np.nanmedian(datai)))
        devi = np.percentile(datai,95) - np.percentile(datai,5)
        meds.append(medi)
        mads.append(madi)
        devs.append(devi)
    
    mins = []
    maxs = []
    for i in range(0,len(files)):
        mini = meds[i] - 2.0*mads[i]
        #maxi = meds[i] + 10.0*mads[i]
        maxi = meds[i] + 2.0*devs[i]
        mins.append(mini)
        maxs.append(maxi)
    
    #aplpy.make_rgb_image('AllWISE_rgb.fits','AllWISE_rgb.png')
    #aplpy.make_rgb_image('AllWISE_rgb.fits','AllWISE_rgb.png',vmin_r=mins[0],vmin_g=mins[1],vmin_b=mins[2])
    
    aplpy.make_rgb_image('AllWISE_rgb.fits','AllWISE_rgb.png',vmin_r=mins[0],vmin_g=mins[1],vmin_b=mins[2],vmax_r=maxs[0],vmax_g=maxs[1],vmax_b=maxs[2])
    
    im20.show_rgb('AllWISE_rgb.png')
    im20.hide_tick_labels()
    im20.ticks.set_color('k')
    im20.ticks.set_minor_frequency(0)
    im20.ticks.set_xspacing(0.5/60.0)
    im20.ticks.set_yspacing(0.5/60.0)
    im20.hide_xaxis_label()
    im20.hide_yaxis_label()
    im20.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
    im20.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=circle_radius)
    if secondary:
        im20.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=circle_radius)
    
    #Create and plot RGB PSO image only if there's enough space (vertical_spacing >= 2)
    if PSO and vertical_spacing >= 2:
        try:
            files = ['y_PSO_TMP.fits','i_PSO_TMP.fits','g_PSO_TMP.fits']
            aplpy.make_rgb_cube(files,'PSO_rgb.fits')
            impsoc = aplpy.FITSFigure('PSO_rgb_2d.fits',figure=fig,subplot=(nyplot,nxplot,15))
            
            meds = []
            mads = []
            devs = []
            #exptimes = []
            for filei in files:
                datai = pyfits.getdata(filei)
                medi = np.nanmedian(datai)
                madi = np.nanmedian(abs(datai - np.nanmedian(datai)))
                devi = np.percentile(datai,95) - np.percentile(datai,5)
                #hdulist = pyfits.open(filei)
                #exptimei = hdulist[0].header['EXPTIME']
                #hdulist.close()
                meds.append(medi)
                mads.append(madi)
                devs.append(devi)
                #exptimes.append(exptimei)
            
            mins = []
            maxs = []
            for i in range(0,len(files)):
                mini = (meds[i] - 2.0*mads[i])#/exptimes[i]*np.average(exptimes)
                maxi = (meds[i] + 10.0*mads[i])#/exptimes[i]*np.average(exptimes)
                #maxi = meds[i] + 5.0*devs[i]
                mins.append(mini)
                maxs.append(maxi)
            
            #aplpy.make_rgb_image('PSO_rgb.fits','PSO_rgb.png',vmin_r=mins[0],vmin_g=mins[1],vmin_b=mins[2],vmax_r=maxs[0],vmax_g=maxs[1],vmax_b=maxs[2])
            
            aplpy.make_rgb_image('PSO_rgb.fits','PSO_rgb.png')
            
            #aplpy.make_rgb_image('PSO_rgb.fits','PSO_rgb.png',vmin_r=wminpsoy,vmin_g=wminpsoi,vmin_b=wminpsog,vmax_r=wmaxpsoy,vmax_g=wmaxpsoi,vmax_b=wmaxpsog)
            
            impsoc.show_rgb('PSO_rgb.png')
            impsoc.hide_tick_labels()
            impsoc.ticks.set_color('k')
            impsoc.ticks.set_minor_frequency(0)
            impsoc.ticks.set_xspacing(0.5/60.0)
            impsoc.ticks.set_yspacing(0.5/60.0)
            impsoc.hide_xaxis_label()
            impsoc.hide_yaxis_label()
            impsoc.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            impsoc.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=circle_radius)
            if secondary:
                impsoc.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',    radius=circle_radius)
        except:
            pass
    
    bottom_space = 0.05
    top_space = 0.05
    pylab.subplots_adjust(left=0.05,right=0.95,bottom=bottom_space,top=1.0-top_space,wspace=0.05,hspace=0.05)
    # Add Labels
    ras = deg2str(ra)
    des = deg2str(de,dec=1)
    c1 = coord.ICRS(ra,de,unit=(u.degree,u.degree))
    c2 = c1.galactic    
    
    sptext = fig.add_subplot(nyplot,nxplot,9+tmass_spacing*nxplot)
    sptext.axis('off')
    xlabels = .02
    ytoplabels = .89
    ydeltalabels = .12

    sptext.annotate(r'$\alpha$ = '+ras+'\t('+str(round(ra,6))+')',xy=(xlabels,ytoplabels),fontsize=15)
    sptext.annotate(r'$\delta$ = '+des+'\t('+str(round(de,6))+')',xy=(xlabels,ytoplabels-ydeltalabels),fontsize=15)
    sptext.annotate(r'$l$ = '+str(round(c2.l.degree,3)),xy=(xlabels,ytoplabels-ydeltalabels*2),fontsize=15)
    sptext.annotate(r'$b$ = '+str(round(c2.b.degree,3)),xy=(xlabels,ytoplabels-ydeltalabels*3),fontsize=15)
    if addtext:
        for i, addtexti in enumerate(addtext):
            sptext.annotate(addtexti,xy=(xlabels,ytoplabels-ydeltalabels*(4+i)),fontsize=15)
    if addtext2:
        sptext.annotate(addtext2,xy=(xlabels,ytoplabels-ydeltalabels*5),fontsize=14)
    if allwise:
        sptext.annotate('AllWISE catalog sources',xy=(xlabels+1e-3,ytoplabels-ydeltalabels*6-1e-3),fontsize=15,color='k')
        sptext.annotate('AllWISE catalog sources',xy=(xlabels,ytoplabels-ydeltalabels*6),fontsize=15,color=allcolor)
    if rejallwise:
        sptext.annotate('AllWISE reject sources',xy=(xlabels+1e-3,ytoplabels-ydeltalabels*7-1e-3),fontsize=15,color='k')
        sptext.annotate('AllWISE reject sources',xy=(xlabels,ytoplabels-ydeltalabels*7),fontsize=15,color=rejcolor)
    if tmass:
        sptext.annotate('2MASS catalog sources',xy=(xlabels+1e-3,ytoplabels-ydeltalabels*8-1e-3),fontsize=15,color='k')
        sptext.annotate('2MASS catalog sources',xy=(xlabels,ytoplabels-ydeltalabels*8),fontsize=15,color=tm_color)
    if rejtmass:
        sptext.annotate('2MASS reject sources',xy=(xlabels+1e-3,ytoplabels-ydeltalabels*9-1e-3),fontsize=15,color='k')
        sptext.annotate('2MASS reject sources',xy=(xlabels,ytoplabels-ydeltalabels*9),fontsize=15,color=tm_color)
    
    sptext2 = fig.add_subplot(nyplot,nxplot,10+tmass_spacing*nxplot)
    sptext2.axis('off')
    if primarypos_label:
        sptext2.annotate(primarypos_label,xy=(0.99,ytoplabels),fontsize=15,color='r',horizontalalignment='right')
    if secondary and secondarypos_label:
        sptext2.annotate(secondarypos_label,xy=(0.99,ytoplabels-ydeltalabels),fontsize=15,color='b',horizontalalignment='right')

    # Remove files (or not)
    if keepfiles:
        pass
    else:
        print "Removing files..."
        cmdrm1 = "rm source.xml 2MASS*.fits AllWISE*.fits DSS*.fits AllWISE_rgb.png *UKIDSS_TMP.fits.gz *VHS_TMP.fits.gz *PSO_TMP.fits* UKIDSS_rgb*.fits UKIDSS_rgb.png PSO_rgb.png PSO_rgb*.fits"
        os.system(cmdrm1)
        if allwise:
            cmdrm2 = "rm allwise.tbl"
            os.system(cmdrm2)
            if rejallwise:
                cmdrm3 = "rm rejallwise.tbl"
                os.system(cmdrm3)
        if tmass:
            cmdrm4 = "rm tmass.tbl"
            os.system(cmdrm4)
    if savepdf:
        pylab.savefig(source_name+'.pdf')
    
    #Return to initial directory
    if main_dir:
        os.chdir(initial_dir)
    
    t2 = datetime.now()
    tdiff = (t2 - t1)
    print "Finder creation took %s seconds" % (round(tdiff.total_seconds(),0))

def oplotfits(fig,fitsfile,nyplot,nxplot,position,ra,de,label,year='',xlabel=0.05,ra2=None,de2=None,north=False,hdu=0,allwise=False,rejallwise=False,tmass=False,allcolor='#FFFF00',rejcolor='b',tm_color='r',secondary='',allwise_ra=None,allwise_de=None,rejallwise_ra=None,rejallwise_de=None,tmass_ra=None,tmass_de=None,rejtmass_ra=None,rejtmass_de=None,circle_radius=0.0025):
    
    im = aplpy.FITSFigure(fitsfile,figure=fig,subplot=(nyplot,nxplot,position),north=north,hdu=hdu)
    im.hide_tick_labels()
    im.ticks.set_color('k')
    im.ticks.set_minor_frequency(0)
    im.ticks.set_xspacing(0.5/60.0)
    im.ticks.set_yspacing(0.5/60.0)
    im.hide_xaxis_label()
    im.hide_yaxis_label()
    data1 = pyfits.getdata(fitsfile)
    med1 = np.nanmedian(data1)
    mad1 = np.nanmedian(abs(data1 - np.nanmedian(data1)))
    min1 = med1 - 2.0*mad1
    max1 = med1 + 10.0*mad1
    im.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
    im.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
    im.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=circle_radius)
    if secondary:
        im.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=circle_radius)
    im.add_label(xlabel,0.9,label,relative=True,size='medium',color='k',bbox=dict(facecolor='white', alpha=0.5),horizontalalignment='left')
    if year:
        im.add_label(xlabel,0.1,year,relative=True,size='medium',color='k',bbox=dict(facecolor='white', alpha=0.5),horizontalalignment='left')
    if allwise:
        im.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
    if rejallwise:
        im.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
    if tmass:
        im.show_circles(tmass_ra,tmass_de,edgecolor=allcolor,radius=0.0015,linewidth=1.5)
    if (rejtmass_ra is not None) and (rejtmass_de is not None):
        im.show_circles(rejtmass_ra,rejtmass_de,edgecolor=allcolor,radius=0.001,linewidth=0.5,alpha=0.8)
        #If possible, also display 2MASS main catalog sources if rejtmass is set
        if (not tmass) and (tmass_ra is not None) and (tmass_de is not None):
            im.show_circles(tmass_ra,tmass_de,edgecolor=allcolor,radius=0.0015,linewidth=1.5,alpha=0.8)
    return min1, max1, im