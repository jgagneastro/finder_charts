import numpy,pylab,aplpy,pyfits,os
import xml.etree.ElementTree as ET
from deg2str import *
from astropy import coordinates as coord
from astropy import units as u
from simbad import *
from datetime import *
from PIL import Image
from query_wsa_fits import *
import pdb

def finderold(source_name,allwise=False,rejallwise=False,tmass=False,SDSS=False,UKIDSS=False,VHS=False,keepfiles=False,allcolor='#FFFF00',rejcolor='b',tm_color='r',plot=True,savepdf=False,secondary='',addtext='',addtext2=''):
    # allwise: overplot AllWISE catalog positions
    # rejallwise = overplot AllWISE reject positions
    # tmass = overplot 2MASS psc positions
    # SDSS = download and plot SDSS data (adds a significant amount of time to run)
    # keepfiles = keep fits, tbl, and xml files
    # allcolor = color of allwise symbols
    # rejcolor = color of allwise reject symbols
    # tm_color = color of tmass symbols
    # plot = show plot (otherwise, finder is just made)
    # savepdf = save a pdf of the finder    
    #pdb.set_trace()

    os.chdir('/Users/gagne/tmp')

    nxplot = 5
    nyplot = 6
    fig_xsize = 11
    fig_ysize = 13.0
    fig_ysize0 = 8.5
#    if nyplot == 5:
#        fig_ysize = 11.0
#    if nyplot == 6:
#        fig_ysize = 13.0

    t1 = datetime.now()    

    ra,de = simbad(source_name)
    if secondary:
        ra2,de2 = simbad(secondary)

    xmlfile = "source.xml"
    
    cmd = "wget -O "+xmlfile+" 'http://irsa.ipac.caltech.edu/applications/finderchart/servlet/api?locstr="+str(ra)+"+"+str(de)+"&subsetsize=3.0' "

    os.system(cmd)
    print source_name
    print ra
    print de
    return

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
        query_wsa_fits(ra,de,3.0,output_file='UKIDSS_TMP.fits.gz',filter='all',catalog='UKIDSS')
    
    if VHS:
        print "Downloading VHS data"
        query_wsa_fits(ra,de,3.0,output_file='VHS_TMP.fits.gz',filter='all',catalog='VHS')
    
    if SDSS:
        print "Downloading SDSS data"
        for i in range(len(images)):
            if images[i][1] == 'u':
                cmd1 = "wget -O SDSS_u.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'g':
                cmd1 = "wget -O SDSS_g.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'r':
                cmd1 = "wget -O SDSS_r.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'i':
                cmd1 = "wget -O SDSS_i.fits '"+images[i][3]+"'"
                os.system(cmd1)
            if images[i][1] == 'z':
                cmd1 = "wget -O SDSS_z.fits '"+images[i][3]+"'"
                os.system(cmd1)
    
    ##Adapt the window size according to available data
    #nxplot = 5
    #nyplot = 4
    #if UKIDSS or VHS:
    #    nyplot = 5
    #if UKIDSS and VHS:
    #    nyplot = 6
    #    pad_vhs_positions = nxplot
    #
    #fig_xsize = 11
    #fig_ysize = 8.5
    #if nyplot == 5:
    #    fig_ysize = 11.0
    #if nyplot == 6:
    #    fig_ysize = 13.0
    
    # Plot the data
    
    if allwise:
        cmd1 = 'wget -O allwise.tbl "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=Box&size=120.0&radunits=arcsec&objstr='+str(ra)+','+str(de)+'&catalog=allwise_p3as_psd&selcols=ra,dec&outfmt=1"'
        os.system(cmd1)
        try:
            awise = numpy.loadtxt('allwise.tbl',skiprows=27,unpack=True,usecols=(0,1))
            allwise_ra,allwise_de = awise
        except:
            print "No AllWISE sources found!"
            allwise=False
    
    if rejallwise:
        cmd2 = 'wget -O rejallwise.tbl "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=Box&size=120.0&radunits=arcsec&objstr='+str(ra)+','+str(de)+'&catalog=allwise_p3as_psr&selcols=ra,dec&outfmt=1"'
        os.system(cmd2)
        try:
            rejawise = numpy.loadtxt('rejallwise.tbl',skiprows=27,unpack=True,usecols=(0,1))
            rejallwise_ra,rejallwise_de = rejawise
        except:
            print "No AllWISE Reject sources found!"
            rejallwise=False

        if tmass:
            cmd3 = 'wget -O tmass.tbl "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=Box&size=120.0&radunits=arcsec&objstr='+str(ra)+','+str(de)+'&catalog=fp_psc&selcols=ra,dec&outfmt=1"'
            os.system(cmd3)
        try:
            tmass_psc = numpy.loadtxt('tmass.tbl',skiprows=37,unpack=True,usecols=(0,1))
            tmass_ra,tmass_de = tmass_psc
        except:
            print "No 2MASS sources found!"
            tmass=False
    
    if plot:
        pylab.ion()
    fig = pylab.figure(figsize=(fig_xsize,fig_ysize))
    pylab.rcParams['font.family'] = 'serif'
    
    for i in range(len(images)):
        if images[i][1] == 'DSS1 Blue':
            im01 = aplpy.FITSFigure('DSS1_Blue.fits',figure=fig,subplot=(nyplot,nxplot,1))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.hide_xaxis_label()
            im01.hide_yaxis_label()
            data1 = pyfits.getdata('DSS1_Blue.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.2,0.9,'DSS1 B',relative=True,size='medium',color='k')
            im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
            
        if images[i][1] == 'DSS1 Red':
            im01 = aplpy.FITSFigure('DSS1_Red.fits',figure=fig,subplot=(nyplot,nxplot,2))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.hide_xaxis_label()
            im01.hide_yaxis_label()
            data1 = pyfits.getdata('DSS1_Red.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.2,0.9,'DSS1 R',relative=True,size='medium',color='k')
            im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
    
        if images[i][1] == 'DSS2 Blue':
            im01 = aplpy.FITSFigure('DSS2_Blue.fits',figure=fig,subplot=(nyplot,nxplot,3))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.hide_xaxis_label()
            im01.hide_yaxis_label()
            data1 = pyfits.getdata('DSS2_Blue.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.2,0.9,'DSS2 B',relative=True,size='medium',color='k')
            im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
    
        if images[i][1] == 'DSS2 Red':
            im01 = aplpy.FITSFigure('DSS2_Red.fits',figure=fig,subplot=(nyplot,nxplot,4))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.hide_xaxis_label()
            im01.hide_yaxis_label()
            data1 = pyfits.getdata('DSS2_Red.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.2,0.9,'DSS2 R',relative=True,size='medium',color='k')
            im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
            
        if images[i][1] == 'DSS2 IR':
            im01 = aplpy.FITSFigure('DSS2_IR.fits',figure=fig,subplot=(nyplot,nxplot,5))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.hide_xaxis_label()
            im01.hide_yaxis_label()
            data1 = pyfits.getdata('DSS2_IR.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.2,0.9,'DSS2 IR',relative=True,size='medium',color='k')
            im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
    
        if SDSS:
            if images[i][1] == 'u':
                im01 = aplpy.FITSFigure('SDSS_u.fits',figure=fig,subplot=(nyplot,nxplot,6))
                im01.hide_tick_labels()
                im01.ticks.set_color('k')
                im01.ticks.set_minor_frequency(0)
                im01.ticks.set_xspacing(0.5/60.0)
                im01.ticks.set_yspacing(0.5/60.0)
                im01.hide_xaxis_label()
                im01.hide_yaxis_label()
                data1 = pyfits.getdata('SDSS_u.fits')
                med1 = numpy.nanmedian(data1)
                mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
                min1 = med1 - 2.0*mad1
                max1 = med1 + 10.0*mad1
                im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
                im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
                im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
                if secondary:
                    im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
                im01.add_label(0.1,0.9,'$u$',relative=True,size='medium',color='k')
                im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
                if allwise:
                    im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
                if rejallwise:
                    im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
                if tmass:
                    im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
    
            if images[i][1] == 'g':
                im01 = aplpy.FITSFigure('SDSS_g.fits',figure=fig,subplot=(nyplot,nxplot,7))
                im01.hide_tick_labels()
                im01.ticks.set_color('k')
                im01.ticks.set_minor_frequency(0)
                im01.ticks.set_xspacing(0.5/60.0)
                im01.ticks.set_yspacing(0.5/60.0)
                im01.hide_xaxis_label()
                im01.hide_yaxis_label()
                data1 = pyfits.getdata('SDSS_g.fits')
                med1 = numpy.nanmedian(data1)
                mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
                min1 = med1 - 2.0*mad1
                max1 = med1 + 10.0*mad1
                im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
                im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
                im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
                if secondary:
                    im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
                im01.add_label(0.1,0.9,'$g$',relative=True,size='medium',color='k')
                im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
                if allwise:
                    im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
                if rejallwise:
                    im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
                if tmass:
                    im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
    
            if images[i][1] == 'r':
                im01 = aplpy.FITSFigure('SDSS_r.fits',figure=fig,subplot=(nyplot,nxplot,8))
                im01.hide_tick_labels()
                im01.ticks.set_color('k')
                im01.ticks.set_minor_frequency(0)
                im01.ticks.set_xspacing(0.5/60.0)
                im01.ticks.set_yspacing(0.5/60.0)
                im01.hide_xaxis_label()
                im01.hide_yaxis_label()
                data1 = pyfits.getdata('SDSS_r.fits')
                med1 = numpy.nanmedian(data1)
                mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
                min1 = med1 - 2.0*mad1
                max1 = med1 + 10.0*mad1
                im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
                im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
                im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
                if secondary:
                    im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
                im01.add_label(0.1,0.9,'$r$',relative=True,size='medium',color='k')
                im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
                if allwise:
                    im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
                if rejallwise:
                    im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
                if tmass:
                    im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
    
            if images[i][1] == 'i':
                im01 = aplpy.FITSFigure('SDSS_i.fits',figure=fig,subplot=(nyplot,nxplot,9))
                im01.hide_tick_labels()
                im01.ticks.set_color('k')
                im01.ticks.set_minor_frequency(0)
                im01.ticks.set_xspacing(0.5/60.0)
                im01.ticks.set_yspacing(0.5/60.0)
                im01.hide_xaxis_label()
                im01.hide_yaxis_label()
                data1 = pyfits.getdata('SDSS_i.fits')
                med1 = numpy.nanmedian(data1)
                mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
                min1 = med1 - 2.0*mad1
                max1 = med1 + 10.0*mad1
                im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
                im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
                im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
                if secondary:
                    im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
                im01.add_label(0.1,0.9,'$i$',relative=True,size='medium',color='k')
                im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
                if allwise:
                    im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
                if rejallwise:
                    im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
                if tmass:
                    im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
    
            if images[i][1] == 'z':
                im01 = aplpy.FITSFigure('SDSS_z.fits',figure=fig,subplot=(nyplot,nxplot,10))
                im01.hide_tick_labels()
                im01.ticks.set_color('k')
                im01.ticks.set_minor_frequency(0)
                im01.ticks.set_xspacing(0.5/60.0)
                im01.ticks.set_yspacing(0.5/60.0)
                im01.hide_xaxis_label()
                im01.hide_yaxis_label()
                data1 = pyfits.getdata('SDSS_z.fits')
                med1 = numpy.nanmedian(data1)
                mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
                min1 = med1 - 2.0*mad1
                max1 = med1 + 10.0*mad1
                im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
                im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
                im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
                if secondary:
                    im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
                im01.add_label(0.1,0.9,'$z$',relative=True,size='medium',color='k')
                im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
                if allwise:
                    im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
                if rejallwise:
                    im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
                if tmass:
                    im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
    
        if images[i][1] == 'J':
            im01 = aplpy.FITSFigure('2MASS_J.fits',figure=fig,subplot=(nyplot,nxplot,21))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.hide_xaxis_label()
            im01.hide_yaxis_label()
            data1 = pyfits.getdata('2MASS_J.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.25,0.9,'2MASS J',relative=True,size='medium',color='k')
            im01.add_label(0.1,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
    
        if images[i][1] == 'H':
            im01 = aplpy.FITSFigure('2MASS_H.fits',figure=fig,subplot=(nyplot,nxplot,22))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.hide_xaxis_label()
            im01.hide_yaxis_label()
            data1 = pyfits.getdata('2MASS_H.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.25,0.9,'2MASS H',relative=True,size='medium',color='k')
            im01.add_label(0.1,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        
        if images[i][1] == 'K':
            im01 = aplpy.FITSFigure('2MASS_K.fits',figure=fig,subplot=(nyplot,nxplot,23))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.hide_xaxis_label()
            im01.hide_yaxis_label()
            data1 = pyfits.getdata('2MASS_K.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.25,0.9,'2MASS K',relative=True,size='medium',color='k')
            im01.add_label(0.1,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        
        if images[i][1] == 'w1':
            im01 = aplpy.FITSFigure('AllWISE_w1.fits',figure=fig,subplot=(nyplot,nxplot,26))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.axis_labels.set_xtext('RA (arcmin)')
            im01.axis_labels.set_ytext('Dec (arcmin)')
            data1 = pyfits.getdata('AllWISE_w1.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.10,0.9,'W1',relative=True,size='medium',color='k')
            im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            #im01.add_scalebar((0.5/60.0),color='k')
            wmin1 = min1
            wmax1 = max1
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        
        if images[i][1] == 'w2':
            im01 = aplpy.FITSFigure('AllWISE_w2.fits',figure=fig,subplot=(nyplot,nxplot,27))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.hide_xaxis_label()
            im01.hide_yaxis_label()
            data1 = pyfits.getdata('AllWISE_w2.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.10,0.9,'W2',relative=True,size='medium',color='k')
            im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            wmin2 = min1
            wmax2 = max1
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        
        if images[i][1] == 'w3':
            im01 = aplpy.FITSFigure('AllWISE_w3.fits',figure=fig,subplot=(nyplot,nxplot,28))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.hide_xaxis_label()
            im01.hide_yaxis_label()
            data1 = pyfits.getdata('AllWISE_w3.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.10,0.9,'W3',relative=True,size='medium',color='k')
            im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            wmin3 = min1
            wmax3 = max1
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        
        if images[i][1] == 'w4':
            im01 = aplpy.FITSFigure('AllWISE_w4.fits',figure=fig,subplot=(nyplot,nxplot,29))
            im01.hide_tick_labels()
            im01.ticks.set_color('k')
            im01.ticks.set_minor_frequency(0)
            im01.ticks.set_xspacing(0.5/60.0)
            im01.ticks.set_yspacing(0.5/60.0)
            im01.hide_xaxis_label()
            im01.hide_yaxis_label()
            data1 = pyfits.getdata('AllWISE_w4.fits')
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            im01.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            im01.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            im01.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                im01.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            im01.add_label(0.10,0.9,'W4',relative=True,size='medium',color='k')
            im01.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                im01.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                im01.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                im01.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
    
# Plot UKIDSS images
    if UKIDSS:
        try:#Y band
            fitsfile = 'Y_UKIDSS_TMP.fits.gz'
            imuk = aplpy.FITSFigure(fitsfile,figure=fig,subplot=(nyplot,nxplot,11))
            imuk.hide_tick_labels()
            imuk.ticks.set_color('k')
            imuk.ticks.set_minor_frequency(0)
            imuk.ticks.set_xspacing(0.5/60.0)
            imuk.ticks.set_yspacing(0.5/60.0)
            imuk.hide_xaxis_label()
            imuk.hide_yaxis_label()
            data1 = pyfits.getdata(fitsfile)
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            imuk.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            imuk.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            imuk.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                imuk.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            imuk.add_label(0.3,0.9,'UKIDSS Y',relative=True,size='medium',color='k')
            #Need to get year
            #imuk.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                imuk.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                imuk.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                imuk.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        except:
            pass
        try:#J band
            fitsfile = 'J_UKIDSS_TMP.fits.gz'
            imuk = aplpy.FITSFigure(fitsfile,figure=fig,subplot=(nyplot,nxplot,12),north=True,hdu=1)
            imuk.hide_tick_labels()
            imuk.ticks.set_color('k')
            imuk.ticks.set_minor_frequency(0)
            imuk.ticks.set_xspacing(0.5/60.0)
            imuk.ticks.set_yspacing(0.5/60.0)
            imuk.hide_xaxis_label()
            imuk.hide_yaxis_label()
            data1 = pyfits.getdata(fitsfile)
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            imuk.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            imuk.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            imuk.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                imuk.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            imuk.add_label(0.3,0.9,'UKIDSS J',relative=True,size='medium',color='k')
            #Need to get year
            #imuk.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                imuk.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                imuk.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                imuk.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        except:
            pass
        try:#H band
            fitsfile = 'H_UKIDSS_TMP.fits.gz'
            imuk = aplpy.FITSFigure(fitsfile,figure=fig,subplot=(nyplot,nxplot,13),north=True,hdu=1)
            imuk.hide_tick_labels()
            imuk.ticks.set_color('k')
            imuk.ticks.set_minor_frequency(0)
            imuk.ticks.set_xspacing(0.5/60.0)
            imuk.ticks.set_yspacing(0.5/60.0)
            imuk.hide_xaxis_label()
            imuk.hide_yaxis_label()
            data1 = pyfits.getdata(fitsfile)
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            imuk.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            imuk.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            imuk.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                imuk.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            imuk.add_label(0.3,0.9,'UKIDSS H',relative=True,size='medium',color='k')
            #Need to get year
            #imuk.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                imuk.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                imuk.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                imuk.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        except:
            pass
        try:#K band
            fitsfile = 'K_UKIDSS_TMP.fits.gz'
            imuk = aplpy.FITSFigure(fitsfile,figure=fig,subplot=(nyplot,nxplot,14),north=True,hdu=1)
            imuk.hide_tick_labels()
            imuk.ticks.set_color('k')
            imuk.ticks.set_minor_frequency(0)
            imuk.ticks.set_xspacing(0.5/60.0)
            imuk.ticks.set_yspacing(0.5/60.0)
            imuk.hide_xaxis_label()
            imuk.hide_yaxis_label()
            data1 = pyfits.getdata(fitsfile)
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            imuk.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            imuk.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            imuk.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                imuk.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            imuk.add_label(0.3,0.9,'UKIDSS K',relative=True,size='medium',color='k')
            #Need to get year
            #imuk.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                imuk.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                imuk.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                imuk.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        except:
            pass
        # Create and plot RGB UKIDSS image
        # Doesnt work beacuse the UKIDSS images must be read with the second extension not the first one
#        try:
#            aplpy.make_rgb_cube(['K_UKIDSS_TMP.fits.gz','H_UKIDSS_TMP.fits.gz','J_UKIDSS_TMP.fits.gz'],'UKIDSS_rgb.fits')
#            imukc = aplpy.FITSFigure('UKIDSS_rgb_2d.fits',figure=fig,subplot=(nyplot,nxplot,25))
#            imukc.show_rgb('UKIDSS_rgb.png')
#            imukc.hide_tick_labels()
#            imukc.ticks.set_color('k')
#            imukc.ticks.set_minor_frequency(0)
#            imukc.ticks.set_xspacing(0.5/60.0)
#            imukc.ticks.set_yspacing(0.5/60.0)
#            imukc.hide_xaxis_label()
#            imukc.hide_yaxis_label()
#            imukc.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
#        except:
#            pass

# Plot VHS images
    if VHS:
        try:#Y band
            fitsfile = 'Y_VHS_TMP.fits.gz'
            imvhs = aplpy.FITSFigure(fitsfile,figure=fig,subplot=(nyplot,nxplot,15),north=True,hdu=1)
            imvhs.hide_tick_labels()
            imvhs.ticks.set_color('k')
            imvhs.ticks.set_minor_frequency(0)
            imvhs.ticks.set_xspacing(0.5/60.0)
            imvhs.ticks.set_yspacing(0.5/60.0)
            imvhs.hide_xaxis_label()
            imvhs.hide_yaxis_label()
            data1 = pyfits.getdata(fitsfile)
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            imvhs.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            imvhs.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            imvhs.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                imvhs.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            imvhs.add_label(0.2,0.9,'VHS Y',relative=True,size='medium',color='k')
            #Need to get year
            #imvhs.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                imvhs.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                imvhs.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                imvhs.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        except:
            pass
        try:#J band
            fitsfile = 'J_VHS_TMP.fits.gz'
            imvhs = aplpy.FITSFigure(fitsfile,figure=fig,subplot=(nyplot,nxplot,16),north=True,hdu=1)
            imvhs.hide_tick_labels()
            imvhs.ticks.set_color('k')
            imvhs.ticks.set_minor_frequency(0)
            imvhs.ticks.set_xspacing(0.5/60.0)
            imvhs.ticks.set_yspacing(0.5/60.0)
            imvhs.hide_xaxis_label()
            imvhs.hide_yaxis_label()
            data1 = pyfits.getdata(fitsfile)
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            imvhs.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            imvhs.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            imvhs.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                imvhs.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            imvhs.add_label(0.2,0.9,'VHS J',relative=True,size='medium',color='k')
            #Need to get year
            #imvhs.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                imvhs.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                imvhs.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                imvhs.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        except:
            pass
        try:#H band
            fitsfile = 'H_VHS_TMP.fits.gz'
            imvhs = aplpy.FITSFigure(fitsfile,figure=fig,subplot=(nyplot,nxplot,17),north=True,hdu=1)
            imvhs.hide_tick_labels()
            imvhs.ticks.set_color('k')
            imvhs.ticks.set_minor_frequency(0)
            imvhs.ticks.set_xspacing(0.5/60.0)
            imvhs.ticks.set_yspacing(0.5/60.0)
            imvhs.hide_xaxis_label()
            imvhs.hide_yaxis_label()
            data1 = pyfits.getdata(fitsfile)
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            imvhs.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            imvhs.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            imvhs.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                imvhs.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            imvhs.add_label(0.2,0.9,'VHS H',relative=True,size='medium',color='k')
            #Need to get year
            #imvhs.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                imvhs.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                imvhs.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                imvhs.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        except:
            pass
        try:#Ks band
            fitsfile = 'KS_VHS_TMP.fits.gz'
            imvhs = aplpy.FITSFigure(fitsfile,figure=fig,subplot=(nyplot,nxplot,18),north=True,hdu=1)
            imvhs.hide_tick_labels()
            imvhs.ticks.set_color('k')
            imvhs.ticks.set_minor_frequency(0)
            imvhs.ticks.set_xspacing(0.5/60.0)
            imvhs.ticks.set_yspacing(0.5/60.0)
            imvhs.hide_xaxis_label()
            imvhs.hide_yaxis_label()
            data1 = pyfits.getdata(fitsfile)
            med1 = numpy.nanmedian(data1)
            mad1 = numpy.nanmedian(abs(data1 - numpy.nanmedian(data1)))
            min1 = med1 - 2.0*mad1
            max1 = med1 + 10.0*mad1
            imvhs.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
            imvhs.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
            imvhs.show_circles(ra,de,edgecolor='r',linewidth=0.7,facecolor='none',radius=0.0025)
            if secondary:
                imvhs.show_circles(ra2,de2,edgecolor='b',linewidth=0.7,facecolor='none',radius=0.0025)
            imvhs.add_label(0.2,0.9,'VHS Ks',relative=True,size='medium',color='k')
            #Need to get year
            #imvhs.add_label(0.15,0.1,images[i][2][0:4],relative=True,size='medium',color='k')
            if allwise:
                imvhs.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
            if rejallwise:
                imvhs.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
            if tmass:
                imvhs.show_circles(tmass_ra,tmass_de,edgecolor='k',facecolor=tm_color,radius=0.0004,linewidth=0.5)
        except:
            pass
    
    # Create and plot RGB AllWISE image
    aplpy.make_rgb_cube(['AllWISE_w3.fits','AllWISE_w2.fits','AllWISE_w1.fits'],'AllWISE_rgb.fits')
    im20 = aplpy.FITSFigure('AllWISE_rgb_2d.fits',figure=fig,subplot=(nyplot,nxplot,30))
    aplpy.make_rgb_image('AllWISE_rgb.fits','AllWISE_rgb.png',vmin_r=wmin3,vmin_g=wmin2,vmin_b=wmin1,vmax_r=wmax3,vmax_g=wmax2,vmax_b=wmax1)
    
    im20.show_rgb('AllWISE_rgb.png')
    im20.hide_tick_labels()
    im20.ticks.set_color('k')
    im20.ticks.set_minor_frequency(0)
    im20.ticks.set_xspacing(0.5/60.0)
    im20.ticks.set_yspacing(0.5/60.0)
    im20.hide_xaxis_label()
    im20.hide_yaxis_label()
    im20.recenter(ra,de,width=(2.0/60.0),height=(2.0/60.0))
    
    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95,wspace=0.05,hspace=0.05)
    # Add Labels
    ras = deg2str(ra)
    des = deg2str(de,dec=1)
    c1 = coord.ICRS(ra,de,unit=(u.degree,u.degree))
    c2 = c1.galactic
    pylab.figtext(0.60,0.45+(fig_ysize-fig_ysize0)/fig_ysize,r'$\alpha$ = '+ras+'\t('+str(round(ra,6))+')',fontsize=15)
    pylab.figtext(0.60,0.42+(fig_ysize-fig_ysize0)/fig_ysize,r'$\delta$ = '+des+'\t('+str(round(de,6))+')',fontsize=15) 
    pylab.figtext(0.60,0.39+(fig_ysize-fig_ysize0)/fig_ysize,r'$l$ = '+str(round(c2.l.degree,3)),fontsize=15)
    pylab.figtext(0.60,0.36+(fig_ysize-fig_ysize0)/fig_ysize,r'$b$ = '+str(round(c2.b.degree,3)),fontsize=15)
    if addtext:
        pylab.figtext(0.60,0.33+(fig_ysize-fig_ysize0)/fig_ysize,addtext,fontsize=15)
    if addtext2:
        pylab.figtext(0.60,0.30+(fig_ysize-fig_ysize0)/fig_ysize,addtext2,fontsize=14)
    if allwise:
        pylab.figtext(0.741,0.339+(fig_ysize-fig_ysize0)/fig_ysize,'Allwise catalog sources',color='k',fontsize=15)
        pylab.figtext(0.74,0.34+(fig_ysize-fig_ysize0)/fig_ysize,'Allwise catalog sources',color=allcolor,fontsize=15)
    if rejallwise:
        pylab.figtext(0.756,0.309+(fig_ysize-fig_ysize0)/fig_ysize,'Allwise reject sources',color='k',fontsize=15)
        pylab.figtext(0.755,0.31+(fig_ysize-fig_ysize0)/fig_ysize,'Allwise reject sources',color=rejcolor,fontsize=15)
    if tmass:
        pylab.figtext(0.741,0.279+(fig_ysize-fig_ysize0)/fig_ysize,'2MASS catalog sources',color='k',fontsize=15)
        pylab.figtext(0.74,0.28+(fig_ysize-fig_ysize0)/fig_ysize,'2MASS catalog sources',color=tm_color,fontsize=15)


    # Remove files (or not)
    if keepfiles:
        j = 1
    else:
        print "Removing files..."
        cmdrm1 = "rm source.xml 2MASS*.fits AllWISE*.fits DSS*.fits AllWISE_rgb.png *UKIDSS_TMP.fits.gz *VHS_TMP.fits.gz UKIDSS_rgb*.fits UKIDSS_rgb.png"
        os.system(cmdrm1)
        if allwise:
            cmdrm2 = "rm allwise.tbl"
            os.system(cmdrm2)
            if rejallwise:
                cmdrm3 = "rm rejallwise.tbl"
                os.system(cmdrm3)
        if SDSS:
            cmdrm4 = "rm SDSS*.fits"
            os.system(cmdrm4)
        if tmass:
            cmdrm5 = "rm tmass.tbl"
            os.system(cmdrm5)
    if savepdf:
        pylab.savefig(source_name+'.pdf')
    
    t2 = datetime.now()
    tdiff = (t2 - t1)
    print "Finder creation took %s seconds" % (round(tdiff.total_seconds(),0))

finderold('HIP 1')