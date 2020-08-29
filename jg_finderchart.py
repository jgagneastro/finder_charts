import os
import numpy as np
import xml.etree.ElementTree as ET
from deg2str import *
#from astropy import coordinates as coord
from astropy.coordinates import SkyCoord
from astropy import units as u
from simbad import *
from datetime import *
from PIL import Image
from query_wsa_fits import *
from query_pso_fits import *
from jdcal import *
import pdb
import glob
import astropy.io.fits as pyfits

#import astropy.io.fits as aplpy
stop=pdb.set_trace

def finder(source_name,allwise=False,rejallwise=False,tmass=False,rejtmass=False,PSO=True,UKIDSS=True,VHS=True,DES=True,keepfiles=False,allcolor='#FFFF00',rejcolor='b',tm_color='r',plot=False,savepdf=True,secondary='',addtext='',addtext2='',skipdownloads=False,circle_radius=0.0025,size=1.667,override_directory=None,primarypos_label=None,secondarypos_label=None,title=None,filename=None,buffer=False,gnirsacq=False,DSS=True,TMASSIM=True,WISE=True,circle_alpha=.8,labels=True,pos_list_gray_ra=None,pos_list_gray_dec=None,pos_list_gray_sizes=None,pos_list_gray_pmra=None,pos_list_gray_pmdec=None,gray_label=None,pos3=None,pos3_label=None,pos4=None,pos4_label=None,pos5=None,pos5_label=None):
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
	
	#Use buffer if needed
	if buffer:
		import matplotlib
		matplotlib.use('Agg')
		plot=False
	
	#Import those after matplotlib is explicitly set as Agg in case Python is not installed as a Mac OS X framework
	import pylab
	import aplpy

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
	color_green = '#4daf4a'#RGB=[77,175,74]
	color_orange = '#ff7f00'#RGB=[255,127,0]
	color_pink = '#f4d7d7'#RGB=[244,215,215]
	col_yellow = '#ffde02'#RGB=[255,222,2]

	nxplot = 5
	nyplot = 3
	fig_xsize = 11
	fig_ysize = 8.5
	
	if not labels:
		nxplot -= 2
	
	t1 = datetime.now()	

	ra,de = simbad(source_name)
	ra2 = None
	de2 = None
	if secondary:
		ra2,de2 = simbad(secondary)
	ra3 = None
	de3 = None
	if pos3:
		ra3,de3 = simbad(pos3)
	ra4 = None
	de4 = None
	if pos4:
		ra4,de4 = simbad(pos4)
	ra5 = None
	de5 = None
	if pos5:
		ra5,de5 = simbad(pos5)
	if filename is None:
		filename = source_name
	
	#Download xml file from IRSA
	xmlfile = "source.xml"
	if skipdownloads is not True:
		print("Getting xml file...")
		cmd = "wget -O "+xmlfile+" 'http://irsa.ipac.caltech.edu/applications/finderchart/servlet/api?locstr="+str(ra)+"+"+str(de)+"&subsetsize="+str(size)+"' "
		os.system(cmd)
	
	# parse xml file
	print("Parsing xml file...")
	
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
		if not DSS and (surveyname == 'DSS' or surveyname == 'DSS1' or surveyname == 'DSS2'):
			continue
		if not TMASSIM and surveyname == '2MASS':
			continue
		if not WISE and (surveyname == 'WISE' or surveyname == 'WISE (AllWISE)'):
			continue
		images.append([surveyname,band,obsdate,fitsurl])
	
	if skipdownloads is not True:
		if DSS:
			print("Downloading DSS data...")
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
		
		if TMASSIM:
			print("Downloading 2MASS data...")
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
		
		if WISE:
			print("Downloading WISE data")
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
			print("Downloading UKIDSS data")
			#Remove previous data
			os.system("rm *_UKIDSS_TMP.fits.gz")
			query_wsa_fits(ra,de,size=size,output_file='UKIDSS_TMP.fits.gz',filter='all',catalog='UKIDSS')
		
		if VHS:
			print("Downloading VHS data")
			#Remove previous data
			os.system("rm *_VHS_TMP.fits.gz")
			query_wsa_fits(ra,de,size=size,output_file='VHS_TMP.fits.gz',filter='all',catalog='VHS')
		
		if PSO:
			print("Downloading Pan-Starrs data")
			#Remove previous data
			os.system("rm *_PSO_TMP.fits*")
			query_pso_fits(ra,de,size=size,output_file='PSO_TMP.fits')

		if DES:
			print("Downloading DES DR1 data")
			#Remove previous data
			os.system("rm *_DES_TMP.fits*")
			query_des_fits(ra,de,size=size,output_file='DES_TMP.fits')
			#from astropy.coordinates import SkyCoord
			#from hips import WCSGeometry
			#from hips import make_sky_image

			#geometry = WCSGeometry.create(skydir=SkyCoord(82.418457, -46.987488, unit='deg', frame='icrs'),width=500, height=500, fov="0.03 deg",coordsys='icrs', projection='AIT')
			#hips_survey = 'CDS/P/DES-DR1/Y'
			#result = make_sky_image(geometry, hips_survey, 'fits')
			#result.write_image('my_image3.fits')
	
	#If no UKIDSS data could be downloaded, turn off the UKIDSS option
	if len(glob.glob('*_UKIDSS_TMP.fits*')) == 0:
		UKIDSS = False
	
	#If no VHS data could be downloaded, turn off the VHS option
	if len(glob.glob('*_VHS_TMP.fits*')) == 0:
		VHS = False
	
	#If no PSO data could be downloaded, turn off the PSO option
	if len(glob.glob('*_PSO_TMP.fits*')) == 0:
		PSO = False
	
	if PSO:
		nxplot = np.maximum(nxplot,5)
	if UKIDSS:
		nxplot = np.maximum(nxplot,4)
	
	#Determine the amount of additional rows needed
	vertical_spacing = 0
	ukidss_spacing = 0
	vhs_spacing = 0
	tmass_spacing = 0
	allwise_spacing = 1
	dss_negspacing = 0
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
	if not DSS:
		vertical_spacing -= 1
		#dss_negspacing -= 1
		ukidss_spacing -= 1
		vhs_spacing -= 1
		tmass_spacing -= 1
		allwise_spacing -= 1
	if not WISE:
		vertical_spacing -= 1
	
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
			print("No AllWISE sources found!")
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
			print("No AllWISE Reject sources found!")
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
			print("No 2MASS sources found!")
			tmass=False
	
	rejtmass_ra = None
	rejtmass_de = None
	rejtmass_ph = None
	rejtmass_rel = None
	if rejtmass:
		if skipdownloads is not True:
			cmd4 = 'curl -o rejtmass.tbl "http://irsa.ipac.caltech.edu/TAP/sync?FORMAT=IPAC_TABLE&QUERY=SELECT+ra,dec,rel,ph_qual+FROM+pt_src_rej+WHERE+CONTAINS(POINT(\'J2000\',ra,dec),CIRCLE(\'J2000\','+str(ra)+','+str(de)+','+str(size/60.0)+'))=1"'
			#This is a temporary work-out, but only fetches rel='A' entries
			#cmd4 = 'wget -O rejtmass.tbl "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=Box&size='+str(size*60.0)+'&radunits=arcsec&objstr='+str(ra)+','+str(de)+'&catalog=pt_src_rej&selcols=ra,dec&outfmt=1"'
			os.system(cmd4)
		try:
			rejtmass_psc = np.loadtxt('rejtmass.tbl',skiprows=16,unpack=True,usecols=(0,1))
			rejtmass_psc2 = np.loadtxt('rejtmass.tbl',skiprows=16,unpack=True,usecols=(2,3),dtype='str')
			rejtmass_ra,rejtmass_de = rejtmass_psc
			rejtmass_rel,rejtmass_ph = rejtmass_psc2
		except:
			print("No 2MASS-Reject sources found!")
			rejtmass=False
	
	if plot:
		pylab.ion()
	fig = pylab.figure(figsize=(fig_xsize,fig_ysize))
	pylab.rcParams['font.family'] = 'serif'
	
	for i in range(len(images)):
		if images[i][1] == 'DSS1 Blue' and DSS:
			min1, max1, im = oplotfits(fig,'DSS1_Blue.fits',nyplot,nxplot,1,ra,de,'DSS1 B',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)

		if images[i][1] == 'DSS1 Red' and DSS:
			oplotfits(fig,'DSS1_Red.fits',nyplot,nxplot,2,ra,de,'DSS1 R',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
	
		if images[i][1] == 'DSS2 Blue' and DSS:
			min1, max1, im = oplotfits(fig,'DSS2_Blue.fits',nyplot,nxplot,3,ra,de,'DSS2 B',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
			if title:
				im.show_xaxis_label()
				im.axis_labels.set_xposition('top')
				im.axis_labels.set_xtext(title)
	
		if images[i][1] == 'DSS2 Red' and DSS:
			oplotfits(fig,'DSS2_Red.fits',nyplot,nxplot,4,ra,de,'DSS2 R',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		
		if images[i][1] == 'DSS2 IR' and DSS:
			oplotfits(fig,'DSS2_IR.fits',nyplot,nxplot,5,ra,de,'DSS2 IR',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		
		#The "gray" positions are overrided by 2MASS-Reject
		rejtmass_ra_sub = None
		rejtmass_de_sub = None
		rejtmass_pmra_sub = None
		rejtmass_pmde_sub = None
		rejtmass_rel_sub = None
		rejtmass_sizes = None
		rejtmass_sizes_sub = None
		if pos_list_gray_ra is not None and pos_list_gray_dec is not None:
			rejtmass_ra_sub = np.array(pos_list_gray_ra)
			rejtmass_de_sub = np.array(pos_list_gray_dec)
		if pos_list_gray_sizes is not None:
			rejtmass_sizes_sub = np.array(pos_list_gray_sizes)
		if pos_list_gray_pmra is not None and pos_list_gray_pmdec is not None:
			rejtmass_pmra_sub = np.array(pos_list_gray_pmra)
			rejtmass_pmde_sub = np.array(pos_list_gray_pmdec)

		if images[i][1] == 'J' and TMASSIM:
			
			#If 2MASS-Reject sources are to be displayed, only display those detected in the appropriate band
			if rejtmass:
				qual = [t[0] for t in rejtmass_ph]
				goodqual = np.where(np.array(qual) != 'U')
				rejtmass_ra_sub = rejtmass_ra[goodqual]
				rejtmass_de_sub = rejtmass_de[goodqual]
				rejtmass_rel_sub = rejtmass_rel[goodqual]
				#if rejtmass_sizes is not None:
				#	rejtmass_sizes_sub = rejtmass_sizes[goodqual]
			
			void1, void2, im = oplotfits(fig,'2MASS_J.fits',nyplot,nxplot,nxplot+1+tmass_spacing*nxplot+dss_negspacing*nxplot,ra,de,'2MASS $J$',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,rejtmass_ra=rejtmass_ra_sub,rejtmass_de=rejtmass_de_sub,rejtmass_sizes=rejtmass_sizes_sub,rejtmass_pmra=rejtmass_pmra_sub,rejtmass_pmde=rejtmass_pmde_sub,circle_radius=circle_radius,size=size,rejtmass_rel=rejtmass_rel_sub,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
			
			if not WISE:
				im.show_xaxis_label()
				im.show_yaxis_label()
				im.axis_labels.set_xtext('RA (East left)')
				im.axis_labels.set_ytext('Dec (North up)')
		
		if images[i][1] == 'H' and TMASSIM:
			
			#If 2MASS-Reject sources are to be displayed, only display those detected in the appropriate band
			if rejtmass:
				qual = [t[1] for t in rejtmass_ph]
				goodqual = np.where(np.array(qual) != 'U')
				rejtmass_ra_sub = rejtmass_ra[goodqual]
				rejtmass_de_sub = rejtmass_de[goodqual]
				rejtmass_rel_sub = rejtmass_rel[goodqual]
				#if rejtmass_sizes is not None:
				#	rejtmass_sizes_sub = rejtmass_sizes[goodqual]

			oplotfits(fig,'2MASS_H.fits',nyplot,nxplot,nxplot+2+tmass_spacing*nxplot+dss_negspacing*nxplot,ra,de,'2MASS $H$',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,rejtmass_ra=rejtmass_ra_sub,rejtmass_de=rejtmass_de_sub,rejtmass_sizes=rejtmass_sizes_sub,rejtmass_pmra=rejtmass_pmra_sub,rejtmass_pmde=rejtmass_pmde_sub,rejtmass_rel=rejtmass_rel_sub,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		
		if images[i][1] == 'K' and TMASSIM:
			
			#If 2MASS-Reject sources are to be displayed, only display those detected in the appropriate band
			if rejtmass:
				qual = [t[2] for t in rejtmass_ph]
				goodqual = np.where(np.array(qual) != 'U')
				rejtmass_ra_sub = rejtmass_ra[goodqual]
				rejtmass_de_sub = rejtmass_de[goodqual]
				rejtmass_rel_sub = rejtmass_rel[goodqual]
				#if rejtmass_sizes is not None:
				#	rejtmass_sizes_sub = rejtmass_sizes[goodqual]

			void1, void2, im = oplotfits(fig,'2MASS_K.fits',nyplot,nxplot,nxplot+3+tmass_spacing*nxplot+dss_negspacing*nxplot,ra,de,'2MASS $K_S$',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,rejtmass_ra=rejtmass_ra_sub,rejtmass_de=rejtmass_de_sub,rejtmass_sizes=rejtmass_sizes_sub,rejtmass_pmra=rejtmass_pmra_sub,rejtmass_pmde=rejtmass_pmde_sub,rejtmass_rel=rejtmass_rel_sub,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
			
			if not labels:
				xlabels = .95
				ytoplabels = .89
				ydeltalabels = .12
				if primarypos_label:	
					im.add_label(xlabels,ytoplabels,primarypos_label,relative=True,size='medium',color='r',bbox=dict(facecolor='white', alpha=0.5),horizontalalignment='right')
				if secondary and secondarypos_label:
					im.add_label(xlabels,ytoplabels-ydeltalabels,secondarypos_label,relative=True,size='medium',color='b',bbox=dict(facecolor='white', alpha=0.5),horizontalalignment='right')
				if gnirsacq:
					im.add_label(xlabels,ytoplabels-2*ydeltalabels,'GNIRS Acq',relative=True,size='medium',color='g',bbox=dict(facecolor='white', alpha=0.5),horizontalalignment='right')
		# else:
		# 	im = fig.add_subplot(nyplot,nxplot,nxplot+3+tmass_spacing*nxplot+dss_negspacing*nxplot)
		# 	im.axis('off')
		# 	if not labels:
		# 		xlabels = .95
		# 		ytoplabels = .89
		# 		ydeltalabels = .12
		# 		if primarypos_label:	
		# 			im.annotate(primarypos_label,xy=(xlabels,ytoplabels),fontsize=15,color='r',horizontalalignment='right')
		# 		if secondary and secondarypos_label:
		# 			im.annotate(secondarypos_label,xy=(xlabels,ytoplabels-ydeltalabels),fontsize=15,color='b',horizontalalignment='right')
		# 		if gnirsacq:
		# 			im.annotate('GNIRS Acq',xy=(xlabels,ytoplabels-2*ydeltalabels),fontsize=15,color='g',horizontalalignment='right')
		
		if images[i][1] == 'w1' and WISE:
			wmin1, wmax1, im = oplotfits(fig,'AllWISE_w1.fits',nyplot,nxplot,nxplot+1+allwise_spacing*nxplot+dss_negspacing*nxplot,ra,de,'W1',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
			im.show_xaxis_label()
			im.show_yaxis_label()
			im.axis_labels.set_xtext('RA (East left)')
			im.axis_labels.set_ytext('Dec (North up)')

		if images[i][1] == 'w2' and WISE:
			wmin2, wmax2, void = oplotfits(fig,'AllWISE_w2.fits',nyplot,nxplot,nxplot+2+allwise_spacing*nxplot+dss_negspacing*nxplot,ra,de,'W2',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		
		if images[i][1] == 'w3' and WISE:
			wmin3, wmax3, void = oplotfits(fig,'AllWISE_w3.fits',nyplot,nxplot,nxplot+3+allwise_spacing*nxplot+dss_negspacing*nxplot,ra,de,'W3',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		
		if images[i][1] == 'w4' and WISE:
			wmin4, wmax4, void = oplotfits(fig,'AllWISE_w4.fits',nyplot,nxplot,nxplot+4+allwise_spacing*nxplot+dss_negspacing*nxplot,ra,de,'W4',year=images[i][2][0:4],ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
	
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
			wminpsog, wmaxpsog, void = oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+1+dss_negspacing*nxplot,ra,de,'PSO $g$',year=year,ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
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
			wminpsor, wmaxpsor, void = oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+2+dss_negspacing*nxplot,ra,de,'PSO $r$',year=year,ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
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
			wminpsoi, wmaxpsoi, void = oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+3+dss_negspacing*nxplot,ra,de,'PSO $i$',year=year,ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
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
			wminpsoz, wmaxpsoz, void = oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+4+dss_negspacing*nxplot,ra,de,'PSO $z$',year=year,ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
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
			wminpsoy, wmaxpsoy, void = oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+5+dss_negspacing*nxplot,ra,de,'PSO $y$',year=year,ra2=ra2,de2=de2,north=False,hdu=0,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
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
			oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+1+ukidss_spacing*nxplot+dss_negspacing*nxplot,ra,de,'UKIDSS Y',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		except:
			pass
		try:#J band
			#Get date and then make plot
			fitsfile = 'J_UKIDSS_TMP.fits.gz'
			hdulist = pyfits.open(fitsfile)
			year = hdulist[0].header['UTDATE'][0:4]
			hdulist.close()
			oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+2+ukidss_spacing*nxplot+dss_negspacing*nxplot,ra,de,'UKIDSS J',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		except:
			pass
		try:#H band
			#Get date and then make plot
			fitsfile = 'H_UKIDSS_TMP.fits.gz'
			hdulist = pyfits.open(fitsfile)
			year = hdulist[0].header['UTDATE'][0:4]
			hdulist.close()
			oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+3+ukidss_spacing*nxplot+dss_negspacing*nxplot,ra,de,'UKIDSS H',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		except:
			pass
		try:#K band
			#Get date and then make plot
			fitsfile = 'K_UKIDSS_TMP.fits.gz'
			hdulist = pyfits.open(fitsfile)
			year = hdulist[0].header['UTDATE'][0:4]
			hdulist.close()
			oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+4+ukidss_spacing*nxplot+dss_negspacing*nxplot,ra,de,'UKIDSS K',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
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
			oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+1+vhs_spacing*nxplot+dss_negspacing*nxplot,ra,de,'VHS $Y$',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		except:
			pass
		try:#J band
			#Get date and then make plot
			fitsfile = 'J_VHS_TMP.fits.gz'
			hdulist = pyfits.open(fitsfile)
			year = hdulist[1].header['DATE'][0:4]
			hdulist.close()
			oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+2+vhs_spacing*nxplot+dss_negspacing*nxplot,ra,de,'VHS $J$',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		except:
			pass
		try:#H band
			#Get date and then make plot
			fitsfile = 'H_VHS_TMP.fits.gz'
			hdulist = pyfits.open(fitsfile)
			year = hdulist[1].header['DATE'][0:4]
			hdulist.close()
			oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+3+vhs_spacing*nxplot+dss_negspacing*nxplot,ra,de,'VHS $H$',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		except:
			pass
		try:#Ks band
			#Get date and then make plot
			fitsfile = 'Ks_VHS_TMP.fits.gz'
			hdulist = pyfits.open(fitsfile)
			year = hdulist[1].header['DATE'][0:4]
			hdulist.close()
			oplotfits(fig,fitsfile,nyplot,nxplot,nxplot+4+vhs_spacing*nxplot+dss_negspacing*nxplot,ra,de,'VHS $K_S$',year=year,ra2=ra2,de2=de2,north=True,hdu=1,allwise=allwise,rejallwise=rejallwise,tmass=tmass,allcolor=allcolor,rejcolor=rejcolor,tm_color=tm_color,secondary=secondary,allwise_ra=allwise_ra,allwise_de=allwise_de,rejallwise_ra=rejallwise_ra,rejallwise_de=rejallwise_de,tmass_ra=tmass_ra,tmass_de=tmass_de,circle_radius=circle_radius,size=size,buffer=buffer,gnirsacq=gnirsacq,circle_alpha=circle_alpha,ra3=ra3,de3=de3,ra4=ra4,de4=de4,ra5=ra5,de5=de5)
		except:
			pass
	
	# Create and plot RGB AllWISE image
	if WISE:
		try:
			files = ['AllWISE_w3.fits','AllWISE_w2.fits','AllWISE_w1.fits']
			aplpy.make_rgb_cube(files,'AllWISE_rgb.fits')
			imawrgb = aplpy.FITSFigure('AllWISE_rgb_2d.fits',figure=fig,subplot=(nyplot,nxplot,nxplot+5+allwise_spacing*nxplot+dss_negspacing*nxplot))
			imawrgb.add_label(0.05,0.9,'W3/W2/W1',relative=True,size='medium',color='k',bbox=dict(facecolor='white', alpha=0.5),horizontalalignment='left')
			
			meds = []
			mads = []
			devs = []
			for filei in files:
				datai = pyfits.getdata(filei)
				medi = np.nanmedian(datai)
				madi = np.nanmedian(abs(datai - np.nanmedian(datai)))
				devi = np.nanpercentile(datai,90) - np.nanpercentile(datai,10)
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
			
			aplpy.make_rgb_image('AllWISE_rgb.fits','AllWISE_rgb.png',vmin_r=mins[0],vmin_g=mins[1],vmin_b=mins[2],vmax_r=maxs[0],vmax_g=maxs[1],vmax_b=maxs[2])
			
			imawrgb.show_rgb('AllWISE_rgb.png')
			imawrgb.hide_tick_labels()
			imawrgb.ticks.hide()
			imawrgb.axis_labels.set_xtext('Size = '+str(size)+' arcmin')
			imawrgb.hide_yaxis_label()
			imawrgb.recenter(ra,de,width=(size/60.0),height=(size/60.0))
			circle_width = 1.2
			imawrgb.show_circles(ra,de,edgecolor=color_red,linewidth=circle_width,facecolor='none',radius=circle_radius,alpha=circle_alpha)
			if secondary:
				imawrgb.show_circles(ra2,de2,edgecolor=color_blue,linewidth=circle_width,facecolor='none',radius=circle_radius,alpha=circle_alpha)
		except:
			pass
	
	#Create and plot RGB PSO image only if there's enough space (vertical_spacing >= 2)
	if PSO and vertical_spacing >= 2:
		try:
			files = ['y_PSO_TMP.fits','i_PSO_TMP.fits','g_PSO_TMP.fits']
			aplpy.make_rgb_cube(files,'PSO_rgb.fits')
			impsoc = aplpy.FITSFigure('PSO_rgb_2d.fits',figure=fig,subplot=(nyplot,nxplot,nxplot+10+dss_negspacing*nxplot))
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
			#aplpy.make_rgb_image('PSO_rgb.fits','PSO_rgb.png')
			
			impsoc.show_rgb('PSO_rgb.png')
			impsoc.hide_tick_labels()
			impsoc.ticks.hide()
			impsoc.hide_xaxis_label()
			impsoc.hide_yaxis_label()
			impsoc.recenter(ra,de,width=(size/60.0),height=(size/60.0))
			impsoc.show_circles(ra,de,edgecolor=color_red,linewidth=circle_width,facecolor='none',radius=circle_radius,alpha=circle_alpha)
			if secondary:
				impsoc.show_circles(ra2,de2,edgecolor=color_blue,linewidth=circle_width,facecolor='none',	radius=circle_radius,alpha=circle_alpha)
		except:
			pass
	
	bottom_space = 0.05
	top_space = 0.05
	pylab.subplots_adjust(left=0.05,right=0.95,bottom=bottom_space,top=1.0-top_space,wspace=0.05,hspace=0.05)

	# Add Labels
	if labels:
		ras = deg2str(ra)
		des = deg2str(de,dec=1)
		#c1 = coord.ICRS(ra*u.degree,de*u.degree)
		c1 = SkyCoord(ra=ra*u.degree, dec=de*u.degree, frame='icrs')
		c2 = c1.galactic
		
		sptext = fig.add_subplot(nyplot,nxplot,nxplot+4+tmass_spacing*nxplot+dss_negspacing*nxplot)
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
		
		sptext2 = fig.add_subplot(nyplot,nxplot,nxplot+5+tmass_spacing*nxplot+dss_negspacing*nxplot)
		sptext2.axis('off')
		if primarypos_label:
			sptext2.annotate(primarypos_label,xy=(0.99,ytoplabels-ydeltalabels*2),fontsize=15,color='r',horizontalalignment='right')
		if secondary and secondarypos_label:
			sptext2.annotate(secondarypos_label,xy=(0.99,ytoplabels-ydeltalabels*3),fontsize=15,color='b',horizontalalignment='right')
		if pos3 is not None and pos3_label is not None:
			sptext2.annotate(pos3_label,xy=(0.99,ytoplabels-ydeltalabels*4),fontsize=15,color=color_purple,horizontalalignment='right')
		if pos4 is not None and pos4_label is not None:
			sptext2.annotate(pos4_label,xy=(0.99,ytoplabels-ydeltalabels*5),fontsize=15,color=color_orange,horizontalalignment='right')
		if pos5 is not None and pos5_label is not None:
			sptext2.annotate(pos5_label,xy=(0.99,ytoplabels-ydeltalabels*6),fontsize=15,color=color_pink,horizontalalignment='right')
		if gray_label is not None:
			sptext2.annotate(gray_label,xy=(0.99,ytoplabels-ydeltalabels*7),fontsize=15,color=color_green,horizontalalignment='right')
		if gnirsacq:
			sptext2.annotate('GNIRS Acq',xy=(0.99,ytoplabels-ydeltalabels*7),fontsize=15,color='g',horizontalalignment='right')

	# Remove files (or not)
	if keepfiles:
		pass
	else:
		print("Removing files...")
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
		pylab.savefig(filename+'.pdf')
	
	#Return to initial directory
	if main_dir:
		os.chdir(initial_dir)
	
	t2 = datetime.now()
	tdiff = (t2 - t1)
	print("Finder creation took %s seconds" % (round(tdiff.total_seconds(),0)))

#This function displays a fits image with the appropriate annotations
def oplotfits(fig,fitsfile,nyplot,nxplot,position,ra,de,label,year='',xlabel=0.05,ra2=None,de2=None,north=False,hdu=0,allwise=False,rejallwise=False,tmass=False,allcolor='#FFFF00',rejcolor='b',tm_color='r',secondary='',allwise_ra=None,allwise_de=None,rejallwise_ra=None,rejallwise_de=None,tmass_ra=None,tmass_de=None,rejtmass_ra=None,rejtmass_de=None,rejtmass_sizes=None,rejtmass_pmra=None,rejtmass_pmde=None,circle_radius=0.0025,size=2.0,rejtmass_rel=None,buffer=False,gnirsacq=False,circle_alpha=0.8,ra3=None,de3=None,ra4=None,de4=None,ra5=None,de5=None):
	
	#Use buffer if needed
	if buffer:
		import matplotlib
		matplotlib.use('Agg')
	import pylab,pyfits,aplpy
	
	#Circle parameters
	circle_width = 1.2

	#List of colors
	color_blue = '#377eb8'#RGB=[55,126,184]
	color_red = '#e41a1c'#RGB=[228,26,28]
	color_purple = '#b27bba'#RGB=[178,123,186]
	color_green = '#4daf4a'#RGB=[77,175,74]
	color_orange = '#ff7f00'#RGB=[255,127,0]
	color_pink = '#f4d7d7'#RGB=[244,215,215]
	col_yellow = '#ffde02'#RGB=[255,222,2]

	#Display FITS image
	im = aplpy.FITSFigure(fitsfile,figure=fig,subplot=(nyplot,nxplot,position),north=north,hdu=hdu)

	#Hide ticks and labels
	im.hide_tick_labels()
	im.ticks.hide()
	im.hide_xaxis_label()
	im.hide_yaxis_label()
	
	#Read the FITS data in order to determine appropriate scaling
	data1 = pyfits.getdata(fitsfile)
	med1 = np.nanmedian(data1)
	mad1 = np.nanmedian(abs(data1 - np.nanmedian(data1)))
	min1 = med1 - 2.0*mad1
	max1 = med1 + 10.0*mad1
	
	#Fix the scaling
	im.show_colorscale(cmap='gist_yarg',aspect='equal',vmax=max1,vmin=min1)
	
	#Recenter and zoom FITS image
	im.recenter(ra,de,width=(size/60.0),height=(size/60.0))
	
	#Display a circle at the input (RA,DEC) position
	im.show_circles(ra,de,edgecolor=color_red,linewidth=circle_width,facecolor='none',radius=circle_radius,alpha=circle_alpha)
	
	#Display a secondary circle if needed
	if secondary:
		im.show_circles(ra2,de2,edgecolor=color_blue,linewidth=circle_width,facecolor='none',radius=circle_radius*0.9,alpha=circle_alpha)
	
	if ra3 is not None and de3 is not None:
		im.show_circles(ra3,de3,edgecolor=color_purple,linewidth=circle_width,facecolor='none',radius=circle_radius*0.6,alpha=circle_alpha)
	if ra4 is not None and de4 is not None:
		im.show_circles(ra4,de4,edgecolor=color_orange,linewidth=circle_width,facecolor='none',radius=circle_radius*0.5,alpha=circle_alpha)
	if ra5 is not None and de5 is not None:
		im.show_circles(ra5,de5,edgecolor=color_pink,linewidth=circle_width,facecolor='none',radius=circle_radius*0.4,alpha=circle_alpha)
	
	#Display GNIRS acquisition field if needed
	if gnirsacq:
		gnirs_rect_side = 100.
		gnirs_rect_height = 10.
		gnirs_circle_size = 15.0
		intercept = 2.0*np.sqrt(gnirs_circle_size**2-(gnirs_rect_height/2.0)**2)
		intercept_angle = np.arcsin(gnirs_rect_height/2.0/gnirs_circle_size)
		fra = 1.0/(2.0*3600.0*np.cos(np.deg2rad(de)))
		fde = 1.0/(2.0*3600.0)
		min_angle = intercept_angle
		#Circle pointint up
		max_angle = np.pi-intercept_angle
		circ_angles = np.linspace(max_angle,min_angle,num=50)
		#Circle pointint down (wrong)
		#min_angle = np.pi+intercept_angle
		#max_angle = 2*np.pi-intercept_angle
		#circ_angles = np.linspace(max_angle,min_angle,num=50)
		circ_ra = ra+gnirs_circle_size*np.cos(circ_angles)*2*fra
		circ_de = de+gnirs_circle_size*np.sin(circ_angles)*2*fde
		#Circle pointint up
		gnirs_pol1 = np.array([[ra-gnirs_rect_side*fra,de-gnirs_rect_height*fde],[ra-gnirs_rect_side*fra,de+gnirs_rect_height*fde],[ra-intercept*fra,de+gnirs_rect_height*fde]])
		gnirs_pol2 = np.array([[ra+intercept*fra,de+gnirs_rect_height*fde],[ra+gnirs_rect_side*fra,de+gnirs_rect_height*fde],[ra+gnirs_rect_side*fra,de-gnirs_rect_height*fde],[ra-gnirs_rect_side*fra,de-gnirs_rect_height*fde]])
		#Circle pointint down (wrong)
		#gnirs_pol1 = np.array([[ra-gnirs_rect_side*fra,de-gnirs_rect_height*fde],[ra-gnirs_rect_side*fra,de+gnirs_rect_height*fde],[ra+gnirs_rect_side*fra,de+gnirs_rect_height*fde],[ra+gnirs_rect_side*fra,de-gnirs_rect_height*fde],[ra+intercept*fra,de-gnirs_rect_height*fde]])
		#gnirs_pol2 = np.array([[ra-intercept*fra,de-gnirs_rect_height*fde],[ra-gnirs_rect_side*fra,de-gnirs_rect_height*fde]])
		gnirs_pol = np.array([np.concatenate((gnirs_pol1[:,0],circ_ra,gnirs_pol2[:,0])),np.concatenate((gnirs_pol1[:,1],circ_de,gnirs_pol2[:,1]))])
		im.show_polygons([gnirs_pol], edgecolor='green',alpha=0.6,linewidth=2)
	
	#Add the plot label (usually, the name of the survey)
	im.add_label(xlabel,0.9,label,relative=True,size='medium',color='k',bbox=dict(facecolor='white', alpha=0.5),horizontalalignment='left')
	
	#Add year of image if specified
	if year:
		im.add_label(xlabel,0.1,year,relative=True,size='medium',color='k',bbox=dict(facecolor='white', alpha=0.5),horizontalalignment='left')
	
	#Add AllWISE source positions if specified
	if allwise:
		im.show_circles(allwise_ra,allwise_de,edgecolor='k',facecolor=allcolor,radius=0.0004,linewidth=0.5)
	
	#Add AllWISE-Reject source positions if specified
	if rejallwise:
		im.show_circles(rejallwise_ra,rejallwise_de,edgecolor='k',facecolor=rejcolor,radius=0.0004,linewidth=0.5)
	
	#Add 2MASS source positions if specified
	if tmass:
		im.show_circles(tmass_ra,tmass_de,edgecolor=color_green,radius=0.0015,linewidth=1.5)
	
	#Add 2MASS-Reject source positions if specified
	if (rejtmass_ra is not None) and (rejtmass_de is not None):
		#If source reliabilities are specified, use dashed circles for reliability 'F'
		if (rejtmass_rel is not None):
			#pdb.set_trace()
			ind_F = np.where(rejtmass_rel == 'F')
			ind_notF = np.where(rejtmass_rel != 'F')
			if len(ind_F[0]) != 0:
				im.show_circles(rejtmass_ra[ind_F],rejtmass_de[ind_F],edgecolor=color_green,radius=0.001,linewidth=0.5,alpha=circle_alpha,linestyle='--')
			if len(ind_notF[0]) != 0:
				im.show_circles(rejtmass_ra[ind_notF],rejtmass_de[ind_notF],edgecolor=color_green,radius=0.001,linewidth=0.5,alpha=circle_alpha)
		#Otherwise just show regular circles for all 2MASS-Reject entries
		else:
			#If variable sizes are given
			if (rejtmass_sizes is not None):
				for x, y, size in zip(rejtmass_ra, rejtmass_de, rejtmass_sizes):
					im.show_circles(x, y, edgecolor=color_green, radius=0.001*size,linewidth=0.5,alpha=circle_alpha)
			#Otherwise
			else:
				im.show_circles(rejtmass_ra,rejtmass_de,edgecolor=color_green,radius=0.001,linewidth=0.5,alpha=circle_alpha)
		
		#Display proper motions if they were specified
		nyears = 17.5
		if (rejtmass_pmra is not None) and (rejtmass_pmde is not None):
			for i_ra, i_de, i_pmra, i_pmde in zip(rejtmass_ra, rejtmass_de, rejtmass_pmra, rejtmass_pmde):
				#Calculate RA displacement in degrees during 10 years
				#pmra is in mas/yr
				dra = i_pmra/(3600.*1000.*np.cos(np.deg2rad(i_de)))*nyears
				dde = i_pmde/(3600.*1000.)*nyears
				if dra+dde != 0:
					im.show_arrows(i_ra-dra, i_de-dde, dra, dde, edgecolor=color_green,linewidth=0.5,alpha=circle_alpha, head_length=0.5, head_width=0.5)

		#If possible, also display 2MASS main catalog sources if rejtmass is set
		if (not tmass) and (tmass_ra is not None) and (tmass_de is not None):
			im.show_circles(tmass_ra,tmass_de,edgecolor=color_green,radius=0.0015,linewidth=1.5,alpha=circle_alpha)
	#Pass the scaling and plot reference ID
	return min1, max1, im