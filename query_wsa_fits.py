import subprocess
import numpy as np
import pdb

#Downloads a UKIDSS fits image
#size is in arcmin, default 10 arcmin
def query_wsa_fits(ra,dec,size=10.0,catalog='UKIDSS',output_file='',filter='all'):
	
    #Determine query appropriate for catalog
	if catalog=='UKIDSS':
	    #Just Added &obsType=object
	    curl = 'curl -X POST -d "ra='+str(ra)+'&dec='+str(dec)+'&filterID='+filter+'&sys=J2000&xsize='+str(size)+'&ysize='+str(size)+'&obsType=object&frameType=stack&database=UKIDSSDR10PLUS&programmeID=all&submit=Submit" http://wsa.roe.ac.uk:8080/wsa/GetImage'
	if catalog=='VHS':
	    curl = 'curl -X POST -d "database=VHSDR6&programmeID=110&ra='+str(ra)+'&dec='+str(dec)+'&sys=J&filterID='+filter+'&xsize='+str(size)+'&ysize='+str(size)+'&obsType=object&frameType=tilestack&mfid=&fsid=&archive=VSA" http://horus.roe.ac.uk:8080/vdfs/GetImage'
#   if output_file=='':
#        output_file = catalog+'_tmp.fits.gz'
    #Send initial query to WSA
	proc = subprocess.Popen([curl], stdout=subprocess.PIPE,shell=True)
	(out, err) = proc.communicate()
	
	#Split string in a list and search for fits image URLs & band names
	outlist = out.split('\n')
	noutlist = np.size(outlist)
	URL_list = []
	bands_list = []
	for i in range(0L,noutlist):
		posi = outlist[i].find('href="')
		if posi == -1:
			continue
		#Avoid downloading the confidence maps
		posi = outlist[i].find('_conf.fit')
		if posi != -1:
			continue
		typei = outlist[i].split()
		URL_i = outlist[i].split('href="')[1].split('"')[0]
		URL_list.append(URL_i)
		band_i = outlist[i].split("&band=")[1].split("&ra=")[0]
		bands_list.append(band_i)
	
	#Download fits images
	nURL = np.size(URL_list)
	for i in range(0L,nURL):
		sub_curl = 'curl -X POST -d "'+URL_list[i].split("?")[1]+'" '+URL_list[i].split("?")[0]
		sub_proc = subprocess.Popen([sub_curl], stdout=subprocess.PIPE,shell=True)
		(sub_out, sub_err) = sub_proc.communicate()
		sub_URL = sub_out.split('href="')[1].split('"')[0]
		print("Downloading band "+bands_list[i]+": "+bands_list[i]+'_'+output_file)
		sub_curl2 = 'curl -o '+bands_list[i]+'_'+output_file+' -X POST -d "'+sub_URL.split("?")[1]+'" '+sub_URL.split("?")[0]
		sub_proc2 = subprocess.Popen([sub_curl2], stdout=subprocess.PIPE,shell=True)
		(sub_out2, sub_err2) = sub_proc2.communicate()
    