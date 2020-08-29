import subprocess
import numpy as np
import pdb

#Downloads a Pan-Starrs fits image
#size is in arcmin, default 10 arcmin
def query_pso_fits(ra,dec,size=10.0,output_file='PSO_TMP.fits'):
	
	curl = 'curl -X POST -d "pos='+str(ra)+'+'+str(dec)+'&filter=color&filter=g&filter=r&filter=i&filter=z&filter=y&filetypes=stack&auxiliary=data&size='+str(int(round(size*4*60)))+'&output_size=0&verbose=0&autoscale=99.500000&catlist=&submit=Submit" http://ps1images.stsci.edu/cgi-bin/ps1cutouts'
	
	#Send initial query to STSCI
	proc = subprocess.Popen([curl], stdout=subprocess.PIPE,shell=True)
	(out, err) = proc.communicate()
	
	#Split string in a list and search for fits image URLs & band names
	outlist = out.decode().split('\n')
	noutlist = np.size(outlist)
	URL_list = []
	bands_list = []
	for i in range(0,noutlist):
		posi = outlist[i].find('href="')
		if posi == -1:
			continue
		#Only download the FITS cutouts
		posi = outlist[i].find('"Download FITS cutout"')
		if posi == -1:
			continue
		URL_i = outlist[i].split('"Download FITS cutout" href="')[1].split('">')[0]
		URL_list.append(URL_i)
		band_i = outlist[i].split(".stack.",1)[1].split('">',1)[0]
		bands_list.append(band_i)
	
	#Download fits images
	nURL = np.size(URL_list)
	for i in range(0,nURL):
		sub_curl = 'curl -o '+bands_list[i]+'_'+output_file+' -X POST -d "'+URL_list[i].split("?")[1]+'" http:'+URL_list[i].split("?")[0]
		print("Downloading band "+bands_list[i]+": "+bands_list[i]+'_'+output_file)
		sub_proc = subprocess.Popen([sub_curl], stdout=subprocess.PIPE,shell=True)
		(sub_out, sub_err) = sub_proc.communicate()
