from astropy.coordinates import SkyCoord
from hips import WCSGeometry
from hips import make_sky_image

#Downloads DES DR1 fits images from HiPS
#size is in arcmin, default 10 arcmin
def query_des_fits(ra, dec, size=10.0, output_file='DES_TMP.fits', xpix=500, ypix=500, projection='AIT'):
	
	coord = SkyCoord(ra, dec, unit='deg', frame='icrs')
	fov = str(round(size/60*100)/100)+' deg'
	geometry = WCSGeometry.create(skydir=coord, width=xpix, height=ypix, fov=fov, coordsys='icrs', projection=projection)
	
	hips_survey_base = 'CDS/P/DES-DR1/'

	#Split string in a list and search for fits image URLs & band names
	bands_list = ['g','r','i','z','Y']
	for i in range(0,len(bands_list)):
		try:
			print("Downloading band "+bands_list[i]+": "+bands_list[i]+'_'+output_file)
			hips_survey = hips_survey_base+bands_list[i]
			result = make_sky_image(geometry, hips_survey, 'fits')
			result.write_image(bands_list[i]+'_'+output_file)
		except:
			print("Failed")