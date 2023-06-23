
from jg_finderchart import *
from mocapy import MocaEngine
import numpy as np
import argparse
import os
import warnings
import sys

output_dir = os.getenv('MOCA_FINDER_CHARTS_DIR',None)

# Identify the moca_oid passed to this program
parser = argparse.ArgumentParser(description='Script to generate a finder chart for a specific moca_oid from the MOCA database.')

# Add the arguments
parser.add_argument('moca_oid', 
                    metavar='moca_oid', 
                    type=int, 
                    help='The unique object identifier in the MOCA database')

# Execute the parse_args() method
args = parser.parse_args()

if args.moca_oid is None:
    warnings.warn("A moca_oid argument is required")
    sys.exit()

# Query MOCAdb data about our target
moca = MocaEngine()
moca_oid = args.moca_oid
npos_max = 5#Max number of positions currently limiting fcharts

df = moca.query("SELECT mo.designation, ec.measurement_epoch_yr AS isort, ec.ra, ec.`dec`, ec.measurement_epoch_yr AS epoch, CONCAT_WS(' ',mission_name,data_release) AS mission FROM data_equatorial_coordinates ec JOIN moca_objects mo USING(moca_oid) WHERE moca_oid="+str(moca_oid)+" ORDER BY id LIMIT "+str(npos_max))

if len(df)==0:
    warnings.warn("No MOCAdb data was found for that moca_oid")
    sys.exit()


df['mission'] = df['mission'].replace('WISE WISE', 'WISE')
df['mission'] = df['mission'].replace('WISE AllWISE', 'AllWISE')
df['mission'] = df['mission'].replace('WISE CatWISE', 'CatWISE')
df['mission'] = df['mission'].replace('pan_starrs DR1', 'Pan-STARRS DR1')
df['mission'] = df['mission'].replace('pan_starrs DR2', 'Pan-STARRS DR2')

# Interpret the results
df = df.sort_values('isort')
df['pos'] = df['ra'].astype(str) + ' ' + df['dec'].astype(str)
df['tag'] = df['mission'] + ' (' + df['epoch'].round(1).astype(str) + ')'

# create a list to hold the values

values = df['pos'].tolist() + [None]*npos_max
tags = df['tag'].tolist() + [None]*npos_max

# create the variables
for i in range(1, npos_max+1):
    exec(f'pos{i} = values[{i-1}]')
    exec(f'tag{i} = tags[{i-1}]')

# Query Gaia DR3 positions around the target
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u

coord = SkyCoord(ra=df['ra'].mean(), dec=df['dec'].mean(), unit=(u.degree, u.degree), frame='icrs')
radius = u.Quantity(2.0/60, u.deg)
j = Gaia.cone_search_async(coord, radius)
r = j.get_results()

minsize = 0.2
all_gaia_sizes = np.log10(np.maximum(np.nan_to_num(np.array(r["parallax"])),10**minsize))
addtext = ''

finder(pos1, filename='fchart_moca_oid_'+str(moca_oid), secondary=pos2, pos3=pos3, pos4=pos4, pos5=pos5, addtext=addtext, title=df['designation'].iloc[0], override_directory=output_dir, primarypos_label=tag1, secondarypos_label=tag2, pos3_label=tag3, pos4_label=tag4, pos5_label=tag5, gray_label='All Gaia DR3', pos_list_gray_ra=np.array(r["ra"]).tolist(), pos_list_gray_dec=np.array(r["dec"]).tolist(), pos_list_gray_pmra=np.nan_to_num(np.array(r["pmra"])).tolist(), pos_list_gray_pmdec=np.nan_to_num(np.array(r["pmdec"])).tolist(),  pos_list_gray_sizes=all_gaia_sizes.tolist())
