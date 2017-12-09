import csv
import pdb
import pandas as pd
import numpy as np
import os
stop = pdb.set_trace

def csv_fix(file):
	with open(file, newline='', encoding='ISO-8859-1') as csvfile:
		#Open CSV file for reading
		spamreader = csv.reader(csvfile, delimiter=';', dialect='excel')
		for ii, row in enumerate(spamreader):
			#If first row read line as header for data frame
			# and initiate empty data frame
			if ii == 0:
				header = row
				df = pd.DataFrame()
			#If not first row store it as a line in the data frame but first replace line returns with symbol __N__ so that other programs can read that CSV file
			else:
				rowmod = []
				for rowi in row:
					rowmod.append(rowi.replace('\n','__N__'))
				dfi = pd.DataFrame(data=np.reshape(np.array(rowmod),(1,len(row))),columns=header)
				df = df.append(dfi)
	#Save output file
	outfile = os.path.splitext(file)[0]+'_fixed.csv'
	df.to_csv(outfile)
	print('Done !')