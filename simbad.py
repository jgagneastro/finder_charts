def simbad( source_name=None, verbose=False ):
    """Gets (radeg,dedeg) for a given SIMBAD resolvable name
       
       (radedeg,dedeg) = simbad_coord( source_name )

       output = tuple of (RA,DEC)"""
    
    from urllib2 import urlopen,Request,quote
    from xml.dom.minidom import parse
    from numpy import float

    if source_name==None:
       print("Missing source name....")
       return None
    
    try:
	ra,de = source_name.split()
	ra = float(ra)
	de = float(de)
	return (float(ra),float(de))
    except:
#   cds. server is about 15% faster in response...
#   url ='http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-ox?'
    	url ='http://cds.u-strasbg.fr/cgi-bin/nph-sesame/-ox?'
    	url += quote(source_name)   
    	
    	if verbose: print(url)

    	res = urlopen(Request( url ))
    	xml = parse(res)
    	res.close()

    	if xml.getElementsByTagName('jradeg')==[]:
    	   return (None,None)
    	else:
    	   jradeg=xml.getElementsByTagName('jradeg')[0].childNodes[0].data
    	   jdedeg=xml.getElementsByTagName('jdedeg')[0].childNodes[0].data                                      
    	   return (float(jradeg),float(jdedeg))
