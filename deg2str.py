from numpy import sign

def deg2str(degrees,sp=':',ra=1,dec=0,prec=2):
  """Convert decimal degress to sexagecimal expression.
     Inseok Song (2007)

     >>> deg2str(123.4567)
     Out: '08:13:49.61'
     >>> deg2str(123.4567,prec=3)
     Out: '08:13:49.608'
     >>> deg2str(123.4567,prec=3,sp=' ')
     Out: '08 13 49.608'
     >>> deg2str(-23.4567,dec=1)
     Out: '-23:27.24.12'
     ...
  """

  if (ra & ~dec):
     if ( (degrees < 0.0) or (degrees > 360.0) ):
        str = ' 0>RA<360 '
	return(str)

     format = '%02d%s%02d%s%0'+'%1d.%1df' % (prec+3,prec)
     degrees = degrees/15.0
     rah = int(degrees)
     ram = int((degrees-rah)*60.0)
     ras = ((degrees - rah)*60 - ram)*60.0
     # str = '%02d%s%02d%s%05.2f' % (rah,sp,ram,sp,ras)
     str = format % (rah,sp,ram,sp,ras)
  elif (dec):
     if (degrees >= 0): 
        mysign='+'
     else:
        mysign='-'

     if (abs(degrees) > 90):
        str = ' |DEC|>90  '
	return(str)

     format = '%c%02.2d%s%02d%s%0'+'%1d.%1df' % (prec+3,prec)
     degrees = abs(degrees)
     ded = int(degrees)
     dem = int((degrees-ded)*60.0)
     des = ((degrees - ded)*60 - dem)*60.0
     str = format % (mysign,ded,sp,dem,sp,des)

  return(str)
