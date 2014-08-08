from numpy import sin, cos, degrees, radians,arcsin
import sys
sys.path.append('../astrometry')
from datetime2hourangle import datetime2sidereal

def azel2radec(az_deg,el_deg,lat_deg,lon_deg,dtime):
#from D.Vallado Fundamentals of Astrodynamics and Applications p.258-259
    az = radians(az_deg);   el = radians(el_deg)
    lat = radians(lat_deg); lon = radians(lon_deg)

    #Vallado "algorithm 28" p 268
    dec = arcsin( sin(el) * sin(lat) + cos(el) * cos(lat) * cos(az) )

    lha = arctan2( -(sin(az) * cos(el)) / cos(dec),
                   (sin(el) - sin(lat)*sin(dec)) / (cos(dec) * cos(lat)) )

    lst = datetime2sidereal(dtime,lon,ra) #lon, ra in RADIANS
    ra  = lst - lha

    return degrees(ra), degrees(dec)