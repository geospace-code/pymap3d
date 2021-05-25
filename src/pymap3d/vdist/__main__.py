import argparse

from .. import vincenty

p = argparse.ArgumentParser(description="Vincenty algorithms")
p.add_argument("lat1", help="latitude1 WGS-84 [degrees]", type=float)
p.add_argument("lon1", help="longitude1 WGS-84 [degrees]", type=float)
p.add_argument("lat2", help="latitude2 WGS-84 [degrees]", type=float)
p.add_argument("lon2", help="longitude2 WGS-84 [degrees]", type=float)
P = p.parse_args()

dist_m = vincenty.vdist(P.lat1, P.lon1, P.lat2, P.lon2)

print("{:.3f} meters {:.3f} deg azimuth".format(*dist_m))
