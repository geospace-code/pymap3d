import argparse

from .. import vincenty

p = argparse.ArgumentParser(
    description="Given starting latitude, longitude: find final lat,lon for distance and azimuth"
)
p.add_argument("lat", help="latitude WGS-84 [degrees]", type=float)
p.add_argument("lon", help="longitude WGS-84 [degrees]", type=float)
p.add_argument("range", help="range from start point [meters]", type=float)
p.add_argument("azimuth", help="clockwise from north: azimuth to start [degrees]", type=float)
P = p.parse_args()

lat2, lon2 = vincenty.vreckon(P.lat, P.lon, P.range, P.azimuth)

print(f"lat, lon = ({lat2.item():.4f}, {lon2.item():.4f})")
