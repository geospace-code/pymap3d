from numpy import cos,arcsin,sqrt
# Michael Hirsch
# from "Astronomical Algorithms" by Jean Meeus Ch. 16 p. 111 (16.5)


def angledist(r0,d0,r1,d1):
    dist_rad = 2.*arcsin( sqrt( haversine(d1-d0) + cos(d0)*cos(d1) * haversine(r1-r0) ) )
    return dist_rad

def haversine(theta):
# http://en.wikipedia.org/wiki/Haversine and Meeus p. 111
    return (1-cos(theta)) / 2.
