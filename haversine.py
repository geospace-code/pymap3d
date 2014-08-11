from numpy import cos,arcsin,sqrt,radians,degrees
# Michael Hirsch
# from "Astronomical Algorithms" by Jean Meeus Ch. 16 p. 111 (16.5)


def angledist(r0,d0,r1,d1):
    #assumes degrees input, degrees output
    r0 = radians(r0); r1 = radians(r1)
    d0 = radians(d0); d1 = radians(d1)
    dist_rad = 2. * arcsin( sqrt( haversine(d1-d0) + cos(d0)*cos(d1) * haversine(r1-r0) ) )

    dist_deg = degrees(dist_rad)
    return dist_deg

def haversine(theta):
# http://en.wikipedia.org/wiki/Haversine and Meeus p. 111
    return (1-cos(theta)) / 2.
