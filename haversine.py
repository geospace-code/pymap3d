from numpy import cos
# Michael Hirsch
# based on: http://en.wikipedia.org/wiki/Haversine
def haversine(theta):
    return (1-cos(theta)) / 2.
