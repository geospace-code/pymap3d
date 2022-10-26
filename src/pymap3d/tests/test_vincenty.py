import pymap3d.vincenty as vincenty
from pytest import approx


def test_track2():
    lats, lons = vincenty.track2(40, 80, 65, -148, npts=3)
    assert lats == approx([40, 69.633139886, 65])
    assert lons == approx([80, 113.06849104, -148])
