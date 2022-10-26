from math import inf, radians

import pymap3d.latitude as latitude
import pymap3d.rcurve as rcurve
import pytest
from pytest import approx


@pytest.mark.parametrize(
    "geodetic_lat,alt_m,geocentric_lat",
    [(0, 0, 0), (90, 0, 90), (-90, 0, -90), (45, 0, 44.80757678), (-45, 0, -44.80757678)],
)
def test_geodetic_alt_geocentric(geodetic_lat, alt_m, geocentric_lat):
    assert latitude.geod2geoc(geodetic_lat, alt_m) == approx(geocentric_lat)

    r = rcurve.geocentric_radius(geodetic_lat)
    assert latitude.geoc2geod(geocentric_lat, r) == approx(geodetic_lat)
    assert latitude.geoc2geod(geocentric_lat, 1e5 + r) == approx(
        latitude.geocentric2geodetic(geocentric_lat, 1e5 + alt_m)
    )

    assert latitude.geod2geoc(geodetic_lat, 1e5 + alt_m) == approx(
        latitude.geodetic2geocentric(geodetic_lat, 1e5 + alt_m)
    )


@pytest.mark.parametrize(
    "geodetic_lat,geocentric_lat",
    [(0, 0), (90, 90), (-90, -90), (45, 44.80757678), (-45, -44.80757678)],
)
def test_geodetic_geocentric(geodetic_lat, geocentric_lat):

    assert latitude.geodetic2geocentric(geodetic_lat, 0) == approx(geocentric_lat)
    assert latitude.geodetic2geocentric(radians(geodetic_lat), 0, deg=False) == approx(
        radians(geocentric_lat)
    )

    assert latitude.geocentric2geodetic(geocentric_lat, 0) == approx(geodetic_lat)
    assert latitude.geocentric2geodetic(radians(geocentric_lat), 0, deg=False) == approx(
        radians(geodetic_lat)
    )


def test_numpy_geodetic_geocentric():
    pytest.importorskip("numpy")
    assert latitude.geodetic2geocentric([45, 0], 0) == approx([44.80757678, 0])
    assert latitude.geocentric2geodetic([44.80757678, 0], 0) == approx([45, 0])


@pytest.mark.parametrize(
    "geodetic_lat, isometric_lat",
    [
        (0, 0),
        (90, inf),
        (-90, -inf),
        (45, 50.227466),
        (-45, -50.227466),
        (89, 271.275),
    ],
)
def test_geodetic_isometric(geodetic_lat, isometric_lat):
    isolat = latitude.geodetic2isometric(geodetic_lat)
    assert isolat == approx(isometric_lat)
    assert isinstance(isolat, float)

    assert latitude.geodetic2isometric(radians(geodetic_lat), deg=False) == approx(
        radians(isometric_lat)
    )

    assert latitude.isometric2geodetic(isometric_lat) == approx(geodetic_lat)
    assert latitude.isometric2geodetic(radians(isometric_lat), deg=False) == approx(
        radians(geodetic_lat)
    )


def test_numpy_geodetic_isometric():
    pytest.importorskip("numpy")
    assert latitude.geodetic2isometric([45, 0]) == approx([50.227466, 0])
    assert latitude.isometric2geodetic([50.227466, 0]) == approx([45, 0])


@pytest.mark.parametrize(
    "geodetic_lat,conformal_lat",
    [(0, 0), (90, 90), (-90, -90), (45, 44.80768406), (-45, -44.80768406), (89, 88.99327)],
)
def test_geodetic_conformal(geodetic_lat, conformal_lat):
    clat = latitude.geodetic2conformal(geodetic_lat)
    assert clat == approx(conformal_lat)
    assert isinstance(clat, float)

    assert latitude.geodetic2conformal(radians(geodetic_lat), deg=False) == approx(
        radians(conformal_lat)
    )

    assert latitude.conformal2geodetic(conformal_lat) == approx(geodetic_lat)
    assert latitude.conformal2geodetic(radians(conformal_lat), deg=False) == approx(
        radians(geodetic_lat)
    )


def test_numpy_geodetic_conformal():
    pytest.importorskip("numpy")
    assert latitude.geodetic2conformal([45, 0]) == approx([44.80768406, 0])
    assert latitude.conformal2geodetic([44.80768406, 0]) == approx([45, 0])


@pytest.mark.parametrize(
    "geodetic_lat,rectifying_lat",
    [(0, 0), (90, 90), (-90, -90), (45, 44.855682), (-45, -44.855682)],
)
def test_geodetic_rectifying(geodetic_lat, rectifying_lat):
    assert latitude.geodetic2rectifying(geodetic_lat) == approx(rectifying_lat)
    assert latitude.geodetic2rectifying(radians(geodetic_lat), deg=False) == approx(
        radians(rectifying_lat)
    )

    assert latitude.rectifying2geodetic(rectifying_lat) == approx(geodetic_lat)
    assert latitude.rectifying2geodetic(radians(rectifying_lat), deg=False) == approx(
        radians(geodetic_lat)
    )


def test_numpy_geodetic_rectifying():
    pytest.importorskip("numpy")
    assert latitude.geodetic2rectifying([45, 0]) == approx([44.855682, 0])
    assert latitude.rectifying2geodetic([44.855682, 0]) == approx([45, 0])


@pytest.mark.parametrize(
    "geodetic_lat,authalic_lat",
    [(0, 0), (90, 90), (-90, -90), (45, 44.87170288), (-45, -44.87170288)],
)
def test_geodetic_authalic(geodetic_lat, authalic_lat):
    assert latitude.geodetic2authalic(geodetic_lat) == approx(authalic_lat)
    assert latitude.geodetic2authalic(radians(geodetic_lat), deg=False) == approx(
        radians(authalic_lat)
    )

    assert latitude.authalic2geodetic(authalic_lat) == approx(geodetic_lat)
    assert latitude.authalic2geodetic(radians(authalic_lat), deg=False) == approx(
        radians(geodetic_lat)
    )


def test_numpy_geodetic_authalic():
    pytest.importorskip("numpy")
    assert latitude.geodetic2authalic([45, 0]) == approx([44.87170288, 0])
    assert latitude.authalic2geodetic([44.87170288, 0]) == approx([45, 0])


@pytest.mark.parametrize(
    "geodetic_lat,parametric_lat",
    [(0, 0), (90, 90), (-90, -90), (45, 44.9037878), (-45, -44.9037878)],
)
def test_geodetic_parametric(geodetic_lat, parametric_lat):
    assert latitude.geodetic2parametric(geodetic_lat) == approx(parametric_lat)
    assert latitude.geodetic2parametric(radians(geodetic_lat), deg=False) == approx(
        radians(parametric_lat)
    )

    assert latitude.parametric2geodetic(parametric_lat) == approx(geodetic_lat)
    assert latitude.parametric2geodetic(radians(parametric_lat), deg=False) == approx(
        radians(geodetic_lat)
    )


def test_numpy_geodetic_parametric():
    pytest.importorskip("numpy")
    assert latitude.geodetic2parametric([45, 0]) == approx([44.9037878, 0])
    assert latitude.parametric2geodetic([44.9037878, 0]) == approx([45, 0])


@pytest.mark.parametrize("lat", [91, -91])
def test_badvals(lat):
    # geodetic_isometric is not included on purpose
    with pytest.raises(ValueError):
        latitude.geodetic2geocentric(lat, 0)
    with pytest.raises(ValueError):
        latitude.geocentric2geodetic(lat, 0)
    with pytest.raises(ValueError):
        latitude.geodetic2conformal(lat)
    with pytest.raises(ValueError):
        latitude.conformal2geodetic(lat)
    with pytest.raises(ValueError):
        latitude.geodetic2rectifying(lat)
    with pytest.raises(ValueError):
        latitude.rectifying2geodetic(lat)
    with pytest.raises(ValueError):
        latitude.geodetic2authalic(lat)
    with pytest.raises(ValueError):
        latitude.authalic2geodetic(lat)
    with pytest.raises(ValueError):
        latitude.geodetic2parametric(lat)
    with pytest.raises(ValueError):
        latitude.parametric2geodetic(lat)
