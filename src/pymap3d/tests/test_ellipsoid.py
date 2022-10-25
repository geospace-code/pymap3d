import pytest
from pytest import approx
import pymap3d as pm

xyz0 = (660e3, -4700e3, 4247e3)


@pytest.mark.parametrize(
    "model,f",
    [
        ("maupertuis", 0.005235602050865236),
        ("plessis", 0.003240020729165458),
        ("everest1830", 0.003324448922118313),
        ("everest1830m", 0.003324449295589469),
        ("everest1967", 0.003324449343845343),
        ("airy", 0.00334085067870327),
        ("bessel", 0.0033427731536659813),
        ("clarke1866", 0.0033900753039287908),
        ("clarke1878", 0.003407549790771363),
        ("clarke1860", 0.003407561308111843),
        ("helmert", 0.0033523298109184524),
        ("hayford", 0.003367003387062615),
        ("international1924", 0.003367003387062615),
        ("krassovsky1940", 0.0033523298336767685),
        ("wgs66", 0.0033528919458556804),
        ("australian", 0.003352891899858333),
        ("international1967", 0.003352896192983603),
        ("grs67", 0.0033529237272191623),
        ("sa1969", 0.003352891899858333),
        ("wgs72", 0.0033527794541680267),
        ("grs80", 0.0033528106811816882),
        ("wgs84", 0.0033528106647473664),
        ("wgs84_mean", 0.0),
        ("iers1989", 0.0033528131102879993),
        ("iers2003", 0.0033528131084554157),
        ("mercury", 0.0009014546199549272),
        ("venus", 0.0),
        ("moon", 0.0012082158679017317),
        ("mars", 0.006123875928193323),
        ("jupyter", 0.06604858798757626),
        ("io", 0.0075968738044488665),
        ("uranus", 0.022927344575296372),
        ("neptune", 0.01708124697141011),
        ("pluto", 0.0),
    ],
)
def test_reference(model, f):
    assert pm.Ellipsoid.from_name(model).flattening == approx(f)


def test_ellipsoid():

    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("maupertuis")) == approx(
        [42.123086280313906, -82.00647850636021, -13462.822154350226]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("plessis")) == approx(
        [42.008184833614905, -82.00647850636021, 1566.9219075104988]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("everest1830")) == approx(
        [42.01302648557789, -82.00647850636021, 1032.4153744896425]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("everest1830m")) == approx(
        [42.0130266467127, -82.00647850636021, 1027.7254294115853]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("everest1967")) == approx(
        [42.01302648557363, -82.00647850636021, 1033.2243733811288]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("airy")) == approx(
        [42.01397060398504, -82.00647850636021, 815.5499438015993]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("bessel")) == approx(
        [42.01407537004288, -82.00647850636021, 987.0246149983182]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("clarke1866")) == approx(
        [42.01680003414445, -82.00647850636021, 313.90267925120395]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("clarke1878")) == approx(
        [42.0177971504227, -82.00647850636021, 380.12002203958457]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("clarke1860")) == approx(
        [42.017799612218326, -82.00647850636021, 321.0980872430816]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("helmert")) == approx(
        [42.01464497456125, -82.00647850636021, 212.63680219872765]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("hayford")) == approx(
        [42.01548834310426, -82.00647850636021, 66.77070154259877]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("international1924")) == approx(
        [42.01548834310426, -82.00647850636021, 66.77070154259877]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("krassovsky1940")) == approx(
        [42.01464632634865, -82.00647850636021, 167.7043859419633]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("wgs66")) == approx(
        [42.014675415414274, -82.00647850636021, 269.1575142686737]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("australian")) == approx(
        [42.01467586302664, -82.00647850636021, 254.17989315657786]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("international1967")) == approx(
        [42.01467603307557, -82.00647850636021, 256.6883857005818]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("grs67")) == approx(
        [42.01467768000789, -82.00647850636021, 254.27066653452297]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("sa1969")) == approx(
        [42.01467586302664, -82.00647850636021, 254.17989315657786]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("wgs72")) == approx(
        [42.01466869328149, -82.00647850636021, 278.8216763935984]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("grs80")) == approx(
        [42.01467053601299, -82.00647850636021, 276.9137384511387]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("wgs84")) == approx(
        [42.01467053507479, -82.00647850636021, 276.91369158042767]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("wgs84_mean")) == approx(
        [41.823366301, -82.0064785, -2.13061272e3]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("iers1989")) == approx(
        [42.01467064467172, -82.00647850636021, 277.9191657339711]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("iers2003")) == approx(
        [42.01467066257621, -82.00647850636021, 277.320060889772]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("mercury")) == approx(
        [41.8430384333997, -82.00647850636021, 3929356.5648451606]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("venus")) == approx(
        [41.82336630167669, -82.00647850636021, 317078.15867127385]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("moon")) == approx(
        [41.842147614909734, -82.00647850636021, 4631711.995926845]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("mars")) == approx(
        [42.00945156056578, -82.00647850636021, 2981246.073616111]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("jupyter")) == approx(
        [75.3013267078341, -82.00647850636021, -61782040.202975556]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("io")) == approx(
        [41.82422244977044, -82.00647850636021, 6367054.626528843]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("uranus")) == approx(
        [47.69837228395133, -82.00647850636021, -18904824.4361074]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("neptune")) == approx(
        [45.931317431546425, -82.00647850636021, -18194050.781948525]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid.from_name("pluto")) == approx(
        [41.82336630167669, -82.00647850636021, 5180878.158671274]
    )
