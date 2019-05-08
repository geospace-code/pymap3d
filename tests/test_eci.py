import pytest
from pytest import approx
import pymap3d as pm

t0 = '2014-04-06T08:00:00'
t1 = '2013-01-15T12:00:05'
eci0 = (-3.977913815668146e6, -2.582332196263046e6, 4.250818828152067e6)


@pytest.mark.parametrize('useastropy', [True, False])
def test_eciecef(useastropy):
    ecef = pm.eci2ecef(*eci0, t1, useastropy=useastropy)
    assert ecef == approx([649012.04640917, -4697980.55129606, 4250818.82815207], rel=0.001)

    assert pm.ecef2eci(*ecef, t1, useastropy=useastropy) == approx(eci0, rel=0.001)


@pytest.mark.parametrize('useastropy', [True, False])
def test_eci_times(useastropy):
    with pytest.raises(ValueError):
        pm.eci2ecef(*eci0, [t0, t0], useastropy=useastropy)

    with pytest.raises(ValueError):
        pm.ecef2eci(*eci0, [t0, t0], useastropy=useastropy)

    x = [eci0[0]] * 2
    y = [eci0[1]] * 2
    z = [eci0[2]] * 2
    t = [t0] * 2
    assert pm.ecef2eci(*pm.eci2ecef(x, y, z, t, useastropy=useastropy),
                       t, useastropy=useastropy) == approx(eci0, rel=0.001)


if __name__ == '__main__':
    pytest.main(['-x', __file__])
