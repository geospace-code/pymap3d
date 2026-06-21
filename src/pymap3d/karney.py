"""
Karney's geodesic algorithms for the ellipsoid.

Pure Python implementation based on:
    C.F.F. Karney, "Algorithms for geodesics",
    J. Geodesy 87, 43-55 (2013)
    https://doi.org/10.1007/s00190-012-0578-z

Algorithm structure and coefficient tables adapted from geographiclib
(MIT license, Copyright (c) Charles Karney).

More robust than Vincenty's method, particularly for near-antipodal
points and nearly-equatorial lines.
"""

from __future__ import annotations

import math as _math

from .ellipsoid import Ellipsoid
from .mathfun import copysign, fmod, hypot, isnan, sqrt

__all__ = ["geodesic_direct", "geodesic_inverse", "geodesic_line", "geodesic_area"]

# --- Constants ---
_tiny = sqrt(2.220446049250313e-16)
_tol0 = 2.220446049250313e-16
_tol1 = 200 * _tol0
_tol2 = sqrt(_tol0)
_xthresh = 1000 * _tol2
_maxit1 = 20
_maxit2 = _maxit1 + 10 + 1

_nC1 = 6
_nC1p = 6
_nC2 = 6
_nC3 = 6
_nC4 = 6

# Outmask flags
_DISTANCE = 1
_REDUCEDLENGTH = 1 << 10
_GEODESICSCALE = 1 << 11
_OUT = _DISTANCE | _REDUCEDLENGTH | _GEODESICSCALE


# --- Low-level helpers ---

def _sum_error(u, v):
    s = u + v
    up = s - v
    vpp = s - up
    up -= u
    vpp -= v
    return s, -(up + vpp)


def _AngNormalize(x):
    y = fmod(x, 360.0)
    if y <= -180:
        y += 360.0
    elif y > 180:
        y -= 360.0
    return y


def _AngRound(x):
    z = 1 / 16.0
    y = abs(x)
    if y < z:
        y = z - y
        y = z - y
    return copysign(y, x)


def _AngDiff(x, y):
    d, t = _sum_error(-x, y)
    d = _AngNormalize(d)
    d, t2 = _sum_error(d if abs(d) != 180 else copysign(180.0, -t), t)
    return d, t2


def _sincosd(x):
    r = fmod(x, 360.0)
    q = 0 if isnan(r) else int(round(r / 90.0))
    r -= 90.0 * q
    r = _math.radians(r)
    s, c = _math.sin(r), _math.cos(r)
    q = q % 4
    if q == 1:
        s, c = c, -s
    elif q == 2:
        s, c = -s, -c
    elif q == 3:
        s, c = -c, s
    if s == 0:
        s = 0.0
    if c == 0:
        c = 0.0
    return s, c


def _atan2d(y, x):
    q = 0
    if abs(y) > abs(x):
        x, y = y, x
        q = 2
    if x < 0:
        x = -x
        q += 1
    ang = _math.degrees(_math.atan2(y, x))
    if q == 1:
        ang = copysign(180.0, y) - ang
    elif q == 2:
        ang = 90.0 - ang
    elif q == 3:
        ang = -90.0 + ang
    return ang


def _norm2(x, y):
    r = hypot(x, y)
    return x / r, y / r


def _polyval(N, p, s, x):
    y = p[s]
    for i in range(1, N + 1):
        y = y * x + p[s + i]
    return y


# --- Series evaluations ---
# Coefficient tables from geographiclib (MIT license).

_A1coeff = (1, 4, 64, 0, 256)
_C1coeff = (
    -1, 6, -16, 32,
    -9, 64, -128, 2048,
    9, -16, 768,
    3, -5, 512,
    -7, 1280,
    -7, 2048,
)
_C1pcoeff = (
    205, -432, 768, 1536,
    4005, -4736, 3840, 12288,
    -225, 116, 384,
    -7173, 2695, 7680,
    3467, 7680,
    38081, 61440,
)
_A2coeff = (-11, -28, -192, 0, 256)
_C2coeff = (
    1, 2, 16, 32,
    35, 64, 384, 2048,
    15, 80, 768,
    7, 35, 512,
    63, 1280,
    77, 2048,
)
_A3coeff = (
    -3, 128,
    -2, -3, 64,
    -1, -3, -1, 16,
    3, -1, -2, 8,
    1, -1, 2,
    1, 1,
)
_C3coeff = (
    3, 128,
    2, 5, 128,
    -1, 3, 3, 64,
    -1, 0, 1, 8,
    -1, 1, 4,
    5, 256,
    1, 3, 128,
    -3, -2, 3, 64,
    1, -3, 2, 32,
    7, 512,
    -10, 9, 384,
    5, -9, 5, 192,
    7, 512,
    -14, 7, 512,
    21, 2560,
)
_C4coeff = (
    97, 15015,
    1088, 156, 45045,
    -224, -4784, 1573, 45045,
    -10656, 14144, -4576, -858, 45045,
    64, 624, -4576, 6864, -3003, 15015,
    100, 208, 572, 3432, -12012, 30030, 45045,
    1, 9009,
    -2944, 468, 135135,
    5792, 1040, -1287, 135135,
    5952, -11648, 9152, -2574, 135135,
    -64, -624, 4576, -6864, 3003, 135135,
    8, 10725,
    1856, -936, 225225,
    -8448, 4992, -1144, 225225,
    -1440, 4160, -4576, 1716, 225225,
    -136, 63063,
    1024, -208, 105105,
    3584, -3328, 1144, 315315,
    -128, 135135,
    -2560, 832, 405405,
    128, 99099,
)


def _A1m1f(eps):
    eps2 = eps * eps
    t = _polyval(3, _A1coeff, 0, eps2) / _A1coeff[4]
    return (t + eps) / (1 - eps)


def _C1f(eps, c):
    eps2 = eps * eps
    d = eps
    o = 0
    for l in range(1, _nC1 + 1):
        m = (_nC1 - l) // 2
        c[l] = d * _polyval(m, _C1coeff, o, eps2) / _C1coeff[o + m + 1]
        o += m + 2
        d *= eps


def _C1pf(eps, c):
    eps2 = eps * eps
    d = eps
    o = 0
    for l in range(1, _nC1p + 1):
        m = (_nC1p - l) // 2
        c[l] = d * _polyval(m, _C1pcoeff, o, eps2) / _C1pcoeff[o + m + 1]
        o += m + 2
        d *= eps


def _A2m1f(eps):
    eps2 = eps * eps
    t = _polyval(3, _A2coeff, 0, eps2) / _A2coeff[4]
    return (t - eps) / (1 + eps)


def _C2f(eps, c):
    eps2 = eps * eps
    d = eps
    o = 0
    for l in range(1, _nC2 + 1):
        m = (_nC2 - l) // 2
        c[l] = d * _polyval(m, _C2coeff, o, eps2) / _C2coeff[o + m + 1]
        o += m + 2
        d *= eps


def _SinCosSeries(sinp, sinx, cosx, c, n):
    ar = 2.0 * (cosx - sinx) * (cosx + sinx)
    if n % 2 == 1:
        y0 = c[n]
        y1 = 0.0
        n -= 1
    else:
        y0 = y1 = 0.0
    while n > 0:
        y1 = ar * y0 - y1 + c[n]
        n -= 1
        y0 = ar * y1 - y0 + c[n]
        n -= 1
    if sinp:
        return 2.0 * sinx * cosx * y0
    else:
        return cosx * (y0 - y1)


def _Astroid(x, y):
    p = x * x
    q = y * y
    r = (p + q - 1.0) / 6.0
    if not (q == 0 and r <= 0):
        S = p * q / 4.0
        r2 = r * r
        r3 = r * r2
        disc = S * (S + 2.0 * r3)
        u = r
        if disc >= 0:
            T3 = S + r3
            T3 += -sqrt(disc) if T3 < 0 else sqrt(disc)
            T = T3 ** (1.0 / 3.0) if T3 >= 0 else -((-T3) ** (1.0 / 3.0))
            u += T + (r2 / T if T != 0 else 0.0)
        else:
            ang = _math.atan2(sqrt(-disc), -(S + r3))
            u += 2.0 * r * _math.cos(ang / 3.0)
        v = sqrt(u * u + q)
        uv = u + v if u >= 0 else q / (v - u)
        w = (uv - q) / (2.0 * v)
        k = uv / (sqrt(uv + w * w) + w)
    else:
        k = 0.0
    return k


# --- Geodesic class ---


class _Geodesic:

    def __init__(self, ell: Ellipsoid | None = None):
        if ell is None:
            ell = Ellipsoid.from_name("wgs84")
        self.a = ell.semimajor_axis
        self.f = ell.flattening
        self._f1 = 1.0 - self.f
        self._e2 = self.f * (2.0 - self.f)
        self._ep2 = self._e2 / (self._f1 ** 2)
        self._n = self.f / (2.0 - self.f)
        self._b = self.a * self._f1
        if self._e2 == 0:
            self._c2 = self.a ** 2
        elif self._e2 > 0:
            e = sqrt(self._e2)
            self._c2 = (self.a ** 2 + self._b ** 2 * _math.atanh(e) / e) / 2.0
        else:
            e = sqrt(-self._e2)
            self._c2 = (self.a ** 2 + self._b ** 2 * _math.atan(e) / e) / 2.0
        self._etol2 = 0.1 * _tol2 / sqrt(
            max(0.001, abs(self.f)) * min(1.0, 1.0 - self.f / 2.0) / 2.0)

        # Precompute A3, C3, C4 coefficient arrays (polynomials in n evaluated)
        self._A3x = self._compute_A3x()
        self._C3x = self._compute_C3x()
        self._C4x = self._compute_C4x()

    def _compute_A3x(self):
        _A3x = [0.0] * _nC3
        o = 0
        k = 0
        for j in range(_nC3 - 1, -1, -1):
            m = min(_nC3 - j - 1, j)
            _A3x[k] = _polyval(m, _A3coeff, o, self._n) / _A3coeff[o + m + 1]
            o += m + 2
            k += 1
        return _A3x

    def _compute_C3x(self):
        _C3x = [0.0] * ((_nC3 * (_nC3 - 1)) // 2)
        o = 0
        k = 0
        for l in range(1, _nC3):
            for j in range(_nC3 - 1, l - 1, -1):
                m = min(_nC3 - j - 1, j)
                _C3x[k] = _polyval(m, _C3coeff, o, self._n) / _C3coeff[o + m + 1]
                o += m + 2
                k += 1
        return _C3x

    def _compute_C4x(self):
        _C4x = [0.0] * ((_nC4 * (_nC4 + 1)) // 2)
        o = 0
        k = 0
        for l in range(_nC4):
            for j in range(_nC4 - 1, l - 1, -1):
                m = _nC4 - j - 1
                _C4x[k] = _polyval(m, _C4coeff, o, self._n) / _C4coeff[o + m + 1]
                o += m + 2
                k += 1
        return _C4x

    def _A3f(self, eps):
        return _polyval(_nC3 - 1, self._A3x, 0, eps)

    def _C3f(self, eps, c):
        mult = 1.0
        o = 0
        for l in range(1, _nC3):
            mult *= eps
            m = _nC3 - 1 - l
            c[l] = mult * _polyval(m, self._C3x, o, eps)
            o += m + 1

    def _C4f(self, eps, c):
        mult = 1.0
        o = 0
        for l in range(_nC4):
            m = _nC4 - l - 1
            c[l] = mult * _polyval(m, self._C4x, o, eps)
            o += m + 1
            mult *= eps

    def _Lengths(self, eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                 cbet1, cbet2, outmask, C1a, C2a):
        s12b = m12b = m0 = M12 = M21 = _math.nan
        A1 = _A1m1f(eps)
        _C1f(eps, C1a)
        if outmask & (_REDUCEDLENGTH | _GEODESICSCALE):
            A2 = _A2m1f(eps)
            _C2f(eps, C2a)
            m0x = A1 - A2
            A2 = 1.0 + A2
        A1 = 1.0 + A1

        if outmask & _DISTANCE:
            B1 = (_SinCosSeries(True, ssig2, csig2, C1a, _nC1) -
                  _SinCosSeries(True, ssig1, csig1, C1a, _nC1))
            s12b = A1 * (sig12 + B1)
            if outmask & (_REDUCEDLENGTH | _GEODESICSCALE):
                B2 = (_SinCosSeries(True, ssig2, csig2, C2a, _nC2) -
                      _SinCosSeries(True, ssig1, csig1, C2a, _nC2))
                J12 = m0x * sig12 + (A1 * B1 - A2 * B2)
        elif outmask & (_REDUCEDLENGTH | _GEODESICSCALE):
            for l in range(1, _nC1 + 1):
                C2a[l] = A1 * C1a[l] - A2 * C2a[l]
            J12 = m0x * sig12 + (
                _SinCosSeries(True, ssig2, csig2, C2a, _nC1) -
                _SinCosSeries(True, ssig1, csig1, C2a, _nC1))

        if outmask & _REDUCEDLENGTH:
            m0 = m0x
            m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * J12

        if outmask & _GEODESICSCALE:
            csig12 = csig1 * csig2 + ssig1 * ssig2
            t = self._ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2)
            M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1
            M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2

        return s12b, m12b, m0, M12, M21

    def _InverseStart(self, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                      lam12, slam12, clam12, C1a, C2a):
        sig12 = -1.0
        salp1 = 2.0
        calp1 = salp2 = calp2 = dnm = _math.nan

        sbet12 = sbet2 * cbet1 - cbet2 * sbet1
        cbet12 = cbet2 * cbet1 + sbet2 * sbet1
        sbet12a = sbet2 * cbet1 + cbet2 * sbet1

        shortline = cbet12 >= 0 and sbet12 < 0.5 and cbet2 * lam12 < 0.5
        if shortline:
            sbetm2 = (sbet1 + sbet2) ** 2
            sbetm2 /= sbetm2 + (cbet1 + cbet2) ** 2
            dnm = sqrt(1.0 + self._ep2 * sbetm2)
            omg12 = lam12 / (self._f1 * dnm)
            somg12 = _math.sin(omg12)
            comg12 = _math.cos(omg12)
        else:
            somg12 = slam12
            comg12 = clam12

        salp1 = cbet2 * somg12
        calp1 = (
            sbet12 + cbet2 * sbet1 * somg12 ** 2 / (1.0 + comg12)
            if comg12 >= 0
            else sbet12a - cbet2 * sbet1 * somg12 ** 2 / (1.0 - comg12))

        ssig12 = hypot(salp1, calp1)
        csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12

        if shortline and ssig12 * self._ep2 / sqrt(1.0 + self._ep2) < 0.01:
            # Short line case
            pass

        if not (csig12 >= 0 or ssig12 >= 6 * abs(self.f) * _math.pi * cbet1 ** 2):
            if self.f >= 0:
                k2 = sbet1 ** 2 * self._ep2
                eps = k2 / (2.0 * (1.0 + sqrt(1.0 + k2)) + k2)
                lamscale = self.f * cbet1 * self._A3f(eps) * _math.pi
                betscale = lamscale * cbet1
                x = (lam12 - _math.pi) / lamscale
                y = sbet12a / betscale
            else:
                cbet12a = cbet2 * cbet1 - sbet2 * sbet1
                bet12a = _math.atan2(sbet12a, cbet12a)
                _, m12b, m0, _, _ = self._Lengths(
                    self._n, _math.pi + bet12a,
                    sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
                    cbet1, cbet2, _REDUCEDLENGTH, C1a, C2a)
                x = -1.0 + m12b / (cbet1 * cbet2 * m0 * _math.pi)
                betscale = (sbet12a / x if x < -0.01
                            else -self.f * cbet1 ** 2 * _math.pi)
                lamscale = betscale / cbet1
                y = (lam12 - _math.pi) / lamscale

            if y < -_tol1 and x > -1.0 - _xthresh:
                if self.f >= 0:
                    salp1 = min(1.0, -x)
                    calp1 = -sqrt(1.0 - salp1 ** 2)
                else:
                    calp1 = max(-1.0 if x > -_tol1 else -1.0, x)
                    salp1 = sqrt(1.0 - calp1 ** 2)
            else:
                k = _Astroid(x, y)
                omg12a = lamscale * (
                    -x * k / (1.0 + k) if self.f >= 0
                    else -y * (1.0 + k) / k)
                somg12 = _math.sin(omg12a)
                comg12 = -_math.cos(omg12a)
                salp1 = cbet2 * somg12
                calp1 = sbet12a - cbet2 * sbet1 * somg12 ** 2 / (1.0 - comg12)

        if salp1 > 0 or isnan(salp1):
            salp1, calp1 = _norm2(salp1, calp1)
        else:
            salp1 = 1.0
            calp1 = 0.0
        return sig12, salp1, calp1, salp2, calp2, dnm

    def _Lambda12(self, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                  salp1, calp1, slam120, clam120, diffp, C1a, C2a, C3a):
        if sbet1 == 0 and calp1 == 0:
            calp1 = -_tiny

        salp0 = salp1 * cbet1
        calp0 = hypot(calp1, salp1 * sbet1)

        ssig1 = sbet1
        somg1 = salp0 * sbet1
        csig1 = comg1 = calp1 * cbet1
        ssig1, csig1 = _norm2(ssig1, csig1)

        salp2 = salp0 / cbet2 if cbet2 != cbet1 else salp1
        calp2 = (
            sqrt((calp1 * cbet1) ** 2 +
                 ((cbet2 - cbet1) * (cbet1 + cbet2) if cbet2 != cbet1
                  else (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2
            if cbet2 != cbet1 or abs(sbet2) != -sbet1
            else abs(calp1))

        ssig2 = sbet2
        somg2 = salp0 * sbet2
        csig2 = comg2 = calp2 * cbet2
        ssig2, csig2 = _norm2(ssig2, csig2)

        sig12 = _math.atan2(
            max(0.0, csig1 * ssig2 - ssig1 * csig2),
            csig1 * csig2 + ssig1 * ssig2)

        somg12 = max(0.0, comg1 * somg2 - somg1 * comg2)
        comg12 = comg1 * comg2 + somg1 * somg2
        eta = _math.atan2(
            somg12 * clam120 - comg12 * slam120,
            comg12 * clam120 + somg12 * slam120)

        k2 = calp0 ** 2 * self._ep2
        eps = k2 / (2.0 * (1.0 + sqrt(1.0 + k2)) + k2)
        self._C3f(eps, C3a)
        B312 = (_SinCosSeries(True, ssig2, csig2, C3a, _nC3 - 1) -
                _SinCosSeries(True, ssig1, csig1, C3a, _nC3 - 1))
        A3 = self._A3f(eps)
        domg12 = -self.f * A3 * salp0 * (sig12 + B312)
        lam12 = eta + domg12

        if diffp:
            if calp2 == 0:
                dlam12 = -2.0 * self._f1 * dn1 / sbet1
            else:
                _, dlam12, _, _, _ = self._Lengths(
                    eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                    cbet1, cbet2, _REDUCEDLENGTH, C1a, C2a)
                dlam12 *= self._f1 / (calp2 * cbet2)
        else:
            dlam12 = _math.nan

        return (lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2,
                eps, somg12, comg12, dlam12)

    def _Inverse(self, lat1, lon1, lat2, lon2):
        lon12, lon12s = _AngDiff(lon1, lon2)
        lonsign = 1.0 if lon12 >= 0 else -1.0
        lon12 = lonsign * _AngRound(lon12)
        lon12s = _AngRound((180.0 - lon12) - lonsign * lon12s)
        lam12 = _math.radians(lon12)
        if lon12 > 90:
            slam12, clam12 = _sincosd(lon12s)
            clam12 = -clam12
        else:
            slam12, clam12 = _sincosd(lon12)

        lat1 = _AngRound(_AngNormalize(lat1))
        lat2 = _AngRound(lat2)

        swapp = -1.0 if abs(lat1) < abs(lat2) else 1.0
        if swapp < 0:
            lonsign *= -1.0
            lat1, lat2 = lat2, lat1

        latsign = 1.0 if lat1 < 0 else -1.0
        lat1 *= latsign
        lat2 *= latsign

        sbet1, cbet1 = _sincosd(lat1)
        sbet1 *= self._f1
        sbet1, cbet1 = _norm2(sbet1, cbet1)
        cbet1 = max(_tiny, cbet1)

        sbet2, cbet2 = _sincosd(lat2)
        sbet2 *= self._f1
        sbet2, cbet2 = _norm2(sbet2, cbet2)
        cbet2 = max(_tiny, cbet2)

        if cbet1 < -sbet1:
            if cbet2 == cbet1:
                sbet2 = sbet1 if sbet2 < 0 else -sbet1

        dn1 = sqrt(1.0 + self._ep2 * sbet1 ** 2)
        dn2 = sqrt(1.0 + self._ep2 * sbet2 ** 2)

        meridian = lat1 == -90 or slam12 == 0

        C1a = [0.0] * (_nC1 + 1)
        C2a = [0.0] * (_nC2 + 1)
        C3a = [0.0] * _nC3

        s12x = m12x = a12 = _math.nan
        M12 = M21 = _math.nan

        if meridian:
            calp1 = clam12
            salp1 = slam12
            calp2 = 1.0
            salp2 = 0.0

            ssig1 = sbet1
            csig1 = calp1 * cbet1
            ssig2 = sbet2
            csig2 = calp2 * cbet2

            sig12 = _math.atan2(
                max(0.0, csig1 * ssig2 - ssig1 * csig2),
                csig1 * csig2 + ssig1 * ssig2)

            s12x, m12x, _, M12, M21 = self._Lengths(
                self._n, sig12,
                ssig1, csig1, dn1, ssig2, csig2, dn2,
                cbet1, cbet2, _OUT, C1a, C2a)

            if sig12 < 1 or m12x >= 0:
                if sig12 < 3 * _tiny:
                    sig12 = m12x = 0.0
                m12x *= self._b
                s12x *= self._b
                a12 = _math.degrees(sig12)
            else:
                meridian = False

        # Equatorial case: both points on equator, geodesic runs along equator
        if (not meridian and
            sbet1 == 0 and
            (self.f <= 0 or lon12s >= self.f * 180)):
            calp1 = calp2 = 0.0
            salp1 = salp2 = 1.0
            s12x = self.a * lam12
            sig12 = lam12 / self._f1
            m12x = self._b * _math.sin(sig12)
            M12 = M21 = _math.cos(sig12)
            a12 = lon12 / self._f1

        elif not meridian:
            sig12, salp1, calp1, salp2, calp2, dnm = self._InverseStart(
                sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                lam12, slam12, clam12, C1a, C2a)

            if sig12 >= 0:
                B1 = _SinCosSeries(True, _math.sin(sig12), _math.cos(sig12), C1a, _nC1)
                s12x = self._b * (1 + _A1m1f(self._n)) * (sig12 + B1)
                m12x = 0.0
                M12 = M21 = 1.0
                a12 = _math.degrees(sig12)
            else:
                tripb = False
                salp1a = _tiny
                calp1a = 1.0
                salp1b = _tiny
                calp1b = -1.0

                for numit in range(_maxit2):
                    (v, salp2, calp2, sig12,
                     ssig1, csig1, ssig2, csig2,
                     eps, domg12, comg12, dlam12) = self._Lambda12(
                        sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                        salp1, calp1, slam12, clam12,
                        numit < _maxit1, C1a, C2a, C3a)

                    if tripb or not (abs(v) >= (8 if numit < _maxit1 else 1) * _tol0):
                        break

                    if v > 0 and (numit < _maxit1 or calp1 / salp1 > calp1b / salp1b):
                        salp1b = salp1
                        calp1b = calp1
                    elif v < 0 and (numit < _maxit1 or calp1 / salp1 < calp1a / salp1a):
                        salp1a = salp1
                        calp1a = calp1

                    if numit < _maxit1 and dlam12 != 0:
                        dalp1 = -v / dlam12
                        sdalp1 = _math.sin(dalp1)
                        cdalp1 = _math.cos(dalp1)
                        nsalp1 = salp1 * cdalp1 + calp1 * sdalp1
                        if nsalp1 > 0:
                            calp1 = calp1 * cdalp1 - salp1 * sdalp1
                            salp1 = nsalp1
                            salp1, calp1 = _norm2(salp1, calp1)
                            tripb = abs(v) <= 16 * _tol0
                            continue

                    salp1 = (salp1a + salp1b) / 2.0
                    calp1 = (calp1a + calp1b) / 2.0
                    salp1, calp1 = _norm2(salp1, calp1)
                    tripb = False
                    if abs(salp1a - salp1b) < _tol0 and abs(calp1a - calp1b) < _tol0:
                        break

                s12x, m12x, _, M12, M21 = self._Lengths(
                    eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                    cbet1, cbet2, _OUT, C1a, C2a)
                m12x *= self._b
                s12x *= self._b
                a12 = _math.degrees(sig12)

        if swapp < 0:
            salp1, salp2 = salp2, salp1
            calp1, calp2 = calp2, calp1
            M12, M21 = M21, M12

        salp1 *= swapp * lonsign
        calp1 *= swapp * latsign
        salp2 *= swapp * lonsign
        calp2 *= swapp * latsign

        azi1 = _atan2d(salp1, calp1)
        azi2 = _atan2d(salp2, calp2)
        return s12x, azi1, azi2

    def _Direct(self, lat1, lon1, azi1, s12):
        salp1, calp1 = _sincosd(_AngNormalize(azi1))
        sbet1, cbet1 = _sincosd(_AngRound(lat1))
        sbet1 *= self._f1
        sbet1, cbet1 = _norm2(sbet1, cbet1)
        cbet1 = max(_tiny, cbet1)

        dn1 = sqrt(1.0 + self._ep2 * sbet1 ** 2)
        salp0 = salp1 * cbet1
        calp0 = hypot(calp1, salp1 * sbet1)

        ssig1 = sbet1
        somg1 = salp0 * sbet1
        csig1 = comg1 = calp1 * cbet1
        # Handle equatorial case: avoid norm2(0,0)
        if ssig1 != 0 or csig1 != 0:
            ssig1, csig1 = _norm2(ssig1, csig1)
        else:
            ssig1 = 0.0
            csig1 = 1.0

        k2 = calp0 ** 2 * self._ep2
        eps = k2 / (2.0 * (1.0 + sqrt(1.0 + k2)) + k2)

        A1 = _A1m1f(eps)
        C1a = [0.0] * (_nC1 + 1)
        _C1f(eps, C1a)
        C1pa = [0.0] * (_nC1p + 1)
        _C1pf(eps, C1pa)

        A1p1 = 1.0 + A1
        B11 = _SinCosSeries(True, ssig1, csig1, C1a, _nC1)
        s_B11 = _math.sin(B11)
        c_B11 = _math.cos(B11)
        stau1 = ssig1 * c_B11 + csig1 * s_B11
        ctau1 = csig1 * c_B11 - ssig1 * s_B11

        tau12 = s12 / (self._b * A1p1)
        s_tau12 = _math.sin(tau12)
        c_tau12 = _math.cos(tau12)

        stau2 = stau1 * c_tau12 + ctau1 * s_tau12
        ctau2 = ctau1 * c_tau12 - stau1 * s_tau12

        B12 = -_SinCosSeries(True, stau2, ctau2, C1pa, _nC1p)
        sig12 = tau12 - (B12 - B11)
        ssig12 = _math.sin(sig12)
        csig12 = _math.cos(sig12)

        if abs(self.f) > 0.01:
            ssig2 = ssig1 * csig12 + csig1 * ssig12
            csig2 = csig1 * csig12 - ssig1 * ssig12
            B12 = _SinCosSeries(True, ssig2, csig2, C1a, _nC1)
            serr = (1.0 + A1) * (sig12 + (B12 - B11)) - s12 / self._b
            sig12 -= serr / sqrt(1.0 + k2 * ssig2 ** 2)
            ssig12 = _math.sin(sig12)
            csig12 = _math.cos(sig12)

        ssig2 = ssig1 * csig12 + csig1 * ssig12
        csig2 = csig1 * csig12 - ssig1 * ssig12
        dn2 = sqrt(1.0 + k2 * ssig2 ** 2)

        sbet2 = calp0 * ssig2
        cbet2 = hypot(salp0, calp0 * csig2)
        if cbet2 == 0:
            cbet2 = csig2 = _tiny

        salp2 = salp0
        calp2 = calp0 * csig2

        somg2 = salp0 * ssig2
        comg2 = csig2

        C3a = [0.0] * _nC3
        self._C3f(eps, C3a)
        A3c = self._A3f(eps)
        B31 = _SinCosSeries(True, ssig1, csig1, C3a, _nC3 - 1)
        B32 = _SinCosSeries(True, ssig2, csig2, C3a, _nC3 - 1)

        # Use LONG_UNROLL method for robustness (handles equatorial case)
        E = copysign(1.0, salp0)
        omg12 = E * (sig12
                     - (_math.atan2(ssig2, csig2) - _math.atan2(ssig1, csig1))
                     + (_math.atan2(E * somg2, comg2) - _math.atan2(E * somg1, comg1)))
        lam12 = omg12 - self.f * A3c * salp0 * (sig12 + (B32 - B31))

        lon12 = _math.degrees(lam12)
        lat2 = _atan2d(sbet2, self._f1 * cbet2)
        lon2 = _AngNormalize(lon1 + lon12)
        azi2 = _atan2d(salp2, calp2)
        return lat2, lon2, azi2


# --- Cache ---

_cache = {}


def _get_geodesic(ell: Ellipsoid | None) -> _Geodesic:
    key = None if ell is None else (ell.semimajor_axis, ell.flattening)
    if key not in _cache:
        _cache[key] = _Geodesic(ell)
    return _cache[key]


# --- Public API ---


def geodesic_inverse(
    lat1: float, lon1: float, lat2: float, lon2: float,
    ell: Ellipsoid | None = None, deg: bool = True,
) -> tuple:
    """
    Solve the inverse geodesic problem using Karney's algorithm.

    Given two points, find the distance and azimuths of the geodesic
    connecting them.

    Parameters
    ----------
    lat1 : float
        Geodetic latitude of first point.
    lon1 : float
        Geodetic longitude of first point.
    lat2 : float
        Geodetic latitude of second point.
    lon2 : float
        Geodetic longitude of second point.
    ell : Ellipsoid, optional
        Reference ellipsoid (default WGS84).
    deg : bool, optional
        If True (default), angles in degrees; if False, radians.

    Returns
    -------
    dist : float
        Geodesic distance in meters.
    azi1 : float
        Forward azimuth at first point (clockwise from north).
    azi2 : float
        Back-azimuth at second point (clockwise from north).
    """
    g = _get_geodesic(ell)
    if not deg:
        lat1, lon1 = _math.degrees(lat1), _math.degrees(lon1)
        lat2, lon2 = _math.degrees(lat2), _math.degrees(lon2)
    dist, azi1, azi2 = g._Inverse(lat1, lon1, lat2, lon2)
    if not deg:
        azi1, azi2 = _math.radians(azi1), _math.radians(azi2)
    return dist, azi1, azi2


def geodesic_direct(
    lat1: float, lon1: float, azi1: float, dist: float,
    ell: Ellipsoid | None = None, deg: bool = True,
) -> tuple:
    """
    Solve the direct geodesic problem using Karney's algorithm.

    Given a starting point, azimuth, and distance, find the endpoint
    and back-azimuth.

    Parameters
    ----------
    lat1 : float
        Geodetic latitude of starting point.
    lon1 : float
        Geodetic longitude of starting point.
    azi1 : float
        Forward azimuth at starting point (clockwise from north).
    dist : float
        Distance along geodesic in meters.
    ell : Ellipsoid, optional
        Reference ellipsoid (default WGS84).
    deg : bool, optional
        If True (default), angles in degrees; if False, radians.

    Returns
    -------
    lat2 : float
        Geodetic latitude of endpoint.
    lon2 : float
        Geodetic longitude of endpoint.
    azi2 : float
        Back-azimuth at endpoint (clockwise from north).
    """
    g = _get_geodesic(ell)
    if not deg:
        lat1, lon1 = _math.degrees(lat1), _math.degrees(lon1)
        azi1 = _math.degrees(azi1)
    lat2, lon2, azi2 = g._Direct(lat1, lon1, azi1, dist)
    if not deg:
        lat2, lon2, azi2 = _math.radians(lat2), _math.radians(lon2), _math.radians(azi2)
    return lat2, lon2, azi2


def geodesic_line(
    lat1: float, lon1: float, lat2: float, lon2: float,
    ell: Ellipsoid | None = None, npts: int = 100, deg: bool = True,
) -> tuple:
    """
    Compute intermediate points along a geodesic.

    Parameters
    ----------
    lat1, lon1 : float
        First endpoint.
    lat2, lon2 : float
        Second endpoint.
    ell : Ellipsoid, optional
        Reference ellipsoid (default WGS84).
    npts : int, optional
        Number of points (default 100).
    deg : bool, optional
        If True (default), angles in degrees; if False, radians.

    Returns
    -------
    lats : list of float
        Latitudes of intermediate points.
    lons : list of float
        Longitudes of intermediate points.
    """
    dist, azi, _ = geodesic_inverse(lat1, lon1, lat2, lon2, ell=ell, deg=deg)
    lats = []
    lons = []
    for i in range(npts):
        frac = i / (npts - 1) if npts > 1 else 0.0
        lat, lon, _ = geodesic_direct(lat1, lon1, azi, dist * frac, ell=ell, deg=deg)
        lats.append(lat)
        lons.append(lon)
    return lats, lons


def geodesic_area(
    lats, lons, ell: Ellipsoid | None = None, deg: bool = True,
) -> tuple:
    """
    Compute the area and perimeter of a geodesic polygon.

    Counterclockwise vertex ordering gives positive area.

    Parameters
    ----------
    lats : array-like
        Latitudes of polygon vertices.
    lons : array-like
        Longitudes of polygon vertices.
    ell : Ellipsoid, optional
        Reference ellipsoid (default WGS84).
    deg : bool, optional
        If True (default), angles in degrees; if False, radians.

    Returns
    -------
    area : float
        Area in square meters (signed: CCW = positive).
    perimeter : float
        Perimeter in meters.
    """
    n = len(lats)
    if n < 3:
        return 0.0, 0.0

    g = _get_geodesic(ell)
    perimeter = 0.0
    area_sum = 0.0
    azi2_prev = None
    azi1_first = None

    for i in range(n):
        j = (i + 1) % n
        dist_ij, azi1_ij, azi2_ij = geodesic_inverse(
            lats[i], lons[i], lats[j], lons[j], ell=ell, deg=deg)
        perimeter += dist_ij
        if i == 0:
            azi1_first = azi1_ij
        else:
            diff, _ = _AngDiff(azi2_prev, azi1_ij)
            area_sum += diff
        azi2_prev = azi2_ij

    diff, _ = _AngDiff(azi2_prev, azi1_first)
    area_sum += diff

    area = g._c2 * _math.radians(area_sum)
    area0 = 4.0 * _math.pi * g._c2
    if abs(area) > area0 / 2:
        area -= copysign(area0, area)
    return area, perimeter
