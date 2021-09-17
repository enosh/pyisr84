"""This code reimplements Joseph Gray's C++ library for conversion between ITM and WGS84 in Python for convinent scripting use.

LICENSE
=======
This code is free software; you can redistribute it and/or modify it at your will.
It is my hope however that if you improve it in any way you will find a way to share it too.

This program is distributed AS-IS in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


Israel Local Grids <==> WGS84 conversions
=========================================
The Israel New Grid (ITM) is a Transverse Mercator projection of the GRS80 ellipsoid.
The Israel Old Grid (ICS) is a Cassini-Soldner projection of the modified Clark 1880 ellipsoid.

To convert from a local grid to WGS84 you first do a "UTM to Lat/Lon" conversion using the
known formulas but with the local grid data (Central Meridian, Scale Factor and False
Easting and Northing). This results in Lat/Long in the local ellipsoid coordinate system.
Afterwards you do a Molodensky transformation from this ellipsoid to WGS84.

To convert from WGS84 to a local grid you first do a Molodensky transformation from WGS84
to the local ellipsoid, after which you do a Lat/Lon to UTM conversion, again with the data of
the local grid instead of the UTM data.

The UTM to Lat/Lon and Lat/Lon to UTM conversion formulas were taken as-is from the
excellent article by Prof. Steven Dutch of the University of Wisconsin at Green Bay:
http://www.uwgb.edu/dutchs/UsefulData/UTMFormulas.htm

The [abridged] Molodensky transformations were taken from
http://home.hiwaay.net/~taylorc/bookshelf/math-science/geodesy/datum/transform/molodensky/
and can be found in many sources on the net.

Additional sources:
===================
1. dX,dY,dZ values:  http://www.geo.hunter.cuny.edu/gis/docs/geographic_transformations.pdf

2. ITM data:  http://www.mapi.gov.il/geodesy/itm_ftp.txt
    for the meridional arc false northing, the value is given at
    http://www.mapi.gov.il/reg_inst/dir2b.doc
    (this doc also gives a different formula for Lat/lon -> ITM, but not the reverse)

3. ICS data:  http://www.mapi.gov.il/geodesy/ics_ftp.txt
   for the meridional arc false northing, the value is given at several places as the
   correction value for Garmin GPS sets, the origin is unknown.
   e.g. http://www.idobartana.com/etrexkb/etrexisr.htm

Notes:
======
1. The conversions between ICS and ITM are
        ITM Lat = ICS Lat - 500000
        ITM Lon = ICS Lon + 50000
    e.g. ITM 678000,230000 <--> ICS 1178000 180000

    Since the formulas for ITM->WGS84 and ICS->WGS84 are different, the results will differ.
    For the above coordinates we get the following results (WGS84)
        ITM->WGS84 32.11'43.945" 35.18'58.782"
        ICS->WGS84 32.11'43.873" 35.18'58.200"
    Difference         ~3m           ~15m

2. If you have, or have seen, formulas that contain the term Sin(1"), I recommend you read
   Prof. Dutch's enlightening explanation about it in his link above.
"""

from typing import Tuple
from math import pi, sqrt, sin, cos


class Datum():
    class WGS84:
        a = 6378137.0
        b = 6356752.3142
        f = 0.00335281066474748     # f = 1/298.257223563
        esq = 0.006694380004260807
        e = 0.0818191909289062
        # deltas to WGS84
        dX = 0
        dY = 0
        dZ = 0

    class GRS80:
        a = 6378137.0
        b = 6356752.3141
        f = 0.0033528106811823      # f = 1/298.257222101
        esq = 0.00669438002290272
        e = 0.0818191910428276
        # deltas to WGS84
        dX = -48
        dY = 55
        dZ = 52

    class CLARK80M:
        a = 6378300.789
        b = 6356566.4116309
        f = 0.003407549767264       # f = 1/293.466
        esq = 0.006803488139112318
        e = 0.08248325975076590
        # deltas to WGS84
        dX = -235
        dY = -85
        dZ = 264


class Grid():
    class ICS:
        lon0 = 0.6145667421719      # central meridian in radians of 35.12'43.490"
        lat0 = 0.55386447682762762  # lat0 = central latitude in radians of 31.44'02.749"
        k0 = 1.00000                # k0 = scale factor
        false_e = 170251.555        # false_easting
        false_n = 2385259.0

    class ITM:
        lon0 = 0.61443473225468920  # central meridian in radians 35.12'16.261"
        lat0 = 0.55386965463774187  # central latitude in radians 31.44'03.817"
        k0 = 1.0000067              # scale factor
        false_e = 219529.584
        false_n = 2885516.9488      # 3512424.3388-626907.390 MAPI says the false northing is 626907.390, and in another place that the meridional arc at the central latitude is 3512424.3388


def itm_to_wgs84(northing: int, easting: int) -> Tuple[float, float]:
    """Israel New Grid (ITM) to WGS84 conversion"""
    lat80, lon80 = grid_to_latlon(northing, easting, Grid.ITM, Datum.GRS80)

    lat84, lon84 = Molodensky(lat80, lon80, Datum.GRS80, Datum.WGS84)

    return (lat84 * 180 / pi, lon84 * 180 / pi)


def wgs84_to_itm(lat: float, lon: float) -> Tuple[int, int]:
    """WGS84 to Israel New Grid (ITM) conversion"""
    latr = lat * pi / 180
    lonr = lon * pi / 180

    lat80, lon80 = Molodensky(latr, lonr, Datum.WGS84, Datum.GRS80)

    return latlon_to_grid(lat80, lon80, Datum.GRS80, Grid.ITM)


def ics_to_wgs84(northing: int, easting: int) -> Tuple[float, float]:
    """Israel Old Grid (ICS) to WGS84 conversion"""
    lat80, lon80 = grid_to_latlon(northing, easting, Grid.ICS, Datum.CLARK80M)

    lat84, lon84 = Molodensky(lat80, lon80, Datum.CLARK80M, Datum.WGS84)

    return (lat84 * 180 / pi, lon84 * 180 / pi)


def wgs84_to_ics(lat: float, lon: float) -> Tuple[int, int]:
    """WGS84 to Israel Old Grid (ICS) conversion"""
    latr = lat * pi / 180
    lonr = lon * pi / 180

    lat80, lon80 = Molodensky(latr, lonr, Datum.WGS84, Datum.CLARK80M)

    return latlon_to_grid(lat80, lon80, Datum.CLARK80M, Grid.ICS)


def grid_to_latlon(northing: int, easting: int, from_grid: Grid, to_datum: Datum) -> Tuple[float, float]:
    """Local Grid to Lat/Lon conversion"""
    y = northing + from_grid.false_n
    x = easting - from_grid.false_e
    M = y / from_grid.k0

    a = to_datum.a
    b = to_datum.b
    e = to_datum.e
    esq = to_datum.esq

    mu = M / (a * (1 - e * e / 4 - 3 * (e**4) / 64 - 5 * (e**6) / 256))

    ee = sqrt(1 - esq)
    e1 = (1 - ee) / (1 + ee)
    j1 = 3 * e1 / 2 - 27 * e1 * e1 * e1 / 32
    j2 = 21 * e1 * e1 / 16 - 55 * e1 * e1 * e1 * e1 / 32
    j3 = 151 * e1 * e1 * e1 / 96
    j4 = 1097 * e1 * e1 * e1 * e1 / 512

    # Footprint Latitude
    fp = mu + j1 * sin(2 * mu) + j2 * sin(4 * mu) + j3 * sin(6 * mu) + j4 * sin(8 * mu)

    sinfp = sin(fp)
    cosfp = cos(fp)
    tanfp = sinfp / cosfp
    eg = e * a / b
    eg2 = eg*eg
    C1 = eg2*cosfp*cosfp
    T1 = tanfp*tanfp
    R1 = a * (1 - e*e) / (1 - (e*sinfp) * (e*sinfp))**1.5
    N1 = a / sqrt(1 - (e*sinfp) * (e*sinfp))
    D = x / (N1 * from_grid.k0)

    Q1 = N1 * tanfp / R1
    Q2 = D * D / 2
    Q3 = (5 + 3*T1 + 10 * C1 - 4*C1*C1 - 9*eg2*eg2) * (D*D*D*D) / 24
    Q4 = (61 + 90*T1 + 298*C1 + 45*T1*T1 - 3*C1*C1 - 252*eg2*eg2)*(D*D*D*D*D*D)/720
    lat = fp - Q1*(Q2 - Q3 + Q4)

    Q5 = D
    Q6 = (1 + 2*T1 + C1)*(D*D*D)/6
    Q7 = (5 - 2*C1 + 28*T1 - 3*C1*C1 + 8*eg2*eg2 + 24*T1*T1)*(D*D*D*D*D)/120
    lon = from_grid.lon0 + (Q5 - Q6 + Q7) / cosfp

    return (lat, lon)


def latlon_to_grid(lat: float, lon: float, from_datum: Datum, to_grid: Grid) -> Tuple[int, int]:
    """Lat/Lon to Local Grid conversion"""
    # Datum data for Lat/Lon to TM conversion
    a = from_datum.a
    e = from_datum.e
    b = from_datum.b

    # Lat/Lon -> TM
    slat1 = sin(lat)
    clat1 = cos(lat)
    clat1sq = clat1 * clat1
    tanlat1sq = slat1 * slat1 / clat1sq
    e2 = e * e
    e4 = e2 * e2
    e6 = e4 * e2
    eg = e * a / b
    eg2 = eg * eg

    l1 = 1 - e2/4 - 3*e4/64 - 5*e6/256
    l2 = 3*e2/8 + 3*e4/32 + 45*e6/1024
    l3 = 15*e4/256 + 45*e6/1024
    l4 = 35*e6/3072
    M = a*(l1*lat - l2*sin(2*lat) + l3*sin(4*lat) - l4*sin(6*lat))
    # double rho = a*(1-e2) / pow((1-(e*slat1)*(e*slat1)),1.5);
    nu = a / sqrt(1 - (e*slat1) * (e*slat1))
    p = lon - to_grid.lon0
    k0 = to_grid.k0
    # y = northing = K1 + K2p2 + K3p4, where
    K1 = M*k0
    K2 = k0*nu*slat1*clat1/2
    K3 = (k0*nu*slat1*clat1*clat1sq/24)*(5 - tanlat1sq + 9*eg2*clat1sq + 4*eg2*eg2*clat1sq*clat1sq)
    # ING north
    Y = K1 + K2*p*p + K3*p*p*p*p - to_grid.false_n

    # x = easting = K4p + K5p3, where
    K4 = k0*nu*clat1
    K5 = (k0*nu*clat1*clat1sq/6)*(1 - tanlat1sq + eg2*clat1*clat1)
    # ING east
    X = K4*p + K5*p*p*p + to_grid.false_e

    return (int(X + 0.5), int(Y + 0.5))


def Molodensky(lat: float, lon: float, from_datum: Datum, to_datum: Datum) -> Tuple[float, float]:
    """Abridged Molodensky transformation between 2 datums"""
    # from->WGS84 - to->WGS84 = from->WGS84 + WGS84->to = from->to
    dX = from_datum.dX - to_datum.dX
    dY = from_datum.dY - to_datum.dY
    dZ = from_datum.dZ - to_datum.dZ

    slat = sin(lat)
    clat = cos(lat)
    slon = sin(lon)
    clon = cos(lon)
    ssqlat = slat * slat

    # dlat = ((-dx * slat * clon - dy * slat * slon + dz * clat)
    #        + (da * rn * from_esq * slat * clat / from_a)
    #        + (df * (rm * adb + rn / adb )* slat * clat))
    #       / (rm + from.h);

    from_f = from_datum.f
    df = to_datum.f - from_f
    from_a = from_datum.a
    da = to_datum.a - from_a
    from_esq = from_datum.esq
    adb = 1.0 / (1.0 - from_f)
    rn = from_a / sqrt(1 - from_esq * ssqlat)
    rm = from_a * (1 - from_esq) / ((1 - from_esq * ssqlat)**1.5)
    from_h = 0.0

    dlat = (-dX*slat*clon - dY*slat*slon + dZ*clat
            + da*rn*from_esq*slat*clat/from_a
            + df*(rm*adb + rn/adb)*slat*clat) / (rm+from_h)

    olat = lat + dlat  # result (radians)

    # dlon = (-dx * slon + dy * clon) / ((rn + from.h) * clat)
    dlon = (-dX * slon + dY * clon) / ((rn + from_h) * clat)  # result (radians)
    olon = lon + dlon

    return (olat, olon)
