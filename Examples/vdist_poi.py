#!/usr/bin/env python
"""
Example of using Google Maps queries and PyMap3D

https://developers.google.com/places/web-service/search

This requires a Google Cloud key, and costs a couple US cents per query.

TODO: Would like to instead query a larger region, would OSM be an option?
"""
import functools
from argparse import ArgumentParser
from pathlib import Path

import pandas
import requests
from pymap3d.vincenty import vdist

URL = "https://maps.googleapis.com/maps/api/place/nearbysearch/json?"


@functools.lru_cache()
def get_place_coords(
    place_type: str, latitude: float, longitude: float, search_radius_km: int, keyfn: Path
) -> pandas.DataFrame:
    """
    Get places using Google Maps Places API
    Requires you to have a Google Cloud account with API key.
    """

    keyfn = Path(keyfn).expanduser()
    key = keyfn.read_text()

    stub = URL + f"location={latitude},{longitude}"

    stub += f"&radius={search_radius_km * 1000}"

    stub += f"&types={place_type}"

    stub += f"&key={key}"

    r = requests.get(stub)
    r.raise_for_status()

    place_json = r.json()["results"]

    places = pandas.DataFrame(
        index=[p["name"] for p in place_json],
        columns=["latitude", "longitude", "distance_km", "vicinity"],
    )
    places["latitude"] = [p["geometry"]["location"]["lat"] for p in place_json]
    places["longitude"] = [p["geometry"]["location"]["lng"] for p in place_json]
    places["vicinity"] = [p["vicinity"] for p in place_json]

    return places


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument(
        "place_type",
        help="Place type to search: https://developers.google.com/places/supported_types",
    )
    p.add_argument(
        "searchloc", help="initial latituude, longitude to search from", nargs=2, type=float
    )
    p.add_argument("radius", help="search radius (kilometers)", type=int)
    p.add_argument("refloc", help="reference location (lat, lon)", nargs=2, type=float)
    p.add_argument("-k", "--keyfn", help="Google Places API key file", default="~/googlemaps.key")
    a = p.parse_args()

    place_coords = get_place_coords(a.place_type, *a.searchloc, a.radius, a.keyfn)

    place_coords["distance_km"] = (
        vdist(place_coords["latitude"], place_coords["longitude"], *a.refloc)[0] / 1e3
    )
