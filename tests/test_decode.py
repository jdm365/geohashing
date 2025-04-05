from geohashing import (
        geohash,
        decode_geohash32,
        decode_geohash64,
        haversine,
        geohash_bulk,
        geohash_around_bulk,
        )
import numpy as np

COORDS = [
    (37.421542, -122.085589),  # Googleplex
    (37.7749, -122.4194),  # San Francisco
    (34.0522, -118.2437),  # Los Angeles
    (40.7128, -74.0060),  # New York
    (51.5074, -0.1278),  # London
    (48.8566, 2.3522),  # Paris
    (55.7558, 37.6176),  # Moscow
    (35.6895, 139.6917),  # Tokyo
    (25.276987, 55.296249),  # Dubai
    (19.4326, -99.1332),  # Mexico City
    ]

def test_decode_32():
    resolution_km = 1.0
    tol = 1e-3

    hashed = geohash(*COORDS[0], resolution_km)
    for lat, lon in COORDS:
        hashed = geohash(lat, lon, resolution_km)

        decoded = decode_geohash32(hashed)

        distance_km = haversine(
                lat1=lat, 
                lon1=lon,
                lat2=decoded[0],
                lon2=decoded[1],
                )
        assert distance_km <= resolution_km + tol, f"Distance: {distance_km} km"

def test_decode_64():
    resolution_km = 0.001

    tol = 1e-2

    hashed = geohash(*COORDS[0], resolution_km)
    for lat, lon in COORDS:
        hashed = geohash(lat, lon, resolution_km)

        decoded = decode_geohash64(hashed)

        distance_km = haversine(
                lat1=lat, 
                lon1=lon,
                lat2=decoded[0],
                lon2=decoded[1],
                )
        assert distance_km <= resolution_km + tol, f"Distance: {distance_km} km"


def test_bulk_equivalence():
    resolution_km = 10.0

    tol = 1e-3

    N = 10_000_000
    lats = np.random.uniform(-90, 90, N).astype(np.float32)
    lons = np.random.uniform(-180, 180, N).astype(np.float32)

    hashed_bulk = geohash_bulk(lats, lons, resolution_km)
    hashed_around = geohash_around_bulk(lats, lons, resolution_km)

    for i in range(N):
        assert hashed_bulk[i] == hashed_around[i, 0], f"Bulk geohash mismatch at index {i}. " \
                f"Bulk: {hashed_bulk[i]}, Around: {hashed_around[i, 0]}. All around geohashes: {hashed_around[i]}"


if __name__ == '__main__':
    test_decode_32()
    test_decode_64()


    COORD_PAIR_1 = (42.227260, -83.321930)
    print(f"COORD_PAIR_1: {COORD_PAIR_1}")
    print(f"geohash: {geohash(COORD_PAIR_1[0], COORD_PAIR_1[1], 10.51)}")

    COORD_PAIR_2 = (42.227103, -83.32168449999)
    print(f"COORD_PAIR_2: {COORD_PAIR_2}")
    print(f"geohash: {geohash(COORD_PAIR_2[0], COORD_PAIR_2[1], 10.51)}")

    test_bulk_equivalence()
