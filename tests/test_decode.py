from geohashing import (
        geohash,
        decode_geohash32,
        decode_geohash64,
        haversine,
        )

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


if __name__ == '__main__':
    test_decode_32()
    test_decode_64()
