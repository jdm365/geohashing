cimport cython

from libc.stdint cimport uint32_t, uint64_t

cimport numpy as np
import numpy as np
np.import_array()



cdef extern from "engine.h":
    ctypedef struct Coord:
        float lat
        float lon

    uint32_t geohash32(
            float lat,
            float lon,
            float resolution_width_km
    )
    uint64_t geohash64(
            float lat,
            float lon,
            float resolution_width_km
    )
    void bulkGeohash32(
            float* lats,
            float* lons,
            uint32_t* hashes,
            size_t num_points,
            float resolution_width_km
    )
    void bulkGeohash64(
            float* lats,
            float* lons,
            uint64_t* hashes,
            size_t num_points,
            float resolution_width_km
    )
    void geohash32Around(
            float lat,
            float lon,
            float resolution_width_km,
            uint32_t* hashes
    )
    void geohash64Around(
            float lat,
            float lon,
            float resolution_width_km,
            uint64_t* hashes
    )
    void bulkGeohash32Around(
            float* lats,
            float* lons,
            uint32_t* hashes,
            size_t num_points,
            float resolution_width_km
    )
    void bulkGeohash64Around(
            float* lats,
            float* lons,
            uint64_t* hashes,
            size_t num_points,
            float resolution_width_km
    )
    Coord decodeGeohash32(uint32_t hash)
    Coord decodeGeohash64(uint64_t hash)
    float haversine(
            float lat1,
            float lon1,
            float lat2,
            float lon2
    )

    float decodeHaversine32(
            uint32_t hash1,
            uint32_t hash2
    )
    float decodeHaversine64(
            uint64_t hash1,
            uint64_t hash2
    )



def geohash_bulk(
        np.ndarray[float, ndim=1] lats,
        np.ndarray[float, ndim=1] lons,
        float resolution_width_km,
        ):
    if resolution_width_km < 0.01:
        return geohash_bulk_64(lats, lons, <double>resolution_width_km)
    else:
        return geohash32_bulk(lats, lons, resolution_width_km)

cdef np.ndarray[uint32_t, ndim=1] geohash32_bulk(
        np.ndarray[float, ndim=1] lats,
        np.ndarray[float, ndim=1] lons,
        float resolution_width_km,
        ):
    cdef np.ndarray[uint32_t, ndim=1] hashes = np.zeros(lats.shape[0], dtype=np.uint32)
    cdef float* lats_ptr = <float*>lats.data
    cdef float* lons_ptr = <float*>lons.data
    cdef uint32_t* hashes_ptr = <uint32_t*>hashes.data
    cdef uint64_t num_points = lats.shape[0]

    bulkGeohash32(lats_ptr, lons_ptr, hashes_ptr, num_points, resolution_width_km)
    return hashes

cdef np.ndarray[uint64_t, ndim=1] geohash_bulk_64(
        np.ndarray[float, ndim=1] lats,
        np.ndarray[float, ndim=1] lons,
        double resolution_width_km,
        ):
    cdef np.ndarray[uint64_t, ndim=1] hashes = np.zeros(lats.shape[0], dtype=np.uint64)
    cdef float* lats_ptr = <float*>lats.data
    cdef float* lons_ptr = <float*>lons.data
    cdef uint64_t* hashes_ptr = <uint64_t*>hashes.data
    cdef uint64_t num_points = lats.shape[0]

    bulkGeohash64(lats_ptr, lons_ptr, hashes_ptr, num_points, resolution_width_km)
    return hashes

def geohash_around_bulk(
        np.ndarray[float, ndim=1] lats,
        np.ndarray[float, ndim=1] lons,
        float resolution_width_km,
        ):
    if resolution_width_km < 0.01:
        return geohash_bulk_64_around(lats, lons, <double>resolution_width_km).reshape(-1, 4)
    else:
        return geohash32_bulk_around(lats, lons, resolution_width_km).reshape(-1, 4)

cdef np.ndarray[uint32_t, ndim=1] geohash32_bulk_around(
        np.ndarray[float, ndim=1] lats,
        np.ndarray[float, ndim=1] lons,
        float resolution_width_km,
        ):
    cdef np.ndarray[uint32_t, ndim=1] hashes = np.zeros(lats.shape[0] * 4, dtype=np.uint32)
    cdef float* lats_ptr = <float*>lats.data
    cdef float* lons_ptr = <float*>lons.data
    cdef uint32_t* hashes_ptr = <uint32_t*>hashes.data
    cdef uint64_t num_points = lats.shape[0]

    bulkGeohash32Around(lats_ptr, lons_ptr, hashes_ptr, num_points, resolution_width_km)
    return hashes

cdef np.ndarray[uint64_t, ndim=1] geohash_bulk_64_around(
        np.ndarray[float, ndim=1] lats,
        np.ndarray[float, ndim=1] lons,
        double resolution_width_km,
        ):
    cdef np.ndarray[uint64_t, ndim=1] hashes = np.zeros(lats.shape[0] * 4, dtype=np.uint64)
    cdef float* lats_ptr = <float*>lats.data
    cdef float* lons_ptr = <float*>lons.data
    cdef uint64_t* hashes_ptr = <uint64_t*>hashes.data
    cdef uint64_t num_points = lats.shape[0]

    bulkGeohash64Around(lats_ptr, lons_ptr, hashes_ptr, num_points, resolution_width_km)
    return hashes

def geohash(
        float lat,
        float lon,
        float resolution_width_km
        ):
    if resolution_width_km < 0.01:
        return geohash64(lat, lon, resolution_width_km)
    else:
        return geohash32(lat, lon, resolution_width_km)

def decode_geohash32(uint32_t _hash):
    cdef Coord coord = decodeGeohash32(_hash)
    return coord.lat, coord.lon

def decode_geohash64(uint64_t _hash):
    cdef Coord coord = decodeGeohash64(_hash)
    return coord.lat, coord.lon

def haversine(
        float lat1,
        float lon1,
        float lat2,
        float lon2
        ):
    return haversine(lat1, lon1, lat2, lon2)


def decode_haversine32(
        uint32_t hash1,
        uint32_t hash2
        ):
    return decodeHaversine32(hash1, hash2)

def decode_haversine64(
        uint64_t hash1,
        uint64_t hash2
        ):
    return decodeHaversine64(hash1, hash2)
