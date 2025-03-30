from geohashing import (
        geohash, 
        decode_geohash32, 
        decode_geohash64,
        geohash_bulk,
        geohash_around_bulk,
        )
import numpy as np
from time import perf_counter


def validate():
    N = 1_000_000_000
    lats = np.random.uniform(-90, 90, N).astype(np.float32)
    lons = np.random.uniform(-180, 180, N).astype(np.float32)

    resolution_km = 1.00

    init = perf_counter()
    geohashes = geohash_bulk(lats, lons, resolution_km)
    ## geohashes = geohash_around_bulk(lats, lons, resolution_km)
    print(f"Time: {perf_counter() - init:.2f} seconds")
    print(f"Million geohashes per second: {(N / (perf_counter() - init) // 1_000_000):.2f}")

    print(geohashes)
    print(geohashes.shape, geohashes.dtype)
    print(f"Num unique geohashes: {len(np.unique(geohashes))}")



if __name__ == "__main__":
    validate()
