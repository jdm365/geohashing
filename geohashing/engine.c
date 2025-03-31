#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "engine.h"

static inline float LON_WRAP_ADD(float lon, float delta) {
	return (lon + delta) > 180.0 ? (lon + delta) - 360.0 : (lon + delta);
}
static inline float LON_WRAP_SUB(float lon, float delta) {
	return (lon - delta) < -180.0 ? (lon - delta) + 360.0 : (lon - delta);
}


inline double intBitsToDouble64(uint64_t i) {
	union {
		double f;
		uint64_t i;
	} u;
	u.i = i;
	return u.f;
}

inline float intBitsToFloat32(uint32_t i) {
	union {
		float f;
		uint32_t i;
	} u;
	u.i = i;
	return u.f;
}

inline void spreadBits64(uint64_t val, uint64_t *result, uint64_t offset) {
	uint64_t x = val;
	x = (x | (x << 16)) & 0x00FFFF0000FFFF00ULL;
	x = (x | (x << 8)) & 0xFF00FF00FF00FF00ULL;
	x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0FULL;
	x = (x | (x << 2)) & 0x3333333333333333ULL;
	x = (x | (x << 1)) & 0x5555555555555555ULL;

	*result |= x << offset;
}

inline uint64_t unspreadBits64(uint64_t x, uint64_t offset) {
	uint64_t result = x >> offset;
	result &= 0x5555555555555555ULL;

	result = (result | (result >> 1)) & 0x3333333333333333ULL;
	result = (result | (result >> 2)) & 0x0F0F0F0F0F0F0F0FULL;
	result = (result | (result >> 4)) & 0x00FF00FF00FF00FFULL;
	result = (result | (result >> 8)) & 0x0000FFFF0000FFFFULL;
	result = (result | (result >> 16)) & 0x0000000FFFFFFFFFULL;
	return result;
}

inline void spreadBits32(uint32_t val, uint32_t *result, int offset) {
	uint32_t x = val;
	x = (x | (x << 8)) & 0x00FF00FF;
	x = (x | (x << 4)) & 0x0F0F0F0F;
	x = (x | (x << 2)) & 0x33333333;
	x = (x | (x << 1)) & 0x55555555;

	*result |= x << offset;
}

inline uint32_t unspreadBits32(uint32_t x, int offset) {
	uint32_t result = x >> offset;
	result &= 0x55555555;

	result = (result | (result >> 1)) & 0x33333333;
	result = (result | (result >> 2)) & 0x0F0F0F0F;
	result = (result | (result >> 4)) & 0x00FF00FF;
	result = (result | (result >> 8)) & 0x0000FFFF;
	return result;
}

inline uint64_t calculate_bits_for_resolution64(double resolution_width_km) {
    // Approximate calculation:
    // Each bit halves the cell size.  We start with roughly half the circumference of the Earth
    // (at the equator), and keep dividing until we reach the desired resolution.

    double current_width = EARTH_RADIUS * PI;
    uint64_t bits = 0;

	if (current_width < resolution_width_km) {
		fprintf(stderr, "Resolution width is too small: %f\n", resolution_width_km);
		return UINT64_MAX;
	}

    while ((current_width > resolution_width_km) && (bits < 32)) {
        current_width *= 0.5;
        ++bits;
    }

    return bits;
}

inline uint64_t calculate_bits_for_resolution32(float resolution_width_km) {
    // Approximate calculation:
    // Each bit halves the cell size.  We start with roughly half the circumference of the Earth
    // (at the equator), and keep dividing until we reach the desired resolution.

    float current_width = (float)(EARTH_RADIUS * PI);
    uint64_t bits = 0;

	if (current_width < resolution_width_km) {
		fprintf(stderr, "Resolution width is too large: %f\n", resolution_width_km);
		fprintf(stderr, "Resolution must be less than or equal to the half of the earth's circumference: %f\n", current_width);
		exit(1);
	}

    while ((current_width > resolution_width_km) && (bits < 16)) {
        current_width *= 0.5;
        ++bits;
    }

    return bits;
}

inline uint64_t _geohash64(float lat, float lon, uint64_t shift_val) {
	double lat_norm = (lat + 90.0) * DOUBLE_DIV_180;
	double lon_norm = (lon + 180.0) * DOUBLE_DIV_360;

	uint64_t result = 0;
	spreadBits64((uint64_t)lon_norm, &result, 0);
	spreadBits64((uint64_t)lat_norm, &result, 1);

	result = (result >> shift_val) << shift_val;
	return result;
}

inline uint32_t _geohash32(float lat, float lon, uint64_t shift_val) {
	float lat_norm = (lat + 90.0) * FLOAT_DIV_180;
	float lon_norm = (lon + 180.0) * FLOAT_DIV_360;

	uint32_t result = 0;
	spreadBits32((uint32_t)lon_norm, &result, 0);
	spreadBits32((uint32_t)lat_norm, &result, 1);

	result = (result >> shift_val) << shift_val;
	return result;
}

inline void _geohash64Around(float lat, float lon, uint64_t shift_val, double resolution_width_km, uint64_t* hashes) {
	if (shift_val < 2) {
		hashes[0] = _geohash64(lat, lon, shift_val);
		return;
	}
	uint64_t next_resolution_hash = _geohash64(lat, lon, shift_val - 2);

	// Determine quadrant
	// -----------
	// | 10 | 11 |
	// |----|----|
	// | 00 | 01 |
	// -----------
	uint32_t last_two_bits = (next_resolution_hash >> (shift_val - 2)) & 0x3;

	float degree_approx = 111.0 / resolution_width_km;

	uint64_t base_hash = (next_resolution_hash >> shift_val) << shift_val;

	// TODO: Figure out faster version.
	// ALSO TODO: Handle wrapping around lat/lon boundaries.
    switch (last_two_bits) {
        case 0:
            hashes[0] = base_hash;
            hashes[1] = _geohash64(lat - degree_approx, lon, shift_val);
            hashes[2] = _geohash64(lat - degree_approx, LON_WRAP_SUB(lon, degree_approx), shift_val);
            hashes[3] = _geohash64(lat, LON_WRAP_SUB(lon, degree_approx), shift_val);
            break;
        case 1:
            hashes[0] = base_hash;
            hashes[1] = _geohash64(lat + degree_approx, lon, shift_val);
            hashes[2] = _geohash64(lat + degree_approx, LON_WRAP_SUB(lon, degree_approx), shift_val);
            hashes[3] = _geohash64(lat, LON_WRAP_SUB(lon, degree_approx), shift_val);
            break;
        case 2:
            hashes[0] = base_hash;
            hashes[1] = _geohash64(lat - degree_approx, lon, shift_val);
            hashes[2] = _geohash64(lat - degree_approx, LON_WRAP_ADD(lon, degree_approx), shift_val);
            hashes[3] = _geohash64(lat, LON_WRAP_ADD(lon, degree_approx), shift_val);
            break;
        case 3:
            hashes[0] = base_hash;
            hashes[1] = _geohash64(lat + degree_approx, lon, shift_val);
            hashes[2] = _geohash64(lat + degree_approx, LON_WRAP_ADD(lon, degree_approx), shift_val);
            hashes[3] = _geohash64(lat, LON_WRAP_ADD(lon, degree_approx), shift_val);
            break;
        default:
            hashes[0] = _geohash64(lat, lon, shift_val);
            break;
    }
}

inline void _geohash32Around(float lat, float lon, uint64_t shift_val, float resolution_width_km, uint32_t* hashes) {
	if (shift_val < 2) {
		hashes[0] = _geohash32(lat, lon, shift_val);
		return;
	}
	uint32_t next_resolution_hash = _geohash32(lat, lon, shift_val - 2);

	// Determine quadrant
	// -----------
	// | 10 | 11 |
	// |----|----|
	// | 00 | 01 |
	// -----------
	uint32_t last_two_bits = (next_resolution_hash >> (shift_val - 2)) & 0x3;

	float degree_approx = 111.0 / resolution_width_km;

	uint32_t base_hash = (next_resolution_hash >> shift_val) << shift_val;

	// TODO: Figure out faster version.
	// ALSO TODO: Handle wrapping around lat/lon boundaries.
    switch (last_two_bits) {
        case 0:
            hashes[0] = base_hash;
            hashes[1] = _geohash32(lat - degree_approx, lon, shift_val);
            hashes[2] = _geohash32(lat - degree_approx, LON_WRAP_SUB(lon, degree_approx), shift_val);
            hashes[3] = _geohash32(lat, LON_WRAP_SUB(lon, degree_approx), shift_val);
            break;
        case 1:
            hashes[0] = base_hash;
            hashes[1] = _geohash32(lat + degree_approx, lon, shift_val);
            hashes[2] = _geohash32(lat + degree_approx, LON_WRAP_SUB(lon, degree_approx), shift_val);
            hashes[3] = _geohash32(lat, LON_WRAP_SUB(lon, degree_approx), shift_val);
            break;
        case 2:
            hashes[0] = base_hash;
            hashes[1] = _geohash32(lat - degree_approx, lon, shift_val);
            hashes[2] = _geohash32(lat - degree_approx, LON_WRAP_ADD(lon, degree_approx), shift_val);
            hashes[3] = _geohash32(lat, LON_WRAP_ADD(lon, degree_approx), shift_val);
            break;
        case 3:
            hashes[0] = base_hash;
            hashes[1] = _geohash32(lat + degree_approx, lon, shift_val);
            hashes[2] = _geohash32(lat + degree_approx, LON_WRAP_ADD(lon, degree_approx), shift_val);
            hashes[3] = _geohash32(lat, LON_WRAP_ADD(lon, degree_approx), shift_val);
            break;
        default:
            hashes[0] = _geohash32(lat, lon, shift_val);
            break;
    }
}

inline uint64_t geohash64(float lat, float lon, double resolution_width_km) {
	if (lat < -90.0 || lat > 90.0) {
		fprintf(stderr, "Latitude out of range: %f\n", lat);
		exit(1);
	}
	if (lon < -180.0 || lon > 180.0) {
		fprintf(stderr, "Longitude out of range: %f\n", lon);
		exit(1);
	}

	// Multiply by 2 to account for lat and lon bits combined.
	uint64_t bits = calculate_bits_for_resolution64(resolution_width_km);
	uint64_t shift_val = 64 - (2 * bits);

	return _geohash64(lat, lon, shift_val);
}

inline uint32_t geohash32(float lat, float lon, float resolution_width_km) {
	if (lat < -90.0 || lat > 90.0) {
	fprintf(stderr, "Latitude out of range: %f\n", lat);
	return 0;
	}
	if (lon < -180.0 || lon > 180.0) {
	fprintf(stderr, "Longitude out of range: %f\n", lon);
	return 0;
	}

	// Multiply by 2 to account for lat and lon bits combined.
	uint64_t bits = calculate_bits_for_resolution32(resolution_width_km);
	uint64_t shift_val = 32 - (2 * bits);

	return _geohash32(lat, lon, shift_val);
}

inline void geohash64Around(float lat, float lon, double resolution_width_km, uint64_t* hashes) {
	if (lat < -90.0 || lat > 90.0) {
		fprintf(stderr, "Latitude out of range: %f\n", lat);
		exit(1);
	}
	if (lon < -180.0 || lon > 180.0) {
		fprintf(stderr, "Longitude out of range: %f\n", lon);
		exit(1);
	}

	// Multiply by 2 to account for lat and lon bits combined.
	uint64_t bits = calculate_bits_for_resolution64(resolution_width_km);
	uint64_t shift_val = 64 - (2 * bits);

	_geohash64Around(lat, lon, shift_val, resolution_width_km, hashes);
}

inline void geohash32Around(float lat, float lon, float resolution_width_km, uint32_t* hashes) {
	if (lat < -90.0 || lat > 90.0) {
		fprintf(stderr, "Latitude out of range: %f\n", lat);
		exit(1);
	}
	if (lon < -180.0 || lon > 180.0) {
		fprintf(stderr, "Longitude out of range: %f\n", lon);
		exit(1);
	}

	// Multiply by 2 to account for lat and lon bits combined.
	uint64_t bits = calculate_bits_for_resolution32(resolution_width_km);
	uint64_t shift_val = 32 - (2 * bits);

	_geohash32Around(lat, lon, shift_val, resolution_width_km, hashes);
}

Coord decodeGeohash64(uint64_t hash) {
	uint64_t lon_bits = unspreadBits64(hash, 0);
	uint64_t lat_bits = unspreadBits64(hash, 1);

	float lat = (float)(intBitsToDouble64(lat_bits) * _180_DIV_DOUBLE - 90.0);
	float lon = (float)(intBitsToDouble64(lon_bits) * _360_DIV_DOUBLE - 180.0);

	return (Coord){.lat = lat, .lon = lon};
}

Coord decodeGeohash32(uint32_t hash) {
	uint32_t lon_bits = unspreadBits32(hash, 0);
	uint32_t lat_bits = unspreadBits32(hash, 1);

	float lat = (float)(intBitsToFloat32(lat_bits) * _180_DIV_DOUBLE - 90.0);
	float lon = (float)(intBitsToFloat32(lon_bits) * _360_DIV_DOUBLE - 180.0);

	return (Coord){.lat = lat, .lon = lon};
}

inline float haversine(float lat1, float lon1, float lat2, float lon2) {
	float phi1 = lat1 * PI_DIV_180;
	float phi2 = lat2 * PI_DIV_180;
	float delta_phi = (lat2 - lat1) * PI_DIV_180;
	float delta_lambda = (lon2 - lon1) * PI_DIV_180;

	float a = sin(delta_phi * 0.5) * sin(delta_phi * 0.5) +
			cos(phi1) * cos(phi2) * sin(delta_lambda * 0.5) *
			sin(delta_lambda * 0.5);

	float c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

	return EARTH_RADIUS * c;
}

float decodeHaversine64(uint64_t hash1, uint64_t hash2) {
	Coord coords1 = decodeGeohash64(hash1);
	Coord coords2 = decodeGeohash64(hash2);

	return haversine(coords1.lat, coords1.lon, coords2.lat, coords2.lon);
}

float decodeHaversine32(uint32_t hash1, uint32_t hash2) {
	Coord coords1 = decodeGeohash32(hash1);
	Coord coords2 = decodeGeohash32(hash2);

	return haversine(coords1.lat, coords1.lon, coords2.lat, coords2.lon);
}

void bulkGeohash64(
		float* lats,
		float* lons,
		uint64_t* hashes,
		uint64_t num_coords,
		double resolution_width_km
		) {
	uint64_t bits = calculate_bits_for_resolution64(resolution_width_km);
	uint64_t shift_val = 64 - (2 * bits);

	uint64_t end_idx = num_coords - (num_coords % 8);

	for (uint64_t i = 0; i < end_idx; i += 8) {
		hashes[i] = _geohash64(lats[i], lons[i], shift_val);
		hashes[i + 1] = _geohash64(lats[i + 1], lons[i + 1], shift_val);
		hashes[i + 2] = _geohash64(lats[i + 2], lons[i + 2], shift_val);
		hashes[i + 3] = _geohash64(lats[i + 3], lons[i + 3], shift_val);
		hashes[i + 4] = _geohash64(lats[i + 4], lons[i + 4], shift_val);
		hashes[i + 5] = _geohash64(lats[i + 5], lons[i + 5], shift_val);
		hashes[i + 6] = _geohash64(lats[i + 6], lons[i + 6], shift_val);
		hashes[i + 7] = _geohash64(lats[i + 7], lons[i + 7], shift_val);
	}

	for (uint64_t i = end_idx; i < num_coords; ++i) {
		hashes[i] = _geohash64(lats[i], lons[i], shift_val);
	}
}


void bulkGeohash32(
		float* lats,
		float* lons,
		uint32_t* hashes,
		uint64_t num_coords,
		float resolution_width_km
		) {
	uint64_t bits = calculate_bits_for_resolution32(resolution_width_km);
	uint64_t shift_val = 32 - (2 * bits);

	printf("shift_val: %llu\n", shift_val);
	uint64_t end_idx = num_coords - (num_coords % 8);

	for (uint64_t i = 0; i < end_idx; i += 8) {
		hashes[i] = _geohash32(lats[i], lons[i], shift_val);
		hashes[i + 1] = _geohash32(lats[i + 1], lons[i + 1], shift_val);
		hashes[i + 2] = _geohash32(lats[i + 2], lons[i + 2], shift_val);
		hashes[i + 3] = _geohash32(lats[i + 3], lons[i + 3], shift_val);
		hashes[i + 4] = _geohash32(lats[i + 4], lons[i + 4], shift_val);
		hashes[i + 5] = _geohash32(lats[i + 5], lons[i + 5], shift_val);
		hashes[i + 6] = _geohash32(lats[i + 6], lons[i + 6], shift_val);
		hashes[i + 7] = _geohash32(lats[i + 7], lons[i + 7], shift_val);
	}

	for (uint64_t i = end_idx; i < num_coords; ++i) {
		hashes[i] = _geohash32(lats[i], lons[i], shift_val);
	}
}

void bulkGeohash64Around(
		float* lats,
		float* lons,
		uint64_t* hashes,
		uint64_t num_coords,
		double resolution_width_km
		) {
	uint64_t bits = calculate_bits_for_resolution32(resolution_width_km);
	uint64_t shift_val = 32 - (2 * bits);

	uint64_t end_idx = num_coords - (num_coords % 8);

	for (uint64_t i = 0; i < end_idx; i += 8) {
		_geohash64Around(lats[i], lons[i], shift_val, resolution_width_km, hashes + i);
		_geohash64Around(lats[i + 1], lons[i + 1], shift_val, resolution_width_km, hashes + (4 * (i + 1)));
		_geohash64Around(lats[i + 2], lons[i + 2], shift_val, resolution_width_km, hashes + (4 * (i + 2)));
		_geohash64Around(lats[i + 3], lons[i + 3], shift_val, resolution_width_km, hashes + (4 * (i + 3)));
		_geohash64Around(lats[i + 4], lons[i + 4], shift_val, resolution_width_km, hashes + (4 * (i + 4)));
		_geohash64Around(lats[i + 5], lons[i + 5], shift_val, resolution_width_km, hashes + (4 * (i + 5)));
		_geohash64Around(lats[i + 6], lons[i + 6], shift_val, resolution_width_km, hashes + (4 * (i + 6)));
		_geohash64Around(lats[i + 7], lons[i + 7], shift_val, resolution_width_km, hashes + (4 * (i + 7)));
	}

	for (uint64_t i = end_idx; i < num_coords; ++i) {
		_geohash64Around(lats[i], lons[i], shift_val, resolution_width_km, hashes + i);
	}
}

void bulkGeohash32Around(
		float* lats,
		float* lons,
		uint32_t* hashes,
		uint64_t num_coords,
		float resolution_width_km
		) {
	uint64_t bits = calculate_bits_for_resolution32(resolution_width_km);
	uint64_t shift_val = 32 - (2 * bits);

	uint64_t end_idx = num_coords - (num_coords % 8);

	for (uint64_t i = 0; i < end_idx; i += 8) {
		_geohash32Around(lats[i], lons[i], shift_val, resolution_width_km, hashes + i);
		_geohash32Around(lats[i + 1], lons[i + 1], shift_val, resolution_width_km, hashes + (4 * (i + 1)));
		_geohash32Around(lats[i + 2], lons[i + 2], shift_val, resolution_width_km, hashes + (4 * (i + 2)));
		_geohash32Around(lats[i + 3], lons[i + 3], shift_val, resolution_width_km, hashes + (4 * (i + 3)));
		_geohash32Around(lats[i + 4], lons[i + 4], shift_val, resolution_width_km, hashes + (4 * (i + 4)));
		_geohash32Around(lats[i + 5], lons[i + 5], shift_val, resolution_width_km, hashes + (4 * (i + 5)));
		_geohash32Around(lats[i + 6], lons[i + 6], shift_val, resolution_width_km, hashes + (4 * (i + 6)));
		_geohash32Around(lats[i + 7], lons[i + 7], shift_val, resolution_width_km, hashes + (4 * (i + 7)));
	}

	for (uint64_t i = end_idx; i < num_coords; ++i) {
		_geohash32Around(lats[i], lons[i], shift_val, resolution_width_km, hashes + i);
	}
}
