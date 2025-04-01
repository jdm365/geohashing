#include <math.h>
#include <stdint.h>
#include <stdio.h>

#define PI 3.14159265358979323846
#define PI_DIV_180 3.14159265358979323846 / 180.0
#define EARTH_RADIUS 6371.0

#define _180_DIV_DOUBLE (double)(180.0 / ((uint64_t)1 << 32))
#define _360_DIV_DOUBLE (double)(360.0 / ((uint64_t)1 << 32))
#define DOUBLE_DIV_180 (double)(((uint64_t)1 << 32) / 180.0)
#define DOUBLE_DIV_360 (double)(((uint64_t)1 << 32) / 360.0)

#define _180_DIV_FLOAT (float)(180.0 / (1 << 16))
#define _360_DIV_FLOAT (float)(360.0 / (1 << 16))
#define FLOAT_DIV_180 (float)(1 << 16) / 180.0
#define FLOAT_DIV_360 (float)(1 << 16) / 360.0

uint64_t doubleToIntBits64(double f);
double intBitsToFloat64(uint64_t i);

uint32_t floatToIntBits32(float f);
float intBitsToFloat32(uint32_t i);

void spreadBits64(uint64_t val, uint64_t *result, uint64_t offset);
uint64_t unspreadBits64(uint64_t x, uint64_t offset);

void spreadBits32(uint32_t val, uint32_t *result, int offset);
uint32_t unspreadBits32(uint32_t x, int offset);

uint64_t calculate_bits_for_resolution64(double resolution_width_km);
uint64_t calculate_bits_for_resolution32(float resolution_width_km);
uint64_t _geohash64(float lat, float lon, uint64_t shift_val);
uint32_t _geohash32(float lat, float lon, uint64_t shift_val);
uint64_t geohash64(float lat, float lon, double resolution_width_km);
uint32_t geohash32(float lat, float lon, float resolution_width_km);

void _geohash64Around(float lat, float lon, uint64_t shift_val, double resolution_width_km, uint64_t* hashes);
void _geohash32Around(float lat, float lon, uint64_t shift_val, float resolution_width_km, uint32_t* hashes);
void geohash64Around(float lat, float lon, double resolution_width_km, uint64_t* hashes);
void geohash32Around(float lat, float lon, float resolution_width_km, uint32_t* hashes);
void bulkGeohash64(float* lats, float* lons, uint64_t* hashes, uint64_t num_coords, double resolution_width_km);
void bulkGeohash32(float* lats, float* lons, uint32_t* hashes, uint64_t num_coords, float resolution_width_km, size_t num_threads);
void bulkGeohash32Around(float* lats, float* lons, uint32_t* hashes, uint64_t num_coords, float resolution_width_km);
void bulkGeohash64Around(float* lats, float* lons, uint64_t* hashes, uint64_t num_coords, double resolution_width_km);

typedef struct {
  float lat;
  float lon;
} Coord;

Coord decodeGeohash64(uint64_t hash);
Coord decodeGeohash32(uint32_t hash);

float _haversine(float lat1, float lon1, float lat2, float lon2);

float decodeHaversine64(uint64_t hash1, uint64_t hash2);
float decodeHaversine32(uint32_t hash1, uint32_t hash2);
