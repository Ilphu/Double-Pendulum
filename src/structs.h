/**
 * @file structs.cc
 * @author Garrett
 * @brief 
 */



#ifndef STRUCT_GUARD
#define STRUCT_GUARD

#include <cmath>

using namespace std;

struct RectPoint {
    double x = 0;
    double y = 0;
};

struct PolarPoint {
    double r = 0;
    double theta = 0;
};

struct Pixel {
    uint8_t r = 0;
    uint8_t g = 0;
    uint8_t b = 0;
};

// struct RectPoint to_rectancular(const struct PolarPoint& pp) {
//     struct RectPoint pt;
//     pt.x = pp.r * cos(pp.theta);
//     pt.y = pp.r * sin(pp.theta);
//     return pt;
// }

// struct PolarPoint to_polar(const struct RectPoint& pt) {
//     struct PolarPoint pp;
//     pp.r = sqrt((pt.x * pt.x) + (pt.y * pt.y));
//     pp.theta = atan(pt.y / pt.x); // [-pi/2, pi/2]
//     return pp;
// }

#endif