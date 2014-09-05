/*
    pbrt source code is Copyright(c) 1998-2014
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file was modified by K. Breeden from the reference implementation for CS348B - Sp2014.
    contact: kbreeden@stanford.edu
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

// cameras/realistic.h
#include "pbrt.h"
#include "camera.h"
#include "film.h"
#include <cfloat>

// a simple 2d point with float vals
struct Point2f {
  Point2f() { x = y = 0; }
  Point2f(float _x, float _y) { x = _x; y = _y; }
  float &operator[](int i) { return (i == 0) ? x : y; }
  float operator[](int i) const { return (i == 0) ? x : y; }
  float x, y;
};

// simple 2d bounding box defined by two Point2fs
struct Bounds2f {
    Bounds2f() {
        float minNum = FLT_MIN; // these values come from the cfloat library
        float maxNum = FLT_MAX;
        pMin = Point2f(maxNum, maxNum);
        pMax = Point2f(minNum, minNum);
    }

    // one point constructor, pmin = pmax = input
    Bounds2f(const Point2f &p) : pMin(p), pMax(p) { }

    // two point constructor
    Bounds2f(const Point2f &p1, const Point2f &p2) {
        pMin = Point2f(std::min(p1.x, p2.x), std::min(p1.y, p2.y));
        pMax = Point2f(std::max(p1.x, p2.x), std::max(p1.y, p2.y));
    }

    // is the given point inside the bounding area?
    bool Inside(const Point2f &pt) const {
        return (pt.x >= pMin.x && pt.x <= pMax.x &&
                pt.y >= pMin.y && pt.y <= pMax.y);
    }

    // returns the area of the bounding rectangle
    float Area() const {
        return ((pMax.x - pMin.x) * (pMax.y - pMin.y));
    }

    // returns a point that is linearly interpolated, componentwise, by (t_x, t_y)
    // between pMin and pMax.
    Point2f Lerp(const Point2f &t) const {
        return Point2f(::Lerp(t.x, pMin.x, pMax.x), ::Lerp(t.y, pMin.y, pMax.y));
    }

    // bump out pMin by (-delta, -delta) and pMax by (delta, delta)
    void Expand(float delta) {
      pMin.x -= delta;
      pMin.y -= delta;
      pMax.x += delta;
      pMax.y += delta;
    }

  // members
  Point2f pMin, pMax;
};

// returns a new Bounds2f which encompasses both of the input bounds
inline Bounds2f Union(const Bounds2f &b, const Bounds2f &c) {
    Bounds2f ret(Point2f(std::min(b.pMin.x, c.pMin.x),
                         std::min(b.pMin.y, c.pMin.y)),
                 Point2f(std::max(b.pMax.x, c.pMax.x),
                         std::max(b.pMax.y, c.pMax.y)));
    return ret;
}



// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
    // RealisticCamera Public Methods
    RealisticCamera(const AnimatedTransform &cam2world,
                    float sopen, float sclose, float apertureDiam,
                    float filmDist, float filmDiagonal,
                    float focusDistance, bool sw,
                    const char *lensFile, Film *film);

    float GenerateRay(const CameraSample &sample, Ray *) const;

private:
    // RealisticCamera Private Declarations
    struct LensElement {
        float curvatureRadius;      // each element is assumed to be spherical; if none is given, the element is planar.
        float thickness;            // distance from element surface to the next element along the z-axis
        float eta;                  // index of refraction
        float apertureRadius;       // diameter of each element
    };

    // RealisticCamera Private Data
    float apertureRadius, focusedFilmDistance, filmDiagonal;
    bool simpleWeighting;
    std::vector<LensElement> elements;
    float filmWidth, filmHeight;
    std::vector<Bounds2f> exitPupilBounds;

    // RealisticCamera Private Methods
    float FocusDistance(float filmDist) const;
    Bounds2f BoundExitPupil(const Point2f &pFilm, float filmDistance) const;
    void RenderExitPupil(float sx, float sy, const char *filename) const;
    bool TraceThroughLenses(Ray *ray, float filmDistance) const;
    Point SampleExitPupil(const Point2f &pFilm, const Point2f &lensSample,
                            float *pdf) const;
    void TestExitPupilBounds() const;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif // PBRT_CAMERAS_REALISTIC_H
