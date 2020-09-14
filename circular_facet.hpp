// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file geoms.hpp
 *
 * This file provides helper functions for testing whether
 * geometric primitives intersect
 */

#ifndef CIRCULAR_FACET_HPP_
#define CIRCULAR_FACET_HPP_

#include "axom/core/Macros.hpp"
#include "axom/core/numerics/Determinants.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "axom/primal/operators/squared_distance.hpp"

#include <vector>
#include <ostream>
#include <cmath>

namespace axom
{
namespace primal
{
namespace detail
{

//---------------------------- FUNCTION DECLARATIONS ---------------------------

typedef Vector2 = primal::Vector< double, 2 >;
typedef Point2 = primal::Point< double, 2 >;
typedef Polygon2 = primal::Polygon< double, 2 >;

//editing here
AXOM_HOST_DEVICE
double getArcArea(Point2& p1, Point2& p2, Point2& center, double radius)
{
    double distance = squared_distance(p1, p2);
    double theta = asin(distance / (2 * radius));
    double area = (theta * pow(radius, 2)) - (distance * cos(theta) * radius / 2);
    if ((p1[0] - center[0]) * (p2[1] - p1[1]) - (p1[1] - center[1]) * (p2[0] - p1[0]) >= 0)
    {
        return area;
    }
    return -area;
}

std::vector<Point2> getCircleLineIntersects(Point2& l1, Point2& l2, Point2& center, double radius)
{
    std::vector<Point2> intersects;
    double a = pow(l2[0] - l1[0], 2) + pow(l2[1] - l1[1], 2);
    double b = 2 * ((l1[0] - center[0]) * (l2[0] - l1[0]) + (l1[1] - center[1]) * (l2[1] - l1[1]));
    double c = pow(l1[0] - center[0], 2) + pow(l1[1] - center[1], 2) - pow(radius, 2);
    double disc = pow(b, 2) - 4 * a * c;
    if (disc > 0)
    {
        double x1 = (-b - sqrt(disc)) / (2 * a);
        if (x1 <= 1 && x1 > 0)
        {
            Point2 inter1 = lerp(l1, l2, x1);
            intersects.push_back(inter1);
        }
        double x2 = (-b + sqrt(disc)) / (2 * a);
        if (x2 < 1 && x2 >= 0)
        {
            Point2 inter2 = lerp(l1, l2, x2);
            intersects.push_back(inter2);
        }
    }
    return intersects;
}

double getCircleIntersectArea(Point2& center, double radius, Polygon2& poly)
{
    int numVertices = poly.numVertices();
    double adjustcorneramount = 1e-14;
    bool notmod = true;
    while (notmod)
    {
        int startAt = 1;
        std::vector<Point2> intersectpoints;
        std::vector<Point2> arcpoints;
        for (int i = 0; i < numVertices; ++i)
        {
            Point2 curpoint = poly[i];
            Point2 nextpoint = poly[(i + 1) % numVertices];
            bool curin = (squared_distance(curpoint, center) <= radius);
            bool nextin = (squared_distance(nextpoint, center) <= radius);
            if (curin)
            {
                intersectpoints.push_back(curpoint);
            }
            std::vector<Point2> lineintersects = getCircleLineIntersects(curpoint, nextpoint, center, radius);
            for (Point2& intersect : lineintersects)
            {
                intersectpoints.push_back(intersect);
                if (arcpoints.size() == 0 && curin && !nextin)
                {
                    startAt = 0;
                }
                arcpoints.push_back(intersect);
            }
        }
        if (arcpoints.size() % 2 == 1)
        {
            for (Point2 p : poly)
            {
                p[0] += adjustcorneramount;
                p[1] += adjustcorneramount;
            }
        }
        else
        {
            notmod = false;
        }
    }

    double area = 0;
    for (int i = 0; i < arcpoints.size(); i += 2)
    {
        area += getArcArea(arcpoints[startAt + i], arcpoints[(startAt + i + 1) % arcpoints.size()], center, radius);
    }
    area += getArea(intersectpoints);

    if (arcpoints.size() == 0 && area == 0)
    {
        if pointInPoly(center, poly)
        {
            return radius * radius * asin(1.0) * 2;
        }
        return 0;
    }
    return area % (radius * radius * asin(1.0) * 2); //arcpoints also needed
}

//editing end

} /* end namespace detail */
} /* end namespace primal */
} /* end namespace axom */

#endif // CIRCULAR_FACET_HPP_
