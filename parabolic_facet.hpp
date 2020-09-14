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

#ifndef PARABOLIC_FACET_HPP_
#define PARABOLIC_FACET_HPP_

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
double getParabolaArcArea(double a, double b, double c, double x1, double x2)
{
    if (x1 <= x2)
    {
        return (x2 - x1) / 2 * (a * (pow(x1, 2) + pow(x2, 2)) + b * (x1 + x2) + 2 * c) - (a * (pow(x2, 3) - pow(x1, 3)) / 3 + b * (pow(x2, 2) - pow(x1, 2)) / 2 + c * (x2 - x1));
    }
    else
    {
        return -getParabolaArcArea(a, b, c, x2, x1);
    }
}

std::vector<Point2> getParabolaLineIntersects(Point2& l1, Point2& l2, double a, double b, double c)
{
    std::vector<Point2> intersects;
    double qa = pow(a * (l2[0] - l1[0]), 2);
    double qb = 2 * a * l1[0] * (l2[0] - l1[0]) + b * (l2[0] - l1[0]) + l1[1] - l2[1];
    double qc = a * l1[0] * *2 + b * l1[0] + c - l1[1];
    if (qa == 0 && qb != 0)
    {
        double x = -qc / qb;
        if (x < 1 && x > 0)
        {
            Point2 inter = Point2::lerp(l1, l2, x);
            intersects.push_back(inter);
        }
        else if (qa != 0)
        {
            double disc = pow(qb, 2) - 4 * qa * qc;
            if (disc > 0)
            {
                double x1 = (-qb - sqrt(disc)) / (2 * qa);
                if (x1 <= 1 && x1 > 0)
                {
                    Point2 inter1 = Point::lerp(l1, l2, x1);
                    intersects.push_back(inter1);
                }
                double x2 = (-qb + sqrt(disc)) / (2 * qa);
                if (x2 < 1 && x2 >= 0)
                {
                    Point2 inter2 = Point::lerp(l1, l2, x2);
                    intersects.push_back(inter2);
                }
            }
        }
    }
    return intersects;
}

double getParabolaFacetArea(Polygon2& poly, double a, double b, double c)
{
    int numVertices = poly.numVertices();
    double adjustcorneramount = 1e-14;
    bool notmod = true;
    while (notmod)
    {
        int startAt = 1;
        std::vector<Point2> intersectpoints;
        std::vector<Point2> parabolapoints;
        for (int i = 0; i < numVertices; ++i)
        {
            Point2 curpoint = poly[i];
            Point2 nextpoint = poly[(i + 1) % numVertices];
            bool curin = (curpoint[1] >= a * pow(curpoint[0], 2) + b * curpoint[0] + c);
            bool nextin = (nextpoint[1] >= a * pow(nextpoint[0], 2) + b * nextpoint[0] + c);
            if (curin)
            {
                intersectpoints.push_back(curpoint);
            }
            std::vector<Point2> lineintersects = getParabolaLineIntersects(curpoint, nextpoint, a, b, c);
            for (Point2& intersect : lineintersects)
            {
                intersectpoints.push_back(intersect);
                if (parabolapoints.size() == 0 && curin && !nextin)
                {
                    startAt = 0;
                }
                parabolapoints.push_back(intersect);
            }
        }
        if (parabolapoints.size() % 2 == 1)
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
    for (int i = 0; i < parabolapoints.size(); i += 2)
    {
        area += getParabolaArcArea(a, b, c, parabolapoints[startAt + i][0], parabolapoints[(startAt + i + 1) % parabolapoints.size()][0]);
    }
    area += geoms::getArea(intersectpoints);

    return area //parabolapoints also needed
}

//editing end

} /* end namespace detail */
} /* end namespace primal */
} /* end namespace axom */

#endif // PARABOLIC_FACET_HPP_
