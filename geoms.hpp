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

#ifndef GEOMS_HPP_
#define GEOMS_HPP_

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
double getArea(Polygon2& poly) const
{
    int numVertices = poly.numVertices()
    if (numVertices < 3) { return 0; }
    double sum;
    for (int i = 0; i < numVertices-1; ++i)
    {
        sum += poly[i + 1][1] * poly[i][0] - poly[i + 1][0] * poly[i][1];
    }
    sum += poly[0][1] * poly[numVertices - 1][0] - poly[0][0] * poly[numVertices - 1][0];
    return sum / 2;
}

Point2 lineIntersect(Point2& l1, Point2& l2, Point2& p1, Point2& p2) const
{
    if ((p2[0] - p1[0]) * (l2[1] - l1[1]) - (p2[1] - p1[1]) * (l2[0] - l1[0]) != 0) //not parallel
    {
        double tp = ((p1[1] - l1[1]) * (l2[0] - l1[0]) - (p1[0] - l1[0]) * (l2[1] - l1[1])) / ((p2[0] - p1[0]) * (l2[1] - l1[1]) - (p2[1] - p1[1]) * (l2[0] - l1[0]));
        return Point2::lerp(p1, p2, tp);
    }
    else //parallel
    {
        return nullptr;
    }
}

std::vector < Point2 > getPolyLineIntersects(Polygon2& poly, Point2& l1, Point2& l2) const
{
    assert(l1 != l2);
    struct intersectDistances
    {
        Point2 intersect;
        double intersectDistance;
    };

    std::vector<intersectDistances> polyLineIntersects;
    int numVertices = poly.numVertices();
    if (l1[0] == l2[0])
    {
        for (int i = 0; i < numVertices; ++i)
        {
            Point2& p1 = poly[i];
            Point2& p2 = poly[(i + 1) % numVertices];
            if ((p1[0] < l1[0] && p2[0] > l1[0]) || (p1[0] > l1[0] && p2[0] < l1[0]))
            {
                t = (l1[0] - p1[0]) / (p2[0] - p1[0]);
                pinter = Point2::lerp(p1, p2, t);
                intersectDistances intersectDistance = {
                    pinter, //point
                    (l2[0]-l1[0])*(pinter[0]-l1[0]) + (l2[1]-l1[1])*(pinter[1]-l1[1]) //distance
                }
                polyLineIntersects.push_back(pinter);
            }
        }
    }
    else
    {
        auto l = [](double x) {
            return l1[1] + (l2[1] - l1[1]) * (x - l1[0]) / (l2[0] - l1[0]);
        }
        for (int i = 0; i < numVertices; ++i)
        {
            Point2& p1 = poly[i];
            Point2& p2 = poly[(i + 1) % numVertices];
            if ((p1[1] > l(p1[0]) && p2[1] < l(p2[0])) || (p1[1] < l(p1[0]) && p2[1] > l(p2[0])))
            {
                t = (p1[1] - l1[1] - (l2[1] - l1[1]) * (p1[0] - l1[0]) / (l2[0] - l1[0])) / ((l2[1] - l1[1]) * (p2[0] - p1[0]) / (l2[0] - l1[0]) - (p2[1] - p1[1));
                pinter = Point2::lerp(p1, p2, t);
                intersectDistances intersectDistance = {
                    pinter, //point
                    (l2[0]-l1[0])*(pinter[0]-l1[0]) + (l2[1]-l1[1])*(pinter[1]-l1[1]), //distance
                }
            }
            polyLineIntersects.push_back(pinter);
        }
    }
    std::sort(polyLineIntersects);
    std::vector<Point2> intersects;
    for (int i = 0; i < polyLineIntersects.size(); ++i)
    {
        intersects.push_back(polyLineIntersects.intersect);
    }
    return intersects;
}

double getPolyLineRegion(Polygon2& poly, Point2& l1, Point2& l2)
{
    assert(l1 != l2);
    std::vector<Point2> intersectRegion;
    int numVertices = poly.numVertices();
    if (l1[0] == l2[0])
    {
        for (int i = 0; i < numVertices; ++i)
        {
            Point2& p1 = poly[i];
            Point2& p2 = poly[(i + 1) % numVertices];
            if ((p1[0] <= l1[0] && l1[1] < l2[1]) || (p1[0] >= l1[0] && l1[1] > l2[1]))
            {
                intersectRegion.push_back(p1);
            }
            if ((p1[0] < l1[0] && p2[0] > l1[0]) || (p1[0] > l1[0] && p2[0] < l1[0]))
            {
                double t = (l1[0] - p1[0]) / (p2[0] - p1[0]);
                Point2 pinter = Point2::lerp(p1, p2, t);
                intersectRegion.push_back(pinter)
            }
        }
    }
    else
    {
        for (int i = 0; i < numVertices; ++i)
        {
            Point2& p1 = poly[i];
            Point2& p2 = poly[(i + 1) % numVertices];
            double line1 = l1[1] + (l2[1] - l1[1]) * (p1[0] - l1[0]) / (l2[0] - l1[0]);
            double line2 = l1[1] + (l2[1] - l1[1]) * (p2[0] - l1[0]) / (l2[0] - l1[0]);
            if ((p1[1] >= line1 && l1[0] < l2[0]) || (p1[1] <= line1 && l1[0] > l2[0]))
            {
                intersectRegion.push_back(p1);
            }
            if ((p1[1] > line1 && p2[1] < line2) || (p1[1] < line1 && p2[1] > line2))
            {
                double t = (p1[1] - l1[1] - (l2[1] - l1[1]) * (p1[0] - l1[0]) / (l2[0] - l1[0])) / ((l2[1] - l1[1]) * (p2[0] - p1[0]) / (l2[0] - l1[0]) - (p2[1] - p1[1]));
                Point2 pinter = Point2::lerp(p1, p2, t);
                intersectRegion.push_back(pinter);
            }
        }
    }
    return intersectRegion
}

bool pointInPoly(Point2& p, Polygon2& poly)
{
    int windNumber = 0;
    int numVertices = poly.numVertices();
    for (int i = 0; i < numVertices; ++i)
    {
        Point2& vcur = poly[i];
        Point2& vnext = poly[(i + 1) % numVertices];
        if (vcur[1] <= p[1] && vnext[1] > p[1] && (vnext[0] - vcur[0]) * (p1[0] - vcur[1]) - (vnext[1] - vcur[1]) * (p[0] - vcur[0]) > 0)
        {
            windNumber += 1
        }
        else if (vcur[1] > p[1] && vnext[1] <= p[1] && (vnext[0] - vcur[0]) * (p[1] - vcur[1]) - (vnext[1] - vcur[1]) * (p[0] - vcur[0]) < 0)
        {
            windNumber -= 1
        }
    }
    if (windNumber == 0)
    {
        return false;
    }
    return true;
}

std::vector<vector<Point2>> getPolyIntersectArea(Polygon2& poly1, Polygon2& poly2)
{
    int numVertices1 = poly1.numVertices();
    int numVertices2 = poly2.numVertices();
    bool testedintersects[numVertices1][numVertices2] = {};
    std::vector<std::vector<Point2>> intersects;

    Point2 intersectpoints[numVertices1][numVertices2] = {};
    double t1s[numVertices1][numVertices2] = {};
    double t2s[numVertices1][numVertices2] = {};

    for (int index1 = 0; index1 < numVertices1; ++index1)
    {
        for (int index2 = 0; index2 < numVertices2; ++index2)
        {
            Point2& p1cur = poly1[index1];
            Point2& p1next = poly1[(index1 + 1) % numVertices1];
            Point2& p2cur = poly1[index2];
            Point2& p2next = poly2[(index2 + 1) % numVertices2];
            if (p2next[0] - p2cur[0])*(p1next[1] - l1[1]) - (p2next[1] - p2cur[1]) * (p1next[0] - l1[0]) != 0:
                testt2 = ((p2cur[1] - p1cur[1]) * (p1next[0] - p1cur[0]) - (p2cur[0] - p1cur[0]) * (p1next[1] - p1cur[1])) / ((p2next[0] - p2cur[0]) * (p1next[1] - p1cur[1]) - (p2next[1] - p2cur[1]) * (p1next[0] - p1cur[0]))
                ttest1 = ((p1cur[1] - p2cur[1]) * (p2next[0] - p2cur[0]) - (p1cur[0] - p2cur[0]) * (p2next[1] - p2cur[1])) / ((p1next[0] - p1cur[0]) * (p2next[1] - p2cur[1]) - (p1next[1] - p1cur[1]) * (p2next[0] - p2cur[0]))
                if (testt1 >= 0 && testt1 <= 1 && testt2 >= 0 && testt2 <= 1)
                {
                    intersectpoints[index1][index2] = Point2::lerp(p2cur, p2next, testt2);
                    t1s[index1][index2] = testt1;
                    t2s[index1][index2] = testt2;
                }
        }
    }

    for (int index1 = 0; index1 < numVertices1; ++index1)
    {
        for (int index2 = 0; index2 < numVertices2; ++index2)
        {
            if (!testedintersects[index1][index2])
            {
                std::vector<Point2> vlist;
                vlist.push_back(intersectpoints[index1][index2]);
                int curi1 = index1;
                int curi2 = index2;
                double v10 = poly1[(index1 + 1) % numVertices1][0] - poly1[index1][0];
                double v11 = poly1[(index1 + 1) % numVertices1][1] - poly1[index1][1]];
                double v20 = poly2[(index2 + 1) % numVertices2][0] - poly2[index2][0];
                double v21 = poly2[(index2 + 1) % numVertices2][1] - poly2[index2][1]];
                if (v10 * v21 - v11 * v20 > 0)
                {
                    bool onOne = false;
                    double t = t2s[index1][index2];
                }
                else
                {
                    bool onOne = true;
                    double t = t1s[index1][index2];
                }

                while (vlist[0] != vlist[vlist.size() - 1] || vlist.size() < 2)
                {
                    if (onOne)
                    {
                        double mint = 2;
                        double mini = curi2;
                        for (int i = 0; i < numVertices2; ++i)
                        {
                            if (t1s[curi1][(curi2 + i + 1) % numVertices2] < mint && (t1s[curi1][(curi2 + i + 1) % numVertices2] > t || (t1s[curi1][(curi2 + i + 1) % numVertices2] == t && i < numVertices2 - 1)))
                            {
                                mini = (curi2 + i + 1) % numVertices2;
                                mint = t1s[curi1][mini];
                            }
                        }
                        if (mint == 2)
                        {
                            vlist.push_back(poly1[(curi1 + 1) % numVertices1]);
                            t = 0;
                            curi1 = (curi1 + 1) % numVertices1;
                        }
                        else
                        {
                            testedintersects[curi1][mini] = true;
                            vlist.push_back(intersectpoints[curi1][mini]);
                            t = t2s[curi1][mini];
                            curi2 = mini;
                            onOne = !onOne;
                        }
                    }
                    else
                    {
                        double mint = 2;
                        double mini = curi1;
                        for (int i = 0; i < numVertices1; ++i)
                        {
                            if (t2s[(curi1 + i + 1) % numVertices1][curi2] < mint && (t2s[(curi1 + i + 1) % numVertices1][curi2] > t || (t2s[(curi1 + i + 1) % numVertices1][curi2] == t && i < numVertices1 - 1)))
                            {
                                mini = (curi1 + i + 1) % numVertices1;
                                mint = t2s[mini][curi2];
                            }
                        }
                        if (mint == 2)
                        {
                            vlist.push_back(poly2[(curi2 + 1) % numVertices2]);
                            t = 0;
                            curi2 = (curi2 + 1) % numVertices2;
                        }
                        else
                        {
                            testedintersects[mini][curi2] = true;
                            vlist.push_back(intersectpoints[mini][curi2]);
                            t = t1s[mini][curi2];
                            curi1 = mini;
                            onOne = !onOne;
                        }
                    }
                    intersects.push_back(vlist);
                }
            }
            else
            {
                testedintersects[index1][index2] = true;
            }
        }
    }
    if (intersects.size() > 0)
    {
        return intersects;
    }
    else
    {
        bool poly1inpoly2 = pointInPoly(poly1[0], poly2);
        bool poly2inpoly1 = pointInPoly(poly2[0], poly1);
        if (poly1inpoly2 && !poly2inpoly1)
        {
            std::vector<Point2> vlist;
            for (int i = 0; i < numVertices1; ++i)
            {
                vlist.push_back(poly1[i]);
            }
            intersects.push_back(vlist);
            return intersects;
        }
        else if (poly2inpoly1 && !poly1inpoly2)
        {
            std::vector<Point2> vlist;
            for (int i = 0; i < numVertices2; ++i)
            {
                vlist.push_back(poly2[i]);
            }
            intersects.push_back(vlist);
            return intersects;
        }
        else
        {
            return intersects;
        }
    }
}

Polygon2 mergePolys(Polygon2& poly1, Polygon2& poly2)
{
    std::vector<Point2> newpoly;
    int numVertices1 = poly1.numVertices();
    int numVertices2 = poly2.numVertices();
    for (int i = 0; i < numVertices1; ++i)
    {
        Point2& vertex1 = poly1[i];
        bool inPoly2 = false;
        int index = nullptr;
        for (int j = 0; j < numVertices2; ++j)
        {
            Point2& checkvertex1 = poly2[j];
            if (checkvertex1 == vertex1)
            {
                inPoly2 = true;
                index = j;
            }
        }
        if (inPoly2 && poly1[(i + 1) % numVertices1] == poly2[(index - 1) % numVertices2])
        {
            for (int j = 0; j < numVertices2 - 2; ++j)
            {
                newpoly.push_back(poly2[(index + 1 + j) % numVertices2]);
            }
            for (int j = 0; j < numVertices1; ++j)
            {
                newpoly.push_back(poly1[(i + 1 + j) % numVertices1]);
            }
            int counter = 0;
            while (counter < newpoly.size())
            {
                if (newpoly[counter] == newpoly[(counter + 2) % newpoly.size()])
                {
                    if (counter == newpoly.size() - 1)
                    {
                        newpoly.erase(0, 2);
                        counter = std::max(counter - 3, 0);
                    }
                    else
                    {
                        newpoly.erase(counter, counter + 2);
                        counter = std::max(counter - 1, 0);
                    }
                }
                else
                {
                    counter += 1;
                }
            }
            Polygon2 retpoly;
            for (int j = 0; j < newpoly.size(); ++j)
            {
                retpoly.add(newpoly[j]);
            }
            return retpoly;
        }
    }
}

//getDistance implemented for points in squared_distance.hpp
//getCentroid implemented for points in poly.Centroid()

//editing end

} /* end namespace detail */
} /* end namespace primal */
} /* end namespace axom */

#endif // GEOMS_HPP_
