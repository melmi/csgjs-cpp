// Original CSG.JS library by Evan Wallace (http://madebyevan.com), under the MIT license.
// GitHub: https://github.com/evanw/csg.js/
// 
// C++ port by Tomasz Dabrowski (http://28byteslater.com), under the MIT license.
// GitHub: https://github.com/dabroz/csgjs-cpp/
// 
// Forked by Mohammad Elmi (http://melmi.ir), under the MIT license.
// GitHub: https://github.com/melmi/csgjs-cpp/
// 
// Constructive Solid Geometry (CSG) is a modeling technique that uses Boolean
// operations like union and intersection to combine 3D solids. This library
// implements CSG operations on meshes elegantly and concisely using BSP trees,
// and is meant to serve as an easily understandable implementation of the
// algorithm. All edge cases involving overlapping coplanar polygons in both
// solids are correctly handled.
//

#include <vector>

#ifndef _CSGJS_H_
#define _CSGJS_H_

namespace csgjs
{

struct Vector
{
	double x, y, z;

	Vector() : x(0.0), y(0.0), z(0.0) {}
	explicit Vector(double x, double y, double z) : x(x), y(y), z(z) {}
};

struct Vertex 
{
	Vector pos;
};

struct Model
{
	std::vector<Vertex> vertices;
	std::vector<int> indices;
};

// public interface - not super efficient, if you use multiple CSG operations you should
// use BSP trees and convert them into model only once. Another optimization trick is
// replacing model with your own class.

Model get_union(const Model & a, const Model & b);
Model get_intersection(const Model & a, const Model & b);
Model get_difference(const Model & a, const Model & b);

}

#endif

