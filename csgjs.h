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

struct Polygon;

// Represents a plane in 3D space.
struct Plane
{
	Vector normal;
	double w;

	Plane();
	Plane(const Vector & a, const Vector & b, const Vector & c);
	bool ok() const;
	void flip();
	void split_polygon(const Polygon & polygon, std::vector<Polygon> & coplanar_front, std::vector<Polygon> & coplanar_back, std::vector<Polygon> & front, std::vector<Polygon> & back) const;
};

// Represents a convex polygon. The vertices used to initialize a polygon must
// be coplanar and form a convex loop. 
struct Polygon
{
	std::vector<Vertex> vertices;
	Plane plane;
	void flip();

	Polygon();
	Polygon(const std::vector<Vertex> & list);
};

struct Model
{
	std::vector<Vertex> vertices;
	std::vector<int> indices;
};

Model model_from_polygons(const std::vector<Polygon> & polygons);
std::vector<Polygon> model_to_polygons(const Model & model);

// public interface - not super efficient, if you use multiple CSG operations you should
// use BSP trees and convert them into model only once. Another optimization trick is
// replacing model with your own class.

Model get_union(const Model & a, const Model & b);
Model get_intersection(const Model & a, const Model & b);
Model get_difference(const Model & a, const Model & b);

}

#endif

