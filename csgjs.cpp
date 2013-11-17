
#include "csgjs.h"

#include <vector>
#include <algorithm>
#include <math.h>

namespace csgjs
{

// `CSG.Plane.EPSILON` is the tolerance used by `splitPolygon()` to decide if a
// point is on the plane.
static const double EPSILON = 0.00001;

struct polygon;

// Represents a plane in 3D space.
struct plane
{
	vector normal;
	double w;

	plane();
	plane(const vector & a, const vector & b, const vector & c);
	bool ok() const;
	void flip();
	void splitPolygon(const polygon & _polygon, std::vector<polygon> & coplanarFront, std::vector<polygon> & coplanarBack, std::vector<polygon> & front, std::vector<polygon> & back) const;
};

// Represents a convex polygon. The vertices used to initialize a polygon must
// be coplanar and form a convex loop. 
struct polygon
{
	std::vector<vertex> vertices;
	plane _plane;
	void flip();

	polygon();
	polygon(const std::vector<vertex> & list);
};

// Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
// by picking a polygon to split along. That polygon (and all other coplanar
// polygons) are added directly to that node and the other polygons are added to
// the front and/or back subtrees. This is not a leafy BSP tree since there is
// no distinction between internal and leaf nodes.
struct csgnode
{
	std::vector<polygon> polygons;
	csgnode * front;
	csgnode * back;
	plane _plane;

	csgnode();
	csgnode(const std::vector<polygon> & list);
	~csgnode();

	csgnode * clone() const;
	void clipTo(const csgnode * other);
	void invert();
	void build(const std::vector<polygon> & polygon);
	std::vector<polygon> clipPolygons(const std::vector<polygon> & list) const;
	std::vector<polygon> allPolygons() const;
};

// Vector implementation

inline static vector operator + (const vector & a, const vector & b) { return vector(a.x + b.x, a.y + b.y, a.z + b.z); }
inline static vector operator - (const vector & a, const vector & b) { return vector(a.x - b.x, a.y - b.y, a.z - b.z); }
inline static vector operator * (const vector & a, double b) { return vector(a.x * b, a.y * b, a.z * b); }
inline static vector operator / (const vector & a, double b) { return a * (1.0f / b); }
inline static double dot(const vector & a, const vector & b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline static vector lerp(const vector & a, const vector & b, double v) { return a + (b - a) * v; }
inline static vector negate(const vector & a) { return a * -1.0f; }
inline static double length(const vector & a) { return sqrt(dot(a, a)); }
inline static vector unit(const vector & a) { return a / length(a); }
inline static vector cross(const vector & a, const vector & b) { return vector(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }

// Vertex implementation

// Create a new vertex between this vertex and `other` by linearly
// interpolating all properties using a parameter of `t`. Subclasses should
// override this to interpolate additional properties.
inline static vertex interpolate(const vertex & a, const vertex & b, double t)
{
	vertex ret;
	ret.pos = lerp(a.pos, b.pos, t);
	return ret;
}

// Plane implementation

plane::plane() : normal(), w(0.0f) 
{
}

bool plane::ok() const 
{
	return length(this->normal) > 0.0f; 
}

void plane::flip()
{
	this->normal = negate(this->normal); 
	this->w *= -1.0f;
}

plane::plane(const vector & a, const vector & b, const vector & c)
{
	this->normal = unit(cross(b - a, c - a));
	this->w = dot(this->normal, a);
}

// Split `polygon` by this plane if needed, then put the polygon or polygon
// fragments in the appropriate lists. Coplanar polygons go into either
// `coplanarFront` or `coplanarBack` depending on their orientation with
// respect to this plane. Polygons in front or in back of this plane go into
// either `front` or `back`.
void plane::splitPolygon(const polygon & _polygon, std::vector<polygon> & coplanarFront, std::vector<polygon> & coplanarBack, std::vector<polygon> & front, std::vector<polygon> & back) const
{
	enum
	{
		COPLANAR = 0,
		FRONT = 1,
		BACK = 2,
		SPANNING = 3
	};

	// Classify each point as well as the entire polygon into one of the above
	// four classes.
	int polygonType = 0;
	std::vector<int> types;

	for (size_t i = 0; i < _polygon.vertices.size(); i++) 
	{
		double t = dot(this->normal, _polygon.vertices[i].pos) - this->w;
		int type = (t < -EPSILON) ? BACK : ((t > EPSILON) ? FRONT : COPLANAR);
		polygonType |= type;
		types.push_back(type);
	}

	// Put the polygon in the correct list, splitting it when necessary.
	switch (polygonType) 
	{
	case COPLANAR:
		{
			if (dot(this->normal, _polygon._plane.normal) > 0)
				coplanarFront.push_back(_polygon);
			else 
				coplanarBack.push_back(_polygon);
			break;
		}
	case FRONT:
		{
			front.push_back(_polygon);
			break;
		}
	case BACK:
		{
			back.push_back(_polygon);
			break;
		}
	case SPANNING:
		{
			std::vector<vertex> f, b;
			for (size_t i = 0; i < _polygon.vertices.size(); i++) 
			{
				int j = (i + 1) % _polygon.vertices.size();
				int ti = types[i], tj = types[j];
				vertex vi = _polygon.vertices[i], vj = _polygon.vertices[j];
				if (ti != BACK) f.push_back(vi);
				if (ti != FRONT) b.push_back(vi);
				if ((ti | tj) == SPANNING) 
				{
					double t = (this->w - dot(this->normal, vi.pos)) / dot(this->normal, vj.pos - vi.pos);
					vertex v = interpolate(vi, vj, t);
					f.push_back(v);
					b.push_back(v);
				}
			}
			if (f.size() >= 3) front.push_back(polygon(f));
			if (b.size() >= 3) back.push_back(polygon(b));
			break;
		}
	}
}

// Polygon implementation

void polygon::flip()
{
	std::reverse(vertices.begin(), vertices.end());
	_plane.flip();
}

polygon::polygon()
{
}

polygon::polygon(const std::vector<vertex> & list) : vertices(list), _plane(vertices[0].pos, vertices[1].pos, vertices[2].pos)
{
}

// Node implementation

// Return a new CSG solid representing space in either this solid or in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline static csgnode * csg_union(const csgnode * a1, const csgnode * b1)
{
	csgnode * a = a1->clone();
	csgnode * b = b1->clone();
	a->clipTo(b);
	b->clipTo(a);
	b->invert();
	b->clipTo(a);
	b->invert();
	a->build(b->allPolygons());
	csgnode * ret = new csgnode(a->allPolygons());
	delete a; a = 0;
	delete b; b = 0;
	return ret;
}

// Return a new CSG solid representing space in this solid but not in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline static csgnode * csg_subtract(const csgnode * a1, const csgnode * b1)
{
	csgnode * a = a1->clone();
	csgnode * b = b1->clone();
	a->invert();
	a->clipTo(b);
	b->clipTo(a);
	b->invert();
	b->clipTo(a);
	b->invert();
	a->build(b->allPolygons());
	a->invert();
	csgnode * ret = new csgnode(a->allPolygons());
	delete a; a = 0;
	delete b; b = 0;
	return ret;
}

// Return a new CSG solid representing space both this solid and in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline static csgnode * csg_intersect(const csgnode * a1, const csgnode * b1)
{
	csgnode * a = a1->clone();
	csgnode * b = b1->clone();
	a->invert();
	b->clipTo(a);
	b->invert();
	a->clipTo(b);
	b->clipTo(a);
	a->build(b->allPolygons());
	a->invert();
	csgnode * ret = new csgnode(a->allPolygons());
	delete a; a = 0;
	delete b; b = 0;
	return ret;
}

// Convert solid space to empty space and empty space to solid space.
void csgnode::invert()
{	
	for (size_t i = 0; i < this->polygons.size(); i++)
		this->polygons[i].flip();
	this->_plane.flip();
	if (this->front) this->front->invert();
	if (this->back) this->back->invert();
	std::swap(this->front, this->back);
}

// Recursively remove all polygons in `polygons` that are inside this BSP
// tree.
std::vector<polygon> csgnode::clipPolygons(const std::vector<polygon> & list) const
{
	if (!this->_plane.ok()) return list;
	std::vector<polygon> list_front, list_back;	
	for (size_t i = 0; i < list.size(); i++)
	{
		this->_plane.splitPolygon(list[i], list_front, list_back, list_front, list_back);
	}
	if (this->front) list_front = this->front->clipPolygons(list_front);
	if (this->back) list_back = this->back->clipPolygons(list_back);
	else list_back.clear();
	
	list_front.insert(list_front.end(), list_back.begin(), list_back.end());
	return list_front;
}

// Remove all polygons in this BSP tree that are inside the other BSP tree
// `bsp`.
void csgnode::clipTo(const csgnode * other)
{
	this->polygons = other->clipPolygons(this->polygons);
	if (this->front) this->front->clipTo(other);
	if (this->back) this->back->clipTo(other);
}

// Return a list of all polygons in this BSP tree.
std::vector<polygon> csgnode::allPolygons() const
{
	std::vector<polygon> list = this->polygons;
	std::vector<polygon> list_front, list_back;
	if (this->front) list_front = this->front->allPolygons();
	if (this->back) list_back = this->back->allPolygons();
	list.insert(list.end(), list_front.begin(), list_front.end());
	list.insert(list.end(), list_back.begin(), list_back.end());
	return list;
}

csgnode * csgnode::clone() const
{
	csgnode * ret = new csgnode();
	ret->polygons = this->polygons;
	ret->_plane = this->_plane;
	if (this->front) ret->front = this->front->clone();
	if (this->back) ret->back = this->back->clone();
	return ret;
}

// Build a BSP tree out of `polygons`. When called on an existing tree, the
// new polygons are filtered down to the bottom of the tree and become new
// nodes there. Each set of polygons is partitioned using the first polygon
// (no heuristic is used to pick a good split).
void csgnode::build(const std::vector<polygon> & list)
{
	if (!list.size()) return;
	if (!this->_plane.ok()) this->_plane = list[0]._plane;
	std::vector<polygon> list_front, list_back;
	for (size_t i = 0; i < list.size(); i++) 
	{
		this->_plane.splitPolygon(list[i], this->polygons, this->polygons, list_front, list_back);
	}
	if (list_front.size()) 
	{
		if (!this->front) this->front = new csgnode;
		this->front->build(list_front);
	}
	if (list_back.size()) 
	{
		if (!this->back) this->back = new csgnode;
		this->back->build(list_back);
	}
}

csgnode::csgnode() : front(0), back(0)
{
}

csgnode::csgnode(const std::vector<polygon> & list) : front(0), back(0)
{
	build(list);
}

csgnode::~csgnode()
{
	delete front;
	delete back;
}

// Public interface implementation

inline static std::vector<polygon> modelToPolygons(const model & model)
{
	std::vector<polygon> list;
	for (size_t i = 0; i < model.indices.size(); i+= 3)
	{
		std::vector<vertex> triangle;
		for (int j = 0; j < 3; j++)
		{
			vertex v = model.vertices[model.indices[i + j]];
			triangle.push_back(v);
		}
		list.push_back(polygon(triangle));
	}
	return list;
}

inline static model modelFromPolygons(const std::vector<polygon> & polygons)
{
	model model;
	int p = 0;
	for (size_t i = 0; i < polygons.size(); i++)
	{
		const polygon & poly = polygons[i];
		for (size_t j = 2; j < poly.vertices.size(); j++)
		{
			model.vertices.push_back(poly.vertices[0]);		model.indices.push_back(p++);
			model.vertices.push_back(poly.vertices[j - 1]);	model.indices.push_back(p++);
			model.vertices.push_back(poly.vertices[j]);		model.indices.push_back(p++);			
		}
	}
	return model;
}

typedef csgnode * csg_function(const csgnode * a1, const csgnode * b1);

inline static model operation(const model & a, const model & b, csg_function fun)
{
	csgnode * A = new csgnode(modelToPolygons(a));
	csgnode * B = new csgnode(modelToPolygons(b));
	csgnode * AB = fun(A, B);
	std::vector<polygon> polygons = AB->allPolygons();
	delete A; A = 0;
	delete B; B = 0;
	delete AB; AB = 0;
	return modelFromPolygons(polygons);
}

model get_union(const model & a, const model & b)
{
	return operation(a, b, csg_union);
}

model get_intersection(const model & a, const model & b)
{
	return operation(a, b, csg_intersect);
}

model get_difference(const model & a, const model & b)
{
	return operation(a, b, csg_subtract);
}

}
