
#include "csgjs.h"

#include <vector>
#include <algorithm>
#include <math.h>

namespace csgjs
{

// `CSG.Plane.EPSILON` is the tolerance used by `split_polygon()` to decide if a
// point is on the plane.
static const double EPSILON = 0.00001;

// Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
// by picking a polygon to split along. That polygon (and all other coplanar
// polygons) are added directly to that node and the other polygons are added to
// the front and/or back subtrees. This is not a leafy BSP tree since there is
// no distinction between internal and leaf nodes.
struct Node
{
	std::vector<Polygon> polygons;
	Node * front;
	Node * back;
	Plane plane;

	Node();
	Node(const std::vector<Polygon> & list);
	~Node();

	Node * clone() const;
	void clip_to(const Node * other);
	void invert();
	void build(const std::vector<Polygon> & polygon);
	std::vector<Polygon> clip_polygons(const std::vector<Polygon> & list) const;
	std::vector<Polygon> all_polygons() const;
};

// Vector implementation

inline static Vector operator + (const Vector & a, const Vector & b) { return Vector(a.x + b.x, a.y + b.y, a.z + b.z); }
inline static Vector operator - (const Vector & a, const Vector & b) { return Vector(a.x - b.x, a.y - b.y, a.z - b.z); }
inline static Vector operator * (const Vector & a, double b) { return Vector(a.x * b, a.y * b, a.z * b); }
inline static Vector operator / (const Vector & a, double b) { return a * (1.0 / b); }
inline static double dot(const Vector & a, const Vector & b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline static Vector lerp(const Vector & a, const Vector & b, double v) { return a + (b - a) * v; }
inline static Vector negate(const Vector & a) { return a * -1.0; }
inline static double length(const Vector & a) { return sqrt(dot(a, a)); }
inline static Vector unit(const Vector & a) { return a / length(a); }
inline static Vector cross(const Vector & a, const Vector & b) { return Vector(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }

// Vertex implementation

// Create a new vertex between this vertex and `other` by linearly
// interpolating all properties using a parameter of `t`. Subclasses should
// override this to interpolate additional properties.
inline static Vertex interpolate(const Vertex & a, const Vertex & b, double t)
{
	Vertex ret;
	ret.pos = lerp(a.pos, b.pos, t);
	return ret;
}

// Plane implementation

Plane::Plane() : normal(), w(0.0) 
{
}

bool Plane::ok() const 
{
	return length(this->normal) > 0.0; 
}

void Plane::flip()
{
	this->normal = negate(this->normal); 
	this->w *= -1.0;
}

Plane::Plane(const Vector & a, const Vector & b, const Vector & c)
{
	this->normal = unit(cross(b - a, c - a));
	this->w = dot(this->normal, a);
}

// Split `polygon` by this plane if needed, then put the polygon or polygon
// fragments in the appropriate lists. Coplanar polygons go into either
// `coplanar_front` or `coplanar_back` depending on their orientation with
// respect to this plane. Polygons in front or in back of this plane go into
// either `front` or `back`.
void Plane::split_polygon(const Polygon & polygon, std::vector<Polygon> & coplanar_front, std::vector<Polygon> & coplanar_back, std::vector<Polygon> & front, std::vector<Polygon> & back) const
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
	int polygon_type = 0;
	std::vector<int> types;

	for (size_t i = 0; i < polygon.vertices.size(); i++) 
	{
		double t = dot(this->normal, polygon.vertices[i].pos) - this->w;
		int type = (t < -EPSILON) ? BACK : ((t > EPSILON) ? FRONT : COPLANAR);
		polygon_type |= type;
		types.push_back(type);
	}

	// Put the polygon in the correct list, splitting it when necessary.
	switch (polygon_type) 
	{
	case COPLANAR:
		{
			if (dot(this->normal, polygon.plane.normal) > 0)
				coplanar_front.push_back(polygon);
			else 
				coplanar_back.push_back(polygon);
			break;
		}
	case FRONT:
		{
			front.push_back(polygon);
			break;
		}
	case BACK:
		{
			back.push_back(polygon);
			break;
		}
	case SPANNING:
		{
			std::vector<Vertex> f, b;
			for (size_t i = 0; i < polygon.vertices.size(); i++) 
			{
				int j = (i + 1) % polygon.vertices.size();
				int ti = types[i], tj = types[j];
				Vertex vi = polygon.vertices[i], vj = polygon.vertices[j];
				if (ti != BACK) f.push_back(vi);
				if (ti != FRONT) b.push_back(vi);
				if ((ti | tj) == SPANNING) 
				{
					double t = (this->w - dot(this->normal, vi.pos)) / dot(this->normal, vj.pos - vi.pos);
					Vertex v = interpolate(vi, vj, t);
					f.push_back(v);
					b.push_back(v);
				}
			}
			if (f.size() >= 3) front.push_back(Polygon(f));
			if (b.size() >= 3) back.push_back(Polygon(b));
			break;
		}
	}
}

// Polygon implementation

void Polygon::flip()
{
	std::reverse(vertices.begin(), vertices.end());
	plane.flip();
}

Polygon::Polygon()
{
}

Polygon::Polygon(const std::vector<Vertex> & list) : vertices(list), plane(vertices[0].pos, vertices[1].pos, vertices[2].pos)
{
}

// Node implementation

// Return a new CSG solid representing space in either this solid or in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline static Node * csg_union(const Node * a1, const Node * b1)
{
	Node * a = a1->clone();
	Node * b = b1->clone();
	a->clip_to(b);
	b->clip_to(a);
	b->invert();
	b->clip_to(a);
	b->invert();
	a->build(b->all_polygons());
	Node * ret = new Node(a->all_polygons());
	delete a; a = 0;
	delete b; b = 0;
	return ret;
}

// Return a new CSG solid representing space in this solid but not in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline static Node * csg_subtract(const Node * a1, const Node * b1)
{
	Node * a = a1->clone();
	Node * b = b1->clone();
	a->invert();
	a->clip_to(b);
	b->clip_to(a);
	b->invert();
	b->clip_to(a);
	b->invert();
	a->build(b->all_polygons());
	a->invert();
	Node * ret = new Node(a->all_polygons());
	delete a; a = 0;
	delete b; b = 0;
	return ret;
}

// Return a new CSG solid representing space both this solid and in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline static Node * csg_intersect(const Node * a1, const Node * b1)
{
	Node * a = a1->clone();
	Node * b = b1->clone();
	a->invert();
	b->clip_to(a);
	b->invert();
	a->clip_to(b);
	b->clip_to(a);
	a->build(b->all_polygons());
	a->invert();
	Node * ret = new Node(a->all_polygons());
	delete a; a = 0;
	delete b; b = 0;
	return ret;
}

// Convert solid space to empty space and empty space to solid space.
void Node::invert()
{	
	for (size_t i = 0; i < this->polygons.size(); i++)
		this->polygons[i].flip();
	this->plane.flip();
	if (this->front) this->front->invert();
	if (this->back) this->back->invert();
	std::swap(this->front, this->back);
}

// Recursively remove all polygons in `polygons` that are inside this BSP
// tree.
std::vector<Polygon> Node::clip_polygons(const std::vector<Polygon> & list) const
{
	if (!this->plane.ok()) return list;
	std::vector<Polygon> list_front, list_back;	
	for (size_t i = 0; i < list.size(); i++)
	{
		this->plane.split_polygon(list[i], list_front, list_back, list_front, list_back);
	}
	if (this->front) list_front = this->front->clip_polygons(list_front);
	if (this->back) list_back = this->back->clip_polygons(list_back);
	else list_back.clear();
	
	list_front.insert(list_front.end(), list_back.begin(), list_back.end());
	return list_front;
}

// Remove all polygons in this BSP tree that are inside the other BSP tree
// `bsp`.
void Node::clip_to(const Node * other)
{
	this->polygons = other->clip_polygons(this->polygons);
	if (this->front) this->front->clip_to(other);
	if (this->back) this->back->clip_to(other);
}

// Return a list of all polygons in this BSP tree.
std::vector<Polygon> Node::all_polygons() const
{
	std::vector<Polygon> list = this->polygons;
	std::vector<Polygon> list_front, list_back;
	if (this->front) list_front = this->front->all_polygons();
	if (this->back) list_back = this->back->all_polygons();
	list.insert(list.end(), list_front.begin(), list_front.end());
	list.insert(list.end(), list_back.begin(), list_back.end());
	return list;
}

Node * Node::clone() const
{
	Node * ret = new Node();
	ret->polygons = this->polygons;
	ret->plane = this->plane;
	if (this->front) ret->front = this->front->clone();
	if (this->back) ret->back = this->back->clone();
	return ret;
}

// Build a BSP tree out of `polygons`. When called on an existing tree, the
// new polygons are filtered down to the bottom of the tree and become new
// nodes there. Each set of polygons is partitioned using the first polygon
// (no heuristic is used to pick a good split).
void Node::build(const std::vector<Polygon> & list)
{
	if (!list.size()) return;
	if (!this->plane.ok()) this->plane = list[0] .plane;
	std::vector<Polygon> list_front, list_back;
	for (size_t i = 0; i < list.size(); i++) 
	{
		this->plane.split_polygon(list[i], this->polygons, this->polygons, list_front, list_back);
	}
	if (list_front.size()) 
	{
		if (!this->front) this->front = new Node;
		this->front->build(list_front);
	}
	if (list_back.size()) 
	{
		if (!this->back) this->back = new Node;
		this->back->build(list_back);
	}
}

Node::Node() : front(0), back(0)
{
}

Node::Node(const std::vector<Polygon> & list) : front(0), back(0)
{
	build(list);
}

Node::~Node()
{
	delete front;
	delete back;
}

// Public interface implementation

std::vector<Polygon> model_to_polygons(const Model & model)
{
	std::vector<Polygon> list;
	for (size_t i = 0; i < model.indices.size(); i+= 3)
	{
		std::vector<Vertex> triangle;
		for (int j = 0; j < 3; j++)
		{
			Vertex v = model.vertices[model.indices[i + j]];
			triangle.push_back(v);
		}
		list.push_back(Polygon(triangle));
	}
	return list;
}

Model model_from_polygons(const std::vector<Polygon> & polygons)
{
	Model model;
	int p = 0;
	for (size_t i = 0; i < polygons.size(); i++)
	{
		const Polygon & poly = polygons[i];
		for (size_t j = 2; j < poly.vertices.size(); j++)
		{
			model.vertices.push_back(poly.vertices[0]);		model.indices.push_back(p++);
			model.vertices.push_back(poly.vertices[j - 1]);		model.indices.push_back(p++);
			model.vertices.push_back(poly.vertices[j]);		model.indices.push_back(p++);			
		}
	}
	return model;
}

typedef Node * csg_function(const Node * a1, const Node * b1);

inline static Model operation(const Model & a, const Model & b, csg_function fun)
{
	Node * A = new Node(model_to_polygons(a));
	Node * B = new Node(model_to_polygons(b));
	Node * AB = fun(A, B);
	std::vector<Polygon> polygons = AB->all_polygons();
	delete A; A = 0;
	delete B; B = 0;
	delete AB; AB = 0;
	return model_from_polygons(polygons);
}

Model get_union(const Model & a, const Model & b)
{
	return operation(a, b, csg_union);
}

Model get_intersection(const Model & a, const Model & b)
{
	return operation(a, b, csg_intersect);
}

Model get_difference(const Model & a, const Model & b)
{
	return operation(a, b, csg_subtract);
}

}
