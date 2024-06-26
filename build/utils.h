#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <random>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#define MATH_EPSILON std::numeric_limits<double>::epsilon()
#define MATH_DOUBLE_MAX std::numeric_limits<double>::max()
static const double Pi = 3.14159265358979323846;
static const double inv_Pi = 1 / Pi;
static const double sqrt_Pi = std::sqrt(Pi);

#define DIRECTION_EQUALITY_THRESHOLD 0.99

typedef CGAL::Simple_cartesian<double>          Kernel;
typedef Kernel::Point_3                         Point;
typedef Kernel::Vector_3						Vector;
typedef Kernel::Plane_3							Plane;
typedef CGAL::Surface_mesh<Point>               Mesh;

// typedef Mesh::Vertex_index          vertex_descriptor;
// typedef Mesh::Halfedge_index        hedge_descriptor;
// typedef Mesh::Edge_index            edge_descriptor;
// typedef Mesh::Face_index            face_descriptor;

typedef boost::graph_traits<Mesh>::vertex_descriptor		vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor		halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor          edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor			face_descriptor;

// для нахождения пересечений с мэшем
typedef Kernel::Ray_3 Ray;

template <typename T>
class TemporaryOverride
{
public:
	TemporaryOverride(T* ptr = nullptr, T new_val = T()) : ptr(ptr)
	{
		if (ptr)
		{
			old_val = *ptr;
			*ptr = new_val;
		}
	}

	~TemporaryOverride()
	{
		if (ptr)
		{
			*ptr = old_val;
		}
	}

	TemporaryOverride(const TemporaryOverride&) = delete;
	TemporaryOverride& operator=(const TemporaryOverride&) = delete;

	TemporaryOverride& operator=(TemporaryOverride&& other) noexcept
	{
		if (this != &other)
		{
			restore();
			ptr = other.ptr;
			old_val = other.old_val;
			other.ptr = nullptr;
		}
		return *this;
	}

private:
	T* ptr;
	T old_val;

	void restore()
	{
		if (ptr)
		{
			*ptr = old_val;
		}
	}
};

double getVectorLength(Vector v);
Vector getVectorDirection(Vector v);
Vector getRayDirection(Ray ray);
double degreesToRadians(double degrees);
double mod(double value, double divisor);
double clamp(double number, double minimum, double maximum);

double getSignedDistanceTo(Vector plane_normal, Point plane_point, Point point);
double getSignedDistanceTo(Plane plane, Point plane_point, Point point);
double getSignedDistanceTo(Plane plane, Point point);

double angleBetween(const Vector& v1, const Vector& v2);
Vector projectToPlane(const Vector& v, const Vector& n);
Point rotatePointAroundAxis(const Point& axis_point, const Vector& axis_direction, const Point& point_to_rotate, double angle);

enum class IRType : unsigned
{
	SPECULAR = 0,
	DIFFRACTION = 1,
	DIFFUSION = 2
};

struct Skip
{
	bool no_skip = false;
	face_descriptor fd;
	Skip()
		: fd(Mesh::null_face())
	{
		no_skip = true;
	}
	Skip(const face_descriptor fd)
		: fd(fd)
	{
	}
	bool operator()(const face_descriptor& t) const
	{
		if (no_skip)
		{
			return false;
		}
		return (t == fd);
	}
};

struct SkipFew
{
	bool no_skip = false;
	std::vector<face_descriptor> fd_list;
	SkipFew()
	{
		no_skip = true;
	}
	SkipFew(const face_descriptor fd)
	{
		fd_list.clear();
		fd_list.push_back(fd);
	}
	void add(const face_descriptor fd)
	{
		fd_list.push_back(fd);
	}
	bool operator()(const face_descriptor& t) const
	{
		if (no_skip)
		{
			return false;
		}
		for (face_descriptor fd : fd_list)
		{
			if (t == fd)
				return true;
		}
		return false;
	}
};


double randomUniformNextVal();

Vector uniformSampleSphere(double t1, double t2);

Vector uniformSampleCircle(double t1);

double uniformSpherePdf();

std::function<Vector(Vector)> getRotateFunc(const Vector& dir);

#endif