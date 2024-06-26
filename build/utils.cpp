#include "utils.h"


double getVectorLength(Vector v)
{
	return std::sqrt(v.squared_length());
}

Vector getVectorDirection(Vector v)
{
	return v / getVectorLength(v);
}

Vector getRayDirection(Ray ray)
{
	return getVectorDirection(ray.to_vector());
}

double degreesToRadians(double degrees)
{
	return degrees * double(0.017453292519943);
}

double mod(double value, double divisor)
{
	return std::fmod(value, divisor);
}

double clamp(double number, double minimum, double maximum)
{
	return std::min(std::max(number, minimum), maximum);
}

double getSignedDistanceTo(Vector plane_normal, Point plane_point, Point point)
{
	return (point - plane_point) * plane_normal;
}

double getSignedDistanceTo(Plane plane, Point plane_point, Point point)
{
	return (point - plane_point) * getVectorDirection(plane.orthogonal_vector());
}

double getSignedDistanceTo(Plane plane, Point point)
{
	Point plane_point = plane.point();
	return (point - plane_point) * getVectorDirection(plane.orthogonal_vector());
}

double angleBetween(const Vector& v1, const Vector& v2)
{
	double angle = getVectorDirection(v1) * getVectorDirection(v2);

	if (angle > 1.)
		return 0;
	else if (angle < -1.)
		return MATH_EPSILON;

	return std::acos(angle);
}

Vector projectToPlane(const Vector& v, const Vector& n)
{
	double projAmount = v * n;

	return v - projAmount * n;
}


double randomUniformNextVal()
{
	static std::random_device rd;
	static std::mt19937 generator = std::mt19937(rd());
	static std::uniform_real_distribution<> distribution = std::uniform_real_distribution<>(0., 1.);
	return distribution(generator);
}

Vector uniformSampleSphere(double t1, double t2)
{
	double z = 1 - 2 * t1;
	double r = std::sqrt(std::max(0., 1. - z * z));
	double phi = 2 * Pi * t2;
	return Vector(r * std::cos(phi), r * std::sin(phi), z);
};

Vector uniformSampleCircle(double t1)
{
	double phi = 2 * Pi * t1;
	return Vector(std::cos(phi), std::sin(phi), 0);
};

double uniformSpherePdf()
{
	// Inv4Pi;
	return 0.25 * inv_Pi;
}

std::function<Vector(Vector)> getRotateFunc(const Vector& dir)
{
	double cos_theta = dir[2];
	double sin_theta = std::sqrt(1 - dir[2] * dir[2]);
	double cos_phi = sin_theta < 1e-5 ? 1 : dir[0] / sin_theta;
	double sin_phi = sin_theta < 1e-5 ? 0 : dir[1] / sin_theta;
	std::array<std::array<double, 3>, 3> matrix{ {
		{cos_theta * cos_phi, -sin_phi, sin_theta * cos_phi},
		{cos_theta * sin_phi, cos_phi, sin_theta * sin_phi},
		{-sin_theta, 0, cos_theta}
		} };
	auto rot = [matrix](Vector x)->Vector {
		std::array<double, 3> result{ {0, 0, 0} };
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				result[i] += matrix[i][j] * x[j];
			}
		}
		return Vector(result[0], result[1], result[2]);
	};
	return rot;
}

Point rotatePointAroundAxis(const Point& axis_point, const Vector& axis_direction, const Point& point_to_rotate, double angle)
{
	// Нормируем направление оси
	Vector norm_axis_direction = axis_direction / std::sqrt(axis_direction.squared_length());

	// Создаем матрицу поворота
	double c = std::cos(angle);
	double s = std::sin(angle);
	double t = 1 - c;
	double x = norm_axis_direction.x();
	double y = norm_axis_direction.y();
	double z = norm_axis_direction.z();

	double rotation_matrix[3][3] = {
		{t * x * x + c,    t * x * y - s * z,  t * x * z + s * y},
		{t * x * y + s * z,  t * y * y + c,    t * y * z - s * x},
		{t * x * z - s * y,  t * y * z + s * x,  t * z * z + c}
	};

	// Вычитаем точку оси, чтобы перенести её в начало координат
	Vector translated_point = point_to_rotate - CGAL::ORIGIN - (axis_point - CGAL::ORIGIN);

	// Поворачиваем точку
	Vector rotated_vector(
		rotation_matrix[0][0] * translated_point.x() + rotation_matrix[0][1] * translated_point.y() + rotation_matrix[0][2] * translated_point.z(),
		rotation_matrix[1][0] * translated_point.x() + rotation_matrix[1][1] * translated_point.y() + rotation_matrix[1][2] * translated_point.z(),
		rotation_matrix[2][0] * translated_point.x() + rotation_matrix[2][1] * translated_point.y() + rotation_matrix[2][2] * translated_point.z()
	);

	// Переносим точку обратно в исходную позицию
	Point rotated_point = CGAL::ORIGIN + rotated_vector + (axis_point - CGAL::ORIGIN);

	return rotated_point;
}