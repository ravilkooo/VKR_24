#include "BigDiffractionEdge.h"

double BigDiffractionEdge::getLength() const
{
	return getVectorLength(m.point(vd_first) - m.point(vd_last));
}

Vector BigDiffractionEdge::getAxis() const
{
	return getVectorDirection(m.point(vd_last) - m.point(vd_first));
}

Point BigDiffractionEdge::getStart() const
{
	return m.point(vd_first);
}

bool BigDiffractionEdge::testOrientation(const Point point, double offset) const
{
	double d1 = getSignedDistanceTo(normal1, m.point(vd_first), point);
	double d2 = getSignedDistanceTo(normal2, m.point(vd_first), point);

	bool front1 = d1 > offset;
	bool front2 = d2 > offset;
	bool behind1 = d1 < -offset;
	bool behind2 = d2 < -offset;

	return !((front1 && front2) || (behind1 && behind2));
}

bool BigDiffractionEdge::testEdgeOrientation(const BigDiffractionEdge& edge2, double offset = 0.0001)
{
	double e11_distance_1 = getSignedDistanceTo(edge2.normal1, m.point(edge2.vd_first), m.point(vd_first));
	double e11_distance_2 = getSignedDistanceTo(edge2.normal2, m.point(edge2.vd_first), m.point(vd_first));
	double e12_distance_1 = getSignedDistanceTo(edge2.normal1, m.point(edge2.vd_first), m.point(vd_last));
	double e12_distance_2 = getSignedDistanceTo(edge2.normal2, m.point(edge2.vd_first), m.point(vd_last));

	bool e11_front_1 = e11_distance_1 > offset;
	bool e11_front_2 = e11_distance_2 > offset;
	bool e12_front_1 = e12_distance_1 > offset;
	bool e12_front_2 = e12_distance_2 > offset;

	bool e11_behind_1 = e11_distance_1 < -offset;
	bool e11_behind_2 = e11_distance_2 < -offset;
	bool e12_behind_1 = e12_distance_1 < -offset;
	bool e12_behind_2 = e12_distance_2 < -offset;

	bool e1Outside = (e11_front_1 && e11_front_2 && e12_front_1 && e12_front_2) ||
		(e11_behind_1 && e11_behind_2 && e12_behind_1 && e12_behind_2);

	if (e1Outside)
		return false;

	double e21_distance_1 = getSignedDistanceTo(normal1, m.point(vd_first), m.point(edge2.vd_first));
	double e21_distance_2 = getSignedDistanceTo(normal2, m.point(vd_first), m.point(edge2.vd_first));
	double e22_distance_1 = getSignedDistanceTo(normal1, m.point(vd_first), m.point(edge2.vd_last));
	double e22_distance_2 = getSignedDistanceTo(normal2, m.point(vd_first), m.point(edge2.vd_last));

	bool e21_front_1 = e21_distance_1 > offset;
	bool e21_front_2 = e21_distance_2 > offset;
	bool e22_front_1 = e22_distance_1 > offset;
	bool e22_front_2 = e22_distance_2 > offset;

	bool e21_behind_1 = e21_distance_1 < -offset;
	bool e21_behind_2 = e21_distance_2 < -offset;
	bool e22_behind_1 = e22_distance_1 < -offset;
	bool e22_behind_2 = e22_distance_2 < -offset;

	bool e2Outside = (e21_front_1 && e21_front_2 && e22_front_1 && e22_front_2) ||
		(e21_behind_1 && e21_behind_2 && e22_behind_1 && e22_behind_2);

	if (e2Outside)
		return false;

	return true;
}

Point BigDiffractionEdge::getFreeVertex1() const
{
	return m.point(m.target(m.next(hed_first)));

}

Point BigDiffractionEdge::getFreeVertex2() const
{
	return m.point(m.target(m.next(hed_last)));
}

Plane BigDiffractionEdge::getPlane1() const
{
	return Plane(getStart(), normal1);
}

Plane BigDiffractionEdge::getPlane2() const
{
	return Plane(getStart(), normal2);
}
