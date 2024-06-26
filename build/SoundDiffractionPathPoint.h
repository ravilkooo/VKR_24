#ifndef SOUNDDIFFRACTIONPATHPOINT_H
#define SOUNDDIFFRACTIONPATHPOINT_H

#include "utils.h"
#include "SoundPathPoint.h"
#include "BigDiffractionEdge.h"

class SoundDiffractionPathPoint : public SoundPathPoint
{
public:
	Point position;

	unsigned diffr_edge_col = 0;

	double distance = 0.;

	Plane listener_plane;

	Plane source_plane;

	const BigDiffractionEdge* diffr_edge;

	SoundDiffractionPathPoint();

	SoundDiffractionPathPoint(Point position);

	SoundDiffractionPathPoint(Point position, BigDiffractionEdge diffr_edge);

	SoundDiffractionPathPoint(Point position, unsigned col);

	SoundDiffractionPathPoint(Point position, unsigned col, BigDiffractionEdge diffr_edge);
	SoundDiffractionPathPoint(Point position, unsigned col, BigDiffractionEdge& diffr_edge);
	SoundDiffractionPathPoint(Point position, unsigned col, BigDiffractionEdge* diffr_edge);
	SoundDiffractionPathPoint(Point position, unsigned col, const BigDiffractionEdge* diffr_edge);
};

#endif