#include "SoundDiffractionPathPoint.h"

SoundDiffractionPathPoint::SoundDiffractionPathPoint()
{
}

SoundDiffractionPathPoint::SoundDiffractionPathPoint(Point position) : position(position)
{
}

SoundDiffractionPathPoint::SoundDiffractionPathPoint(Point position, BigDiffractionEdge diffr_edge)
	: position(position), diffr_edge(&diffr_edge)
{

}

SoundDiffractionPathPoint::SoundDiffractionPathPoint(Point position, unsigned col)
	: position(position), diffr_edge_col(col)
{
}

SoundDiffractionPathPoint::SoundDiffractionPathPoint(Point position, unsigned col,
	BigDiffractionEdge diffr_edge)
	: position(position), diffr_edge_col(col), diffr_edge(&diffr_edge)
{
}

SoundDiffractionPathPoint::SoundDiffractionPathPoint(Point position, unsigned col,
	BigDiffractionEdge& diffr_edge)
	: position(position), diffr_edge_col(col), diffr_edge(&diffr_edge)
{
}

SoundDiffractionPathPoint::SoundDiffractionPathPoint(Point position, unsigned col,
	BigDiffractionEdge* diffr_edge)
	: position(position), diffr_edge_col(col), diffr_edge(diffr_edge)
{
}

SoundDiffractionPathPoint::SoundDiffractionPathPoint(Point position, unsigned col,
	const BigDiffractionEdge* diffr_edge)
	: position(position), diffr_edge_col(col), diffr_edge(diffr_edge)
{
}