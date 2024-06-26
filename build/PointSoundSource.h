#ifndef POINTSOUNDSOURCE_H
#define POINTSOUNDSOURCE_H

#include "utils.h"

class PointSoundSource
{
public:
	// позиция источника
	Point position;

	double power = 1.;
	
	// пускаем луч в рандомном направлении
	Vector sample_ray();

	PointSoundSource();
	
	PointSoundSource(Point pos);
};

#endif