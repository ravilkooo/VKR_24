#ifndef SOUNDPATHPOINT_H
#define SOUNDPATHPOINT_H

#include "utils.h"

class SoundPathPoint
{
public:
	Point position;

	double distance;

	virtual ~SoundPathPoint() {}
	
};

#endif