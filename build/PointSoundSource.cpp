#include "PointSoundSource.h"

Vector PointSoundSource::sample_ray()
{
	return uniformSampleSphere(randomUniformNextVal(), randomUniformNextVal());
}

PointSoundSource::PointSoundSource()
	: position(Point(0, 0, 0))
{
}

PointSoundSource::PointSoundSource(Point pos)
	: position(pos)
{
}
