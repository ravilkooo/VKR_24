#ifndef POINTSOUNDSOURCE_H
#define POINTSOUNDSOURCE_H

#include "utils.h"

class PointSoundSource
{
public:
	// ������� ���������
	Point position;

	double power = 1.;
	
	// ������� ��� � ��������� �����������
	Vector sample_ray();

	PointSoundSource();
	
	PointSoundSource(Point pos);
};

#endif