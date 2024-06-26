#ifndef PATHDATA_H
#define PATHDATA_H

#include "utils.h"
#include "FrequencyBandResponse.h"
#include "SoundDiffractionPathPoint.h"

enum class PathType : unsigned
{
	SPECULAR = 0,
	DIFFRACTION = 1,
	DIFFUSION = 2
};

class PathData
{
public:
	PathType path_type;
	FrequencyBandResponse energy;
	double time;
	std::vector<SoundDiffractionPathPoint> diffraction_points_list;

	PathData();

	PathData(PathType path_type,
		FrequencyBandResponse energy,
		double time,
		std::vector<SoundDiffractionPathPoint> diffraction_points_list);

};

#endif // !PATHDATA_H
