#include "PathData.h"

PathData::PathData() :
	path_type(PathType::SPECULAR),
	energy(FrequencyBandResponse()),
	time(0)
{
}

PathData::PathData(PathType path_type,
	FrequencyBandResponse energy,
	double time,
	std::vector<SoundDiffractionPathPoint> diffraction_points_list) :
	path_type(path_type),
	energy(energy),
	time(time),
	diffraction_points_list(diffraction_points_list)
{
}
