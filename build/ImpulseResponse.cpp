#include "ImpulseResponse.h"

double ImpulseResponse::getTravelingTime(double speed_of_sound)
{
	return distance / speed_of_sound;
}

ImpulseResponse::ImpulseResponse()
{
	diffraction_points_list = std::vector<SoundDiffractionPathPoint>();
	point_responses = std::vector<FrequencyBandResponse>();
}

ImpulseResponse::ImpulseResponse(unsigned num_frequencies)
	: num_frequencies(num_frequencies),
	fbr(num_frequencies)
{
	diffraction_points_list = std::vector<SoundDiffractionPathPoint>();
	point_responses = std::vector<FrequencyBandResponse>();
}

ImpulseResponse::ImpulseResponse(const ImpulseResponse& IR)
	: num_frequencies(IR.num_frequencies),
	fbr(IR.fbr),
	distance(IR.distance),
	ir_type(IR.ir_type),
	diffraction_points_list(IR.diffraction_points_list)
{
	// Копирование объектов из diffraction_points_list
	for (auto point : IR.diffraction_points_list)
	{
		diffraction_points_list.push_back(point);
	}
}
