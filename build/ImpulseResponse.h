#ifndef IMPULSERESPONSE_H
#define IMPULSERESPONSE_H

#include "utils.h"
#include "SoundDiffractionPathPoint.h"
#include "FrequencyBandResponse.h"

class ImpulseResponse
{
public:
	// количество частот
	unsigned num_frequencies = 3;

	// значения затуханий для каждой частоты
	FrequencyBandResponse fbr;

	// расстояние, прошедешее лучом
	double distance = 0;

	// тип луча/пути
	IRType ir_type;

	// путь
	std::vector<SoundDiffractionPathPoint> diffraction_points_list;

	// значения для каждой точки пути
	std::vector<FrequencyBandResponse> point_responses;

	// получить время получения данного луча
	double getTravelingTime(double speed_of_sound);

	size_t last_valid_index;

	ImpulseResponse();

	ImpulseResponse(unsigned num_frequencies);
	
	ImpulseResponse(const ImpulseResponse& IR);
	
	//double getFBR(size_t freq_idx);
};

#endif