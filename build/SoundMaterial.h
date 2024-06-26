#ifndef SOUNDMATERIAL_H
#define SOUNDMATERIAL_H

#include <string>
#include <boost/math/special_functions/beta.hpp>
#include "utils.h"
#include "FrequencyBandResponse.h"

// Для отсеивания отражений почти паралельных отражающей плоскости
#define PHONG_EPSILON 0.01

class SoundMaterial
{
public:

	std::string material_name;
	FrequencyBandResponse reflectivity;
	FrequencyBandResponse scattering;

	SoundMaterial();
	SoundMaterial(FrequencyBandResponse reflectivity, FrequencyBandResponse scattering, std::string material_name);

	std::pair<double, double> getClosestParameters(double frequency);

	// гладкий бетон
	static const SoundMaterial CONCRETE;

	// снег
	static const SoundMaterial SNOW;

	// мой глушащий материал
	static const SoundMaterial SILENCE;

	// shininess
	static const int n_specular = 100;
};

#endif