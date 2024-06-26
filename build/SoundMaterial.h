#ifndef SOUNDMATERIAL_H
#define SOUNDMATERIAL_H

#include <string>
#include <boost/math/special_functions/beta.hpp>
#include "utils.h"
#include "FrequencyBandResponse.h"

// ��� ���������� ��������� ����� ����������� ���������� ���������
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

	// ������� �����
	static const SoundMaterial CONCRETE;

	// ����
	static const SoundMaterial SNOW;

	// ��� �������� ��������
	static const SoundMaterial SILENCE;

	// shininess
	static const int n_specular = 100;
};

#endif