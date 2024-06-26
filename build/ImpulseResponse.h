#ifndef IMPULSERESPONSE_H
#define IMPULSERESPONSE_H

#include "utils.h"
#include "SoundDiffractionPathPoint.h"
#include "FrequencyBandResponse.h"

class ImpulseResponse
{
public:
	// ���������� ������
	unsigned num_frequencies = 3;

	// �������� ��������� ��� ������ �������
	FrequencyBandResponse fbr;

	// ����������, ���������� �����
	double distance = 0;

	// ��� ����/����
	IRType ir_type;

	// ����
	std::vector<SoundDiffractionPathPoint> diffraction_points_list;

	// �������� ��� ������ ����� ����
	std::vector<FrequencyBandResponse> point_responses;

	// �������� ����� ��������� ������� ����
	double getTravelingTime(double speed_of_sound);

	size_t last_valid_index;

	ImpulseResponse();

	ImpulseResponse(unsigned num_frequencies);
	
	ImpulseResponse(const ImpulseResponse& IR);
	
	//double getFBR(size_t freq_idx);
};

#endif