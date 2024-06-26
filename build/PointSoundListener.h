#ifndef POINTSOUNDLISTENER_H
#define POINTSOUNDLISTENER_H

#include "utils.h"
#include "ImpulseResponse.h"
#include "PathData.h"

class PointSoundListener
{
public:
	// позиция источника
	Point position;
	std::vector<double> frequencies;

	size_t x_idx = 0;
	size_t y_idx = 0;

	double sensivity = 1.;

	// пускаем луч в рандомном направлении
	Vector sample_ray(bool is_2d);

	// дошедшие сигналы
	std::vector<PathData> received_diffraction_IRs;
	std::vector<std::pair<double, double>> received_phong_IRs_archived;
	std::vector<std::vector<std::pair<double, double>>> received_diffr_IRs_archived;
	std::vector<std::pair<double, double>> received_direct_IRs_archived;

	// счётчик сигналов каждого типа
	std::vector<unsigned> IR_cnt;

	PointSoundListener();

	PointSoundListener(Point pos);
	PointSoundListener(Point pos, size_t x_idx, size_t y_idx);

	void archiveDiffractionIRs();
	void archivePhongIRs();
	void saveToStream(std::ofstream& out, std::vector<std::pair<double, double>>& received_IRs_archived) const;
	void appendListenerToFile(const std::string& filename);
	void allClear();
};

#endif