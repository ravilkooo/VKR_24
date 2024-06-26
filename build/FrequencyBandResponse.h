#ifndef FREQBANDRESPONSE_H
#define FREQBANDRESPONSE_H

#include <complex>
typedef std::complex<double> complex_double;
#include "utils.h"
#include "SoundDiffractionPathPoint.h"

class FrequencyBandResponse
{
public:
	unsigned num_frequencies = 3;

	std::vector<double> frequencies;
	std::vector<double> attenuations;

	FrequencyBandResponse();
	
	FrequencyBandResponse(const FrequencyBandResponse& fbr);

	FrequencyBandResponse(unsigned num_frequencies);

	void changeAttenuation(size_t idx, double val);

	void setFrequency(double frequency, double attenuation);
	
	std::vector<double> getFrequencies();
	
	double getClosestValue(double frequency);

	void allClear();

	double operator [] (size_t freq_idx);

	FrequencyBandResponse operator * (const FrequencyBandResponse& other) const;

	FrequencyBandResponse& operator *= (const FrequencyBandResponse& other);

	friend FrequencyBandResponse operator * (double val, const FrequencyBandResponse&);
};

FrequencyBandResponse computeUTDAttenuation(
	const Point& source_position, const SoundDiffractionPathPoint& diffraction_point,
	const Point& listener_position, double speed_of_sound, const std::vector<double>& frequencies,
	double epsilon_h);

FrequencyBandResponse computeUTDAttenuation_2nd_order(
	const Point& source_position, const Point& diffraction_point_1,
	const Point& diffraction_point_2, const Point& listener_position,
	const Vector& source_face_normal_1, const Vector& listener_face_normal_1,
	const Vector& source_face_normal_2, const Vector& listener_face_normal_2,
	const Vector& edge_axis_1, const Vector& edge_axis_2,
	double speed_of_sound, const std::vector<double>& frequencies);

FrequencyBandResponse computeUTDAttenuation_N_order(
	const Point& source_position, std::vector<SoundDiffractionPathPoint> diffraction_points,
	const Point& listener_position,
	double speed_of_sound, const std::vector<double>& frequencies, double epsilon_h);

double UTD_coefficient(double n, double k,
	double p, double r,
	double thetaI, double alphaI, double alphaD);

double UTD_alpha(double beta, double n, int nSign);

double UTD_L(double p, double r, double thetaI);

int UTD_N(double beta, double n, int nSign);

double UTD_cotan(double numer, double denom);

complex_double UTD_euler(double x);

complex_double UTD_estimateF(double X);

complex_double UTD_freqTerm(double n, double k, double thetaI);

double UTD_sphereDisKouyoumjian(double r, double p);

double UTD_sphereDis(double r, double p);

double cotangent(double x);

#endif