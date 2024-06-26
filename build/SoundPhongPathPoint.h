#ifndef SOUNDPHONGPATHPOINT_H
#define SOUNDPHONGPATHPOINT_H

#include "utils.h"
#include "SoundPathPoint.h"
#include "SoundMaterial.h"

class SoundPhongPathPoint : public SoundPathPoint
{
public:

	enum class PointType : unsigned
	{
		UNDEFINED = 0,
		// точечный источник
		SOURCE = 1,
		// точечный приёмник
		LISTENER = 2,
		SURFACE = 3
	};

	enum class ReflectionType : unsigned
	{
		DIFFUSE = 0,
		SPECULAR = 1,
		BOTH = 2,
		NOCONTRIBUTION = 3
	};

	PointType type;
	Point position;
	Vector normal = CGAL::Null_vector();
	SoundMaterial* material;

	double distance;
	double frequency;
	double attenuation = 1;
	double beta = 1;

	ReflectionType refl_type;
	double brdf_fwd = 1;
	double brdf_rev = 1;
	double pdf_fwd = 1;
	double pdf_rev = 1;

	SoundPhongPathPoint();
	SoundPhongPathPoint(Point position);

	bool is_connectible() const;

	Vector getPhongReflection(Vector in_dir, ReflectionType& refl_type,
		double& _pdf_fwd, double& _pdf_rev,
		double& _brdf, double& _brdf_rev);
	double calcPhongPdf(const Vector& in_dir, const Vector& out_dir);
	double calcPhongPdf(const SoundPhongPathPoint& prev, const SoundPhongPathPoint& next);

	double geom(const SoundPhongPathPoint& next) const;

	double calcBRDF(const Vector& in_dir, const Vector& out_dir) const;
	double calcBRDF(const SoundPhongPathPoint& prev, const SoundPhongPathPoint& next) const;

	double calcPdfProj(const SoundPhongPathPoint* prev, const SoundPhongPathPoint& next);
	double calcPdfClear(const SoundPhongPathPoint* prev, const SoundPhongPathPoint& next);
	double convertPdf(double pdf, const SoundPhongPathPoint& next) const;

};

double hemispherePhongDiffusePdf(double cos_theta);

double hemispherePhongSpecularPdf(double cos_alpha, int n);

double ibeta(double x, double a, double b);

double gamma_quot(double a, double b);

double calc_I_M(double NdotV, double n);

#endif