#include "SoundPhongPathPoint.h"

SoundPhongPathPoint::SoundPhongPathPoint()
{

}

SoundPhongPathPoint::SoundPhongPathPoint(Point position)
	: position(position)
{

}

bool SoundPhongPathPoint::is_connectible() const
{
	return (refl_type == ReflectionType::BOTH
		|| refl_type == ReflectionType::SPECULAR
		|| refl_type == ReflectionType::SPECULAR);
}

Vector SoundPhongPathPoint::getPhongReflection(Vector in_dir, ReflectionType& refl_type,
	double& _pdf_fwd, double& _pdf_rev,
	double& _brdf_fwd, double& _brdf_rev)
{
	auto pars = material->getClosestParameters(frequency);
	double r = pars.first; // reflectivity
	double s = pars.second; // scaterring
	// ro_d = s*r; ro_s = (1-s)*r
	double ro_d = s * r;
	double ro_s = (1 - s) * r;

	double xi = randomUniformNextVal();

	if (xi > ro_d + ro_s)
	{
		_pdf_fwd = 0;
		_pdf_rev = 0;
		_brdf_fwd = 0;
		_brdf_rev = 0;
		refl_type = ReflectionType::NOCONTRIBUTION;
		return CGAL::Null_vector();
	}

	// std::cout << "ro_d = " << ro_d << "; ro_s = " << ro_s << "\n";

	double t_1, t_2;

	double cos_theta_i = in_dir * normal;
	Vector ideal_spec_refl = 2 * normal * cos_theta_i - in_dir;
	double cos_alpha;
	Vector refl_dir;
	_pdf_fwd = 1;

	double cos_theta_o;
	Vector ideal_spec_refl_rev;
	double cos_alpha_rev;
	_pdf_rev = 1;

	if (xi <= ro_d)
	{
		//diffuse
		refl_type = ReflectionType::DIFFUSE;

		auto rotateToNorm = getRotateFunc(normal);
		do
		{
			t_1 = randomUniformNextVal();
			t_2 = randomUniformNextVal();
			refl_dir = Vector(sqrt(1 - t_1) * cos(2 * Pi * t_2),
				sqrt(1 - t_1) * sin(2 * Pi * t_2),
				sqrt(t_1));
			refl_dir = rotateToNorm(refl_dir);
		} while (refl_dir * normal < PHONG_EPSILON);
	}
	else // if (xi <= ro_d + ro_s)
	{
		//specular
		refl_type = ReflectionType::SPECULAR;

		auto rotToSpecular = getRotateFunc(ideal_spec_refl);
		do
		{
			t_1 = randomUniformNextVal();
			t_2 = randomUniformNextVal();
			refl_dir = Vector(sqrt(1 - pow(t_1, 2. / (material->n_specular + 1))) * cos(2 * Pi * t_2),
				sqrt(1 - pow(t_1, 2. / (material->n_specular + 1))) * sin(2 * Pi * t_2),
				pow(t_1, 1. / (material->n_specular + 1)));
			refl_dir = rotToSpecular(refl_dir);
			if (refl_dir * normal < 0)
			{
				refl_dir = -2 * normal * (refl_dir * normal) + refl_dir;
			}
		} while (refl_dir * normal < PHONG_EPSILON);

	}

	//// _fwd
	//cos_alpha = ideal_spec_refl * refl_dir;
	//cos_alpha = std::max(0., cos_alpha);

	//// _rev
	//cos_theta_o = refl_dir * normal;
	//ideal_spec_refl_rev = 2 * normal * cos_theta_o - refl_dir;
	//cos_alpha_rev = ideal_spec_refl_rev * in_dir;
	//cos_alpha_rev = std::max(0., cos_alpha_rev);

	//// расчёт _brdf_fwd
	//double diffuse_brdf = inv_Pi;
	//double I_M = calc_I_M(cos_alpha, material->n_specular); // коэффициент нормализации
	//double specular_brdf = std::pow(cos_alpha, material->n_specular) / I_M;
	//_brdf_fwd = ro_d* diffuse_brdf + ro_s * specular_brdf;

	//// расчёт brdf_rev
	//double diffuse_brdf_rev = inv_Pi;
	//double I_M_rev = calc_I_M(cos_alpha_rev, material->n_specular); // коэффициент нормализации
	//double specular_brdf_rev = std::pow(cos_alpha_rev, material->n_specular) / I_M_rev;
	//_brdf_rev = ro_d * diffuse_brdf_rev + ro_s * specular_brdf_rev;

	_pdf_fwd = calcPhongPdf(in_dir, refl_dir);
	_brdf_fwd = calcBRDF(in_dir, refl_dir);
	_pdf_rev = calcPhongPdf(refl_dir, in_dir);
	_brdf_rev = calcBRDF(refl_dir, in_dir);
	return refl_dir;
}

double SoundPhongPathPoint::calcPhongPdf(const Vector& in_dir, const Vector& out_dir)
{
	double pdf;

	auto pars = material->getClosestParameters(frequency);
	double r = pars.first; // reflectivity
	double s = pars.second; // scaterring
	// ro_d = s*r; ro_s = (1-s)*r
	double ro_d = s * r;
	double ro_s = (1 - s) * r;

	double cos_theta_i = std::abs(in_dir * normal);

	double pdf_specular = 0;
	double pdf_diffuse = 0;

	// if (refl_type == ReflectionType::SPECULAR || refl_type == ReflectionType::BOTH)
	{
		Vector ideal_spec_refl = 2 * normal * cos_theta_i - in_dir;
		double cos_alpha = std::abs(ideal_spec_refl * out_dir);
		cos_alpha = std::max(0., cos_alpha);

		// расчёт pdf_fwd для specular
		pdf_specular = hemispherePhongSpecularPdf(cos_alpha, material->n_specular);
		pdf_specular *= ro_s; // закоментил, т.к. должно быть затухание / (ro_d + ro_s);
	}
	// закоментил, т.к. одно и тоже направление может быть получено разными типами отражения
	// else if (refl_type == ReflectionType::DIFFUSE || refl_type == ReflectionType::BOTH)
	{
		// расчёт pdf_fwd для diffuse
		pdf_diffuse = hemispherePhongDiffusePdf(cos_theta_i);
		pdf_diffuse *= ro_d; // закоментил, т.к. должно быть затухание / (ro_d + ro_s);
	}
	pdf = pdf_specular + pdf_diffuse;
	assert(pdf > 0);
	return pdf;
}
double SoundPhongPathPoint::calcPhongPdf(const SoundPhongPathPoint& prev, const SoundPhongPathPoint& next)
{
	Vector in_dir = getVectorDirection(prev.position - this->position);
	Vector out_dir = getVectorDirection(next.position - this->position);
	return calcPhongPdf(in_dir, out_dir);	
}

double hemispherePhongDiffusePdf(double cos_theta)
{
	// InvPi * cos(theta)
	return inv_Pi * cos_theta;
}
double hemispherePhongSpecularPdf(double cos_alpha, int n)
{
	// (n+1) Inv2Pi  * cos^n(alpha)
	return (n + 1) * 0.5 * inv_Pi * std::pow(cos_alpha, n);
}

double ibeta(double x, double a, double b)
{
	return boost::math::beta(a, b, x);
}
double gamma_quot(double a, double b)
{
	return std::exp(std::lgamma(a) - std::lgamma(b));
}
double calc_I_M(double NdotV, double n)
{
	double const& costerm = NdotV;
	double        sinterm_sq = 1.0 - costerm * costerm;
	double        halfn = 0.5 * n;

	double negterm = costerm;
	if (n >= 1e-18)
	{
		negterm *= halfn * ibeta(sinterm_sq, halfn, 0.5);
	}

	return (
		2 * Pi * costerm +
		sqrt_Pi * gamma_quot(halfn + 0.5, halfn + 1.0) * (std::pow(sinterm_sq, halfn) - negterm)
		) / (n + 2.0);
}

double SoundPhongPathPoint::geom(const SoundPhongPathPoint& next) const
{
	//std::cout << "geom_in\n";
	auto w = next.position - position;
	double inv_dist2 = 1. / w.squared_length();
	//std::cout << "geom_out\n";
	double geom_val = inv_dist2;
	w *= std::sqrt(inv_dist2);
	if (type == SoundPhongPathPoint::PointType::SURFACE)/*
		|| this->type == VertexType::SphereSound
		|| this->type == VertexType::SphereMicrophone)*/
	{
		geom_val *= std::abs(normal * w);
	}
	if (next.type == SoundPhongPathPoint::PointType::SURFACE)/*
		|| next.type == VertexType::SphereSound
		|| next.type == VertexType::SphereMicrophone)*/
	{
		geom_val *= std::abs(next.normal * w);
	}
	return geom_val;
}

double SoundPhongPathPoint::calcBRDF(const Vector& in_dir, const Vector& out_dir) const
{
	auto pars = material->getClosestParameters(frequency);
	double r = pars.first; // reflectivity
	double s = pars.second; // scaterring
	// ro_d = s*r; ro_s = (1-s)*r
	double ro_d = s * r;
	double ro_s = (1 - s) * r;

	double cos_theta_i = in_dir * normal;

	double pdf_specular = 0;
	double pdf_diffuse = 0;

	double diffuse_brdf = inv_Pi;
	Vector ideal_spec_refl = 2 * normal * cos_theta_i - in_dir;
	double cos_alpha = ideal_spec_refl * out_dir;
	cos_alpha = std::max(0., cos_alpha);
	double I_M = calc_I_M(cos_alpha, material->n_specular); // коэффициент нормализации
	double specular_brdf = std::pow(cos_alpha, material->n_specular) / I_M;

	return ro_d * diffuse_brdf + ro_s * specular_brdf;
}
double SoundPhongPathPoint::calcBRDF(const SoundPhongPathPoint& prev, const SoundPhongPathPoint& next) const
{
	Vector out_dir = getVectorDirection(next.position - this->position);
	Vector in_dir = getVectorDirection(prev.position - this->position);
	return calcBRDF(in_dir, out_dir);
}

double SoundPhongPathPoint::calcPdfClear(const SoundPhongPathPoint* prev, const SoundPhongPathPoint& next)
{
    //std::cout << "pdf_clear_in\n";
    /*if (type == VertexType::SphereSound)
    {
        //std::cout << " || ";
        pdf = CosineHemispherePdf(abs_dot(Normalize(next.p - this->p), this->ng));
        return pdf;
    }
    else if (type == VertexType::SphereMicrophone)
    {
        //std::cout << " || ";
        pdf = CosineHemispherePdf(abs_dot(Normalize(next.p - this->p), this->ng));
        return pdf;
    }
    else*/
	if (type == SoundPhongPathPoint::PointType::SURFACE)
	{
		assert(prev);
		return calcPhongPdf(*prev, next);
	}
    else if (type == SoundPhongPathPoint::PointType::SOURCE)
    {
        return uniformSpherePdf();
    }
    else if (type == SoundPhongPathPoint::PointType::LISTENER)
    {
		return uniformSpherePdf();
    }
	return 0;
}		   
double SoundPhongPathPoint::calcPdfProj(const SoundPhongPathPoint* prev, const SoundPhongPathPoint& next)
{
	double pdf = calcPdfClear(prev, next);
	//std::cout << " (pdf_clear[" << ((type == VertexType::SphereSound ? "SpSn" : "-")) << "] = " << pdf << ") ";
	if (type == SoundPhongPathPoint::PointType::SURFACE)
		//|| type == VertexType::SphereSound
		//|| type == VertexType::SphereMicrophone)
	{
		/*if (abs_dot(this->ng, Normalize(next.p - this->p)) < EPSILON)
		{
			std::cout << "GOTCHA! this->p=" << this->p << ";\t next.p=" << next.p <<"; this->ng=" << this->ng << "\n";
			std::cout << "       dot( , )=" << dot(this->ng, next.p - this->p) << "\n";
		}*/
		pdf /= normal * getVectorDirection(next.position - this->position);
		//std::cout << " (pdf_corr = " << pdf << ") ";
	}
	else if (type == SoundPhongPathPoint::PointType::SOURCE
		|| type == SoundPhongPathPoint::PointType::LISTENER)
	{
		// они напрямую соединяются, потому что не поверхность, а точка
	}
	return pdf;
}
double SoundPhongPathPoint::convertPdf(double pdf, const SoundPhongPathPoint& next) const
{
	auto w = next.position - position;
	double inv_dist2 = 1. / w.squared_length();
	if (next.type == SoundPhongPathPoint::PointType::SURFACE)
	{
		pdf *= std::abs(next.normal * (w * std::sqrt(inv_dist2)));
	}
	return pdf * inv_dist2;
}
