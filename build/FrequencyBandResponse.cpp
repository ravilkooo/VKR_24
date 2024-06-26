#include "FrequencyBandResponse.h"

FrequencyBandResponse::FrequencyBandResponse() :
	attenuations(std::vector<double>(num_frequencies, 1.))
{
	const double min_freq = 20.0f;
	const double max_freq = 20000.0f;
	const double ln_min_freq = std::log(min_freq);
	const double ln_max_freq = std::log(max_freq);

	frequencies = std::vector<double>(num_frequencies);

	for (int i = 0; i < num_frequencies; i++)
	{
		const double ln_freq = ln_min_freq + (i + 0.5) * (ln_max_freq - ln_min_freq) / num_frequencies;
		frequencies[i] = std::exp(ln_freq);
	}
}

FrequencyBandResponse::FrequencyBandResponse(const FrequencyBandResponse& fbr):
	num_frequencies(fbr.num_frequencies),
	attenuations(fbr.attenuations),
	frequencies(fbr.frequencies)
{
}

FrequencyBandResponse::FrequencyBandResponse(unsigned num_frequencies) :
	num_frequencies(num_frequencies),
	attenuations(std::vector<double>(num_frequencies, 1.))
{
	if (num_frequencies == 1)
	{
		frequencies = std::vector<double>(num_frequencies);
		frequencies[0] = 100;
	}
	else if (num_frequencies == 3)
	{
		frequencies = std::vector<double>(num_frequencies);
		frequencies[0] = 84;
		frequencies[1] = 369;
		frequencies[2] = 1625;
	}
	else
	{
		const double min_freq = 20.;
		const double max_freq = 20000.;
		const double ln_min_freq = std::log(min_freq);
		const double ln_max_freq = std::log(max_freq);

		frequencies = std::vector<double>(num_frequencies);

		for (int i = 0; i < num_frequencies; i++)
		{
			const double ln_freq = ln_min_freq + (i + 0.5) * (ln_max_freq - ln_min_freq) / num_frequencies;
			frequencies[i] = std::exp(ln_freq);
		}

	}

}

void FrequencyBandResponse::changeAttenuation(size_t idx, double val)
{
	attenuations[idx] = val;
}

double FrequencyBandResponse::operator[](size_t freq_idx)
{
	return attenuations[freq_idx];
}

FrequencyBandResponse FrequencyBandResponse::operator * (const FrequencyBandResponse & other) const
{
	FrequencyBandResponse new_fbr = FrequencyBandResponse(other.num_frequencies);
	for (int i = 0; i < num_frequencies; i++)
	{
		new_fbr.attenuations[i] = attenuations[i] * other.attenuations[i];
	}
	return new_fbr;
}

FrequencyBandResponse& FrequencyBandResponse::operator *= (const FrequencyBandResponse& other)
{
	for (int i = 0; i < num_frequencies; i++)
	{
		attenuations[i] *= other.attenuations[i];
	}
	return *this;
}

FrequencyBandResponse operator*(double val, const FrequencyBandResponse& response)
{
	FrequencyBandResponse new_fbr = FrequencyBandResponse(response.num_frequencies);
	for (int i = 0; i < response.num_frequencies; i++)
	{
		new_fbr.attenuations[i] = response.attenuations[i] * val;
	}
	return new_fbr;
}


void FrequencyBandResponse::setFrequency(double frequency, double attenuation)
{
	if (num_frequencies == 0 || frequency > frequencies.back())
	{
		frequencies.push_back(frequency);
		attenuations.push_back(attenuation);
		num_frequencies = frequencies.size();
		return;
	}

	// -------freq1]-------freq2]--...-------freqN]
	for (unsigned i = 0; i < num_frequencies; i++)
	{
		if (frequency < frequencies[i])
		{
			frequencies.insert(frequencies.begin() + i, frequency);
			attenuations.insert(attenuations.begin() + i, attenuation);
			num_frequencies = frequencies.size();
			return;
		}
		else if (frequency == frequencies[i])
		{
			attenuations[i] = attenuation;
			return;
		}
	}
}

std::vector<double> FrequencyBandResponse::getFrequencies()
{
	return frequencies;
}

void FrequencyBandResponse::allClear()
{
	num_frequencies = 0;
	frequencies.clear();
	attenuations.clear();
}

double FrequencyBandResponse::getClosestValue(double frequency)
{
	if (frequency > frequencies.back())
	{
		return attenuations.back();
	}

	// -------freq1]-------freq2]--...-------freqN]------(freqINF=freqN)
	for (unsigned i = 0; i < frequencies.size(); i++)
	{
		if (frequency <= frequencies[i])
		{
			return attenuations[i];
		}
	}
	return attenuations.back();
}

FrequencyBandResponse computeUTDAttenuation(
	const Point& source_position, const SoundDiffractionPathPoint& diffraction_point,
	const Point& listener_position, double speed_of_sound, const std::vector<double>& frequencies,
	double epsilon_h)
{
	Vector s_f_n_curr = getVectorDirection(diffraction_point.source_plane.orthogonal_vector());
	Vector l_f_n_curr = getVectorDirection(diffraction_point.listener_plane.orthogonal_vector());

	Point _on_curr_edge = diffraction_point.diffr_edge->getStart();
	/*double _signed_dist = getSignedDistanceTo(Plane(_on_curr_edge, l_f_n_curr), _on_curr_edge,
		listener_position);
	if (_signed_dist < -epsilon_h)
	{
		s_f_n_curr = diffraction_point.diffr_edge->normal2;
		l_f_n_curr = diffraction_point.diffr_edge->normal1;
	}*/
	Vector edge_axis_curr = diffraction_point.diffr_edge->getAxis();
	if (CGAL::cross_product(s_f_n_curr, l_f_n_curr) * edge_axis_curr < 0)
	{
		edge_axis_curr = -diffraction_point.diffr_edge->getAxis();
	}

	//std::cout << "*********** edge_axis source_face_normal " << edge_axis << "; " << source_face_normal << ";\n";
	const Vector source_face_vector = CGAL::cross_product(edge_axis_curr, s_f_n_curr);

	const double n = (2. * Pi - angleBetween(-s_f_n_curr, l_f_n_curr)) / Pi;
	Vector source_direction = source_position - diffraction_point.position;  //from apex point to source
	Vector listener_direction = listener_position - diffraction_point.position; // from apex to receiver
	const double p = getVectorLength(source_direction);
	const double r = getVectorLength(listener_direction);

	if (p < MATH_EPSILON || r < MATH_EPSILON)
		return FrequencyBandResponse(frequencies.size());

	source_direction /= p;
	listener_direction /= r;

	double thetaI = angleBetween(source_direction, edge_axis_curr);

	if (thetaI > Pi * 0.5)
		thetaI = Pi - thetaI;

	const Vector s_dir = projectToPlane(source_direction, edge_axis_curr);
	const Vector r_dir = projectToPlane(listener_direction, edge_axis_curr);

	//std::cout << "*********** source_direction listener_direction " << source_direction << "; " << listener_direction << "; " << ";\n";
	//std::cout << "*********** s_dir r_dir " << s_dir << "; " << r_dir << "; " << source_face_vector << "; " << n << ";\n";

	const double alpha_i = angleBetween(-s_dir, source_face_vector);
	const double alpha_d = angleBetween(r_dir, source_face_vector) + Pi;
	const double alpha_SB = alpha_i + Pi + 0.001;
	const double lerp = clamp((n * Pi - alpha_d) / (n * Pi - alpha_SB), 0, 1);

	//std::cout << "*********** alpha " << alpha_i << " " << alpha_d << " " << alpha_SB << ";\n";

	FrequencyBandResponse result(frequencies.size());

	for (int i = 0; i < frequencies.size(); i++)
	{
		double lambda = speed_of_sound / frequencies[i];
		double k = 2. * Pi / lambda;

		double UTD_coeff = UTD_coefficient(n, k, p, r, thetaI, alpha_i, alpha_d);
		double SB_coeff = UTD_coefficient(n, k, p, r, thetaI, alpha_i, alpha_SB);
		double norm_coeff = clamp(UTD_coeff / SB_coeff, 0, 1);
		double final_coeff = (1.0f - lerp) * UTD_coeff + (lerp)*norm_coeff;
		/*std::cout << "*********** UTD_coeff " << UTD_coeff << ";\n";
		std::cout << "*********** SB_coeff " << SB_coeff << ";\n";
		std::cout << "*********** final_coeff " << final_coeff << ";\n";*/
		// Square to convert to intensity from pressure.
		result.changeAttenuation(i, clamp(final_coeff * final_coeff, 0., 1.));
		/*if (i == 0)
		{
			std::cout << "1st order: (" << listener_position.x() << ", " << listener_position.y() << ", " << listener_position.z() << ")\n\t";
			std::cout << UTD_coeff << "; " << SB_coeff << "; " << norm_coeff << "; " << final_coeff << " || ";
			std::cout << n << "; " << p << "; " << r << "; " << alpha_i << "; " << alpha_SB << "; " << alpha_d << "\n";
			if (listener_position.x() == 2.2 && listener_position.y() == 2.8 && listener_position.z() == 0)
			{
				std::cout << "\n";
			}
		}*/
		// ????
		//result[i] = std::clamp( norm_coeff*norm_coeff, Real(0), Real(1) );
	}

	return result;
}

FrequencyBandResponse computeUTDAttenuation_2nd_order(
	const Point& source_position, const Point& diffraction_point_1,
	const Point& diffraction_point_2, const Point& listener_position,
	const Vector& source_face_normal_1, const Vector& listener_face_normal_1,
	const Vector& source_face_normal_2, const Vector& listener_face_normal_2,
	const Vector& edge_axis_1, const Vector& edge_axis_2,
	double speed_of_sound, const std::vector<double>& frequencies)
{
	/*			..
			 .. SB_1
		 ..	D1. . . . .	D2 . . . . L_imaginary
	   ..  /  \		   / | .SB_2
	 S	  /    \	  /  |	 .	
		 /      \	 /   |	   L
	
	*/

	//std::cout << "*********** edge_axis source_face_normal " << edge_axis << "; " << source_face_normal << ";\n";
	const Vector source_face_vector_2 = CGAL::cross_product(edge_axis_2, source_face_normal_2);

	const double n_2 = (2. * Pi - angleBetween(-source_face_normal_2, listener_face_normal_2)) / Pi;
	Vector source_direction_2 = diffraction_point_1 - diffraction_point_2;  //from apex point to source
	Vector listener_direction_2 = listener_position - diffraction_point_2; // from apex to receiver
	const double p_2 = getVectorLength(source_direction_2);
	const double r_2 = getVectorLength(listener_direction_2);

	if (p_2 < MATH_EPSILON || r_2 < MATH_EPSILON)
		return FrequencyBandResponse(frequencies.size());

	source_direction_2 /= p_2;
	listener_direction_2 /= r_2;

	double theta_I_2 = angleBetween(source_direction_2, edge_axis_2);

	if (theta_I_2 > Pi * 0.5)
		theta_I_2 = Pi - theta_I_2;

	const Vector s_dir_2 = projectToPlane(source_direction_2, edge_axis_2);
	const Vector r_dir_2 = projectToPlane(listener_direction_2, edge_axis_2);

	//std::cout << "*********** source_direction listener_direction " << source_direction << "; " << listener_direction << "; " << ";\n";
	//std::cout << "*********** s_dir r_dir " << s_dir << "; " << r_dir << "; " << source_face_vector << "; " << n << ";\n";

	const double alpha_i_2 = angleBetween(-s_dir_2, source_face_vector_2);
	const double alpha_d_2 = angleBetween(r_dir_2, source_face_vector_2) + Pi;
	const double alpha_SB_2 = alpha_i_2 + Pi + 0.001;

	const double lerp_2 = (n_2 * Pi - alpha_d_2) / (n_2 * Pi - alpha_SB_2);

	// ....................
	/*			..
			 .. SB_1 (не угол, а плоскость)
		 ..	D1. . . . .	D2 . . . . L_imaginary
	   ..  /  \		   / | . SB_2 (не угол, а плоскость)
	 S	  /    \	  /  |	 .
		 /      \	 /   |	   L

	*/
	Point shadow_listener_position;
	if (CGAL::cross_product(diffraction_point_1 - diffraction_point_2, diffraction_point_2 - listener_position) * edge_axis_2 > 0)
	{
		shadow_listener_position = rotatePointAroundAxis(diffraction_point_2, edge_axis_2, listener_position, alpha_SB_2 - alpha_d_2);
	}
	else
	{
		shadow_listener_position = rotatePointAroundAxis(diffraction_point_2, -edge_axis_2, listener_position, alpha_SB_2 - alpha_d_2);
	}
	//std::cout << "*********** edge_axis source_face_normal " << edge_axis << "; " << source_face_normal << ";\n";
	const Vector source_face_vector_1 = CGAL::cross_product(edge_axis_1, source_face_normal_1);

	const double n_1 = (2. * Pi - angleBetween(-source_face_normal_1, listener_face_normal_1)) / Pi;
	Vector source_direction_1 = source_position - diffraction_point_1;  //from apex point to source
	Vector listener_direction_1 = shadow_listener_position - diffraction_point_1; // from apex to receiver
	const double p_1 = getVectorLength(source_direction_1);
	const double r_1 = getVectorLength(listener_direction_1);

	if (p_1 < MATH_EPSILON || r_1 < MATH_EPSILON)
		return FrequencyBandResponse(frequencies.size());

	source_direction_1 /= p_1;
	listener_direction_1 /= r_1;

	double theta_I_1 = angleBetween(source_direction_1, edge_axis_1);

	if (theta_I_1 > Pi * 0.5)
		theta_I_1 = Pi - theta_I_1;

	const Vector s_dir_1 = projectToPlane(source_direction_1, edge_axis_1);
	const Vector r_dir_1 = projectToPlane(listener_direction_1, edge_axis_1);

	//std::cout << "*********** source_direction listener_direction " << source_direction << "; " << listener_direction << "; " << ";\n";
	//std::cout << "*********** s_dir r_dir " << s_dir << "; " << r_dir << "; " << source_face_vector << "; " << n << ";\n";

	const double alpha_i_1 = angleBetween(-s_dir_1, source_face_vector_1);
	const double alpha_d_1 = angleBetween(r_dir_1, source_face_vector_1) + Pi;
	const double alpha_SB_1 = alpha_i_1 + Pi + 0.001;

	const double lerp_1 = (n_1 * Pi - alpha_d_1) / (n_1 * Pi - alpha_SB_1);

	//std::cout << "*********** alpha " << alpha_i << " " << alpha_d << " " << alpha_SB << ";\n";

	FrequencyBandResponse result(frequencies.size());

	for (int i = 0; i < frequencies.size(); i++)
	{
		double lambda = speed_of_sound / frequencies[i];
		double k = 2. * Pi / lambda;

		double UTD_1_incident_coeff = UTD_coefficient(n_1, k, p_1, r_1, theta_I_1, alpha_i_1, alpha_d_1);
		double SB_1_incident_coeff = UTD_coefficient(n_1, k, p_1, r_1, theta_I_1, alpha_i_1, alpha_SB_1);
		double norm_coeff_1 = UTD_1_incident_coeff / SB_1_incident_coeff;
		double final_coeff_1 = (1.0f - lerp_1) * UTD_1_incident_coeff + (lerp_1)*norm_coeff_1;

		double UTD_2_coeff = UTD_coefficient(n_2, k, p_2, r_2, theta_I_2, alpha_i_2, alpha_d_2);
		double SB_2_coeff = UTD_coefficient(n_2, k, p_2, r_2, theta_I_2, alpha_i_2, alpha_SB_2);
		double norm_coeff_2 = UTD_2_coeff / SB_2_coeff;
		double final_coeff_2 = (1.0f - lerp_2) * UTD_2_coeff + (lerp_2)*norm_coeff_2;


		//double final_coeff = (1.0f - lerp) * UTD_coeff + (lerp)*norm_coeff;
		//double final_coeff = SB_2_coeff > 0 ? UTD_1_incident_coeff * UTD_2_coeff / SB_2_coeff : 0;
		//double final_coeff = final_coeff_1 * final_coeff_2;
		/*std::cout << "*********** UTD_coeff " << UTD_coeff << ";\n";
		std::cout << "*********** SB_coeff " << SB_coeff << ";\n";
		std::cout << "*********** final_coeff " << final_coeff << ";\n";*/
		// Square to convert to intensity from pressure.
		result.changeAttenuation(i, clamp(final_coeff_1 * final_coeff_1, 0., 1.) * clamp(final_coeff_2 * final_coeff_2, 0., 1.));
		if (i == 0)
		{
			//std::cout << "first_corr: " << UTD_2_coeff << "; " << SB_2_coeff << "; " << norm_coeff_2 << "; " << final_coeff_2 << '\n';
			//std::cout << "last_corr: " << UTD_1_incident_coeff << "; " << SB_1_incident_coeff << "; " << norm_coeff_1 << "; " << final_coeff_1 << '\n';
		}
		//result[i] = std::clamp( norm_coeff*norm_coeff, Real(0), Real(1) );
	}

	return result;
}


FrequencyBandResponse computeUTDAttenuation_N_order(
	const Point& source_position, std::vector<SoundDiffractionPathPoint> diffraction_points,
	const Point& listener_position,
	double speed_of_sound, const std::vector<double>& frequencies, double epsilon_h)
{
	std::vector<double> total_attenuations(frequencies.size());
	Point shadow_listener_position = listener_position;
	/*			..
			 .. SB_1
		 ..	D1. . . . .	D2 . . . . L_imaginary
	   ..  /  \		   / | .SB_2
	 S	  /    \	  /  |	 .
		 /      \	 /   |	   L

	*/
	// первый коэф затухания
	{

		SoundDiffractionPathPoint& diffr_point_curr = diffraction_points[1];
		SoundDiffractionPathPoint& diffr_point_prevcurr = diffraction_points[2];

		Vector s_f_n_curr = getVectorDirection(diffr_point_curr.source_plane.orthogonal_vector());
		Vector l_f_n_curr = getVectorDirection(diffr_point_curr.listener_plane.orthogonal_vector());

		Point _on_curr_edge = diffr_point_curr.diffr_edge->getStart();
		/*double _signed_dist = getSignedDistanceTo(Plane(_on_curr_edge, l_f_n_curr), _on_curr_edge,
			listener_position);
		if (_signed_dist < -epsilon_h)
		{
			s_f_n_curr = diffr_point_curr.diffr_edge->normal2;
			l_f_n_curr = diffr_point_curr.diffr_edge->normal1;
		}*/
		Vector edge_axis_curr = diffr_point_curr.diffr_edge->getAxis();
		if (CGAL::cross_product(s_f_n_curr, l_f_n_curr) * edge_axis_curr < 0)
		{
			edge_axis_curr = -diffr_point_curr.diffr_edge->getAxis();
		}


		//std::cout << "*********** edge_axis source_face_normal " << edge_axis << "; " << source_face_normal << ";\n";
		const Vector source_face_vector_curr = CGAL::cross_product(edge_axis_curr, s_f_n_curr);

		const double n_curr = (2. * Pi - angleBetween(-s_f_n_curr, l_f_n_curr)) / Pi;
		Vector source_direction_curr = diffr_point_prevcurr.position - diffr_point_curr.position;  //from apex point to source
		Vector listener_direction_curr = listener_position - diffr_point_curr.position; // from apex to receiver
		const double p_curr = getVectorLength(source_direction_curr);
		const double r_curr = getVectorLength(listener_direction_curr);

		if (p_curr < MATH_EPSILON || r_curr < MATH_EPSILON)
		{
			for (int i = 0; i < frequencies.size(); i++)
			{
				total_attenuations[i] = 1;
			}
		}
		else
		{
			source_direction_curr /= p_curr;
			listener_direction_curr /= r_curr;

			double theta_I_curr = angleBetween(source_direction_curr, edge_axis_curr);

			if (theta_I_curr > Pi * 0.5)
				theta_I_curr = Pi - theta_I_curr;

			const Vector s_dir_curr = projectToPlane(source_direction_curr, edge_axis_curr);
			const Vector r_dir_curr = projectToPlane(listener_direction_curr, edge_axis_curr);

			//std::cout << "*********** source_direction listener_direction " << source_direction << "; " << listener_direction << "; " << ";\n";
			//std::cout << "*********** s_dir r_dir " << s_dir << "; " << r_dir << "; " << source_face_vector << "; " << n << ";\n";

			const double alpha_i_curr = angleBetween(-s_dir_curr, source_face_vector_curr);
			const double alpha_d_curr = angleBetween(r_dir_curr, source_face_vector_curr) + Pi;
			const double alpha_SB_curr = alpha_i_curr + Pi + 0.001;

			const double lerp_curr = clamp((n_curr * Pi - alpha_d_curr) / (n_curr * Pi - alpha_SB_curr), 0, 1);

			for (int i = 0; i < frequencies.size(); i++)
			{
				double lambda = speed_of_sound / frequencies[i];
				double k = 2. * Pi / lambda;

				double UTD_coeff_curr = UTD_coefficient(n_curr, k, p_curr, r_curr, theta_I_curr, alpha_i_curr, alpha_d_curr);
				double SB_coeff_curr = UTD_coefficient(n_curr, k, p_curr, r_curr, theta_I_curr, alpha_i_curr, alpha_SB_curr);
				double norm_coeff_curr = clamp(UTD_coeff_curr / SB_coeff_curr, 0, 1);
				double final_coeff_curr = (1.0f - lerp_curr) * UTD_coeff_curr + (lerp_curr)*norm_coeff_curr;

				/*if (i == 0)
				{
					std::cout << "2nd order: (" << listener_position.x() << ", " << listener_position.y() << ", " << listener_position.z() << ")\n\t";
					std::cout << UTD_coeff_curr << "; " << SB_coeff_curr << "; " << norm_coeff_curr << "; " << final_coeff_curr << " || ";
					std::cout << n_curr << "; " << p_curr << "; " << r_curr << "; " << alpha_i_curr << "; " << alpha_SB_curr << "; " << alpha_d_curr << "\n";
				}*/

				total_attenuations[i] = clamp(final_coeff_curr * final_coeff_curr, 0., 1.);
			}

			if (CGAL::cross_product(diffr_point_prevcurr.position - diffr_point_curr.position, diffr_point_curr.position - listener_position) * edge_axis_curr > 0)
			{
				shadow_listener_position = rotatePointAroundAxis(_on_curr_edge, edge_axis_curr, listener_position, alpha_SB_curr - alpha_d_curr);
			}
			else
			{
				shadow_listener_position = rotatePointAroundAxis(_on_curr_edge, -edge_axis_curr, listener_position, alpha_SB_curr - alpha_d_curr);
			}
		}
	}

	// остальные коэфф-ты затухания кроме последнего
	for (int i = 2; i < diffraction_points.size()-1; i++)
	{
		SoundDiffractionPathPoint& diffr_point_nextcurr = diffraction_points[i - 1];
		SoundDiffractionPathPoint& diffr_point_curr = diffraction_points[i];
		SoundDiffractionPathPoint& diffr_point_prevcurr = diffraction_points[i + 1];

		Vector s_f_n_curr = getVectorDirection(diffr_point_curr.source_plane.orthogonal_vector());
		Vector l_f_n_curr = getVectorDirection(diffr_point_curr.listener_plane.orthogonal_vector());

		Point _on_curr_edge = diffr_point_curr.diffr_edge->getStart();
		/*double _signed_dist = getSignedDistanceTo(Plane(_on_curr_edge, l_f_n_curr), _on_curr_edge,
			shadow_listener_position);
		if (_signed_dist < -epsilon_h)
		{
			s_f_n_curr = diffr_point_curr.diffr_edge->normal2;
			l_f_n_curr = diffr_point_curr.diffr_edge->normal1;
		}*/
		Vector edge_axis_curr = diffr_point_curr.diffr_edge->getAxis();
		if (CGAL::cross_product(s_f_n_curr, l_f_n_curr) * edge_axis_curr < 0)
		{
			edge_axis_curr = -diffr_point_curr.diffr_edge->getAxis();
		}


		//std::cout << "*********** edge_axis source_face_normal " << edge_axis << "; " << source_face_normal << ";\n";
		const Vector source_face_vector_curr = CGAL::cross_product(edge_axis_curr, s_f_n_curr);

		const double n_curr = (2. * Pi - angleBetween(-s_f_n_curr, l_f_n_curr)) / Pi;
		Vector source_direction_curr = diffr_point_prevcurr.position - diffr_point_curr.position;  //from apex point to source
		Vector listener_direction_curr = shadow_listener_position - diffr_point_curr.position; // from apex to receiver
		const double p_curr = getVectorLength(source_direction_curr);
		const double r_curr = getVectorLength(listener_direction_curr);

		if (p_curr < MATH_EPSILON || r_curr < MATH_EPSILON)
		{
			for (int i = 0; i < frequencies.size(); i++)
			{
				total_attenuations[i] = 1;
			}
		}
		else
		{
			source_direction_curr /= p_curr;
			listener_direction_curr /= r_curr;

			double theta_I_curr = angleBetween(source_direction_curr, edge_axis_curr);

			if (theta_I_curr > Pi * 0.5)
				theta_I_curr = Pi - theta_I_curr;

			const Vector s_dir_curr = projectToPlane(source_direction_curr, edge_axis_curr);
			const Vector r_dir_curr = projectToPlane(listener_direction_curr, edge_axis_curr);

			//std::cout << "*********** source_direction listener_direction " << source_direction << "; " << listener_direction << "; " << ";\n";
			//std::cout << "*********** s_dir r_dir " << s_dir << "; " << r_dir << "; " << source_face_vector << "; " << n << ";\n";

			const double alpha_i_curr = angleBetween(-s_dir_curr, source_face_vector_curr);
			const double alpha_d_curr = angleBetween(r_dir_curr, source_face_vector_curr) + Pi;
			const double alpha_SB_curr = alpha_i_curr + Pi + 0.001;

			const double lerp_curr = (n_curr * Pi - alpha_d_curr) / (n_curr * Pi - alpha_SB_curr);

			for (int i = 0; i < frequencies.size(); i++)
			{
				double lambda = speed_of_sound / frequencies[i];
				double k = 2. * Pi / lambda;

				double UTD_coeff_curr = UTD_coefficient(n_curr, k, p_curr, r_curr, theta_I_curr, alpha_i_curr, alpha_d_curr);
				double SB_coeff_curr = UTD_coefficient(n_curr, k, p_curr, r_curr, theta_I_curr, alpha_i_curr, alpha_SB_curr);
				double norm_coeff_curr = UTD_coeff_curr / SB_coeff_curr;
				double final_coeff_curr = (1.0f - lerp_curr) * UTD_coeff_curr + (lerp_curr)*norm_coeff_curr;

				//if (i == 0)
					//std::cout << final_coeff_curr << '\n';

				total_attenuations[i] *= clamp(final_coeff_curr * final_coeff_curr, 0., 1.);
			}

			if (CGAL::cross_product(diffr_point_prevcurr.position - diffr_point_curr.position, diffr_point_curr.position - shadow_listener_position) * edge_axis_curr > 0)
			{
				shadow_listener_position = rotatePointAroundAxis(_on_curr_edge, edge_axis_curr, shadow_listener_position, alpha_SB_curr - alpha_d_curr);
			}
			else
			{
				shadow_listener_position = rotatePointAroundAxis(_on_curr_edge, -edge_axis_curr, shadow_listener_position, alpha_SB_curr - alpha_d_curr);
			}
		}
	}
	// ....................
	/*			..
			 .. SB_1 (не угол, а плоскость)
		 ..	D1. . . . .	D2 . . . . L_imaginary
	   ..  /  \		   / | . SB_2 (не угол, а плоскость)
	 S	  /    \	  /  |	 .
		 /      \	 /   |	   L

	*/

	// последний коэф затухания
	{
		SoundDiffractionPathPoint& diffr_point_nextcurr = diffraction_points[diffraction_points.size() - 2];
		SoundDiffractionPathPoint& diffr_point_curr = diffraction_points[diffraction_points.size() - 1];

		Vector s_f_n_curr = getVectorDirection(diffr_point_curr.source_plane.orthogonal_vector());
		Vector l_f_n_curr = getVectorDirection(diffr_point_curr.listener_plane.orthogonal_vector());

		Point _on_curr_edge = diffr_point_curr.diffr_edge->getStart();
		/*double _signed_dist = getSignedDistanceTo(Plane(_on_curr_edge, l_f_n_curr), _on_curr_edge,
			shadow_listener_position);
		if (_signed_dist < -epsilon_h)
		{
			s_f_n_curr = diffr_point_curr.diffr_edge->normal2;
			l_f_n_curr = diffr_point_curr.diffr_edge->normal1;
		}*/
		Vector edge_axis_curr = diffr_point_curr.diffr_edge->getAxis();
		if (CGAL::cross_product(s_f_n_curr, l_f_n_curr) * edge_axis_curr < 0)
		{
			edge_axis_curr = -diffr_point_curr.diffr_edge->getAxis();
		}

		//std::cout << "*********** edge_axis source_face_normal " << edge_axis << "; " << source_face_normal << ";\n";
		const Vector source_face_vector_curr = CGAL::cross_product(edge_axis_curr, s_f_n_curr);

		const double n_curr = (2. * Pi - angleBetween(-s_f_n_curr, l_f_n_curr)) / Pi;
		Vector source_direction_curr = source_position - diffr_point_curr.position;  //from apex point to source
		Vector listener_direction_curr = shadow_listener_position - diffr_point_curr.position; // from apex to receiver
		const double p_curr = getVectorLength(source_direction_curr);
		const double r_curr = getVectorLength(listener_direction_curr);

		if (p_curr < MATH_EPSILON || r_curr < MATH_EPSILON)
		{
			for (int i = 0; i < frequencies.size(); i++)
			{
				total_attenuations[i] = 1;
			}
		}
		else
		{
			source_direction_curr /= p_curr;
			listener_direction_curr /= r_curr;

			double theta_I_curr = angleBetween(source_direction_curr, edge_axis_curr);

			if (theta_I_curr > Pi * 0.5)
				theta_I_curr = Pi - theta_I_curr;

			const Vector s_dir_curr = projectToPlane(source_direction_curr, edge_axis_curr);
			const Vector r_dir_curr = projectToPlane(listener_direction_curr, edge_axis_curr);

			//std::cout << "*********** source_direction listener_direction " << source_direction << "; " << listener_direction << "; " << ";\n";
			//std::cout << "*********** s_dir r_dir " << s_dir << "; " << r_dir << "; " << source_face_vector << "; " << n << ";\n";

			const double alpha_i_curr = angleBetween(-s_dir_curr, source_face_vector_curr);
			const double alpha_d_curr = angleBetween(r_dir_curr, source_face_vector_curr) + Pi;
			const double alpha_SB_curr = alpha_i_curr + Pi + 0.001;

			const double lerp_curr = clamp((n_curr * Pi - alpha_d_curr) / (n_curr * Pi - alpha_SB_curr), 0, 1);

			for (int i = 0; i < frequencies.size(); i++)
			{
				double lambda = speed_of_sound / frequencies[i];
				double k = 2. * Pi / lambda;

				double UTD_coeff_curr = UTD_coefficient(n_curr, k, p_curr, r_curr, theta_I_curr, alpha_i_curr, alpha_d_curr);
				double SB_coeff_curr = UTD_coefficient(n_curr, k, p_curr, r_curr, theta_I_curr, alpha_i_curr, alpha_SB_curr);
				double norm_coeff_curr = clamp(UTD_coeff_curr / SB_coeff_curr, 0, 1);
				double final_coeff_curr = (1.0f - lerp_curr) * UTD_coeff_curr + (lerp_curr)*norm_coeff_curr;

				/*if (i == 0)
				{
					std::cout << "\t" << UTD_coeff_curr << "; " << SB_coeff_curr << "; " << norm_coeff_curr << "; " << final_coeff_curr << " || ";
					std::cout << n_curr << "; " << p_curr << "; " << r_curr << "; " << alpha_i_curr << "; " << alpha_SB_curr << "; " << alpha_d_curr << "\n";
					if (listener_position.x() == 2.2 && listener_position.y() == 2.7 && listener_position.z() == 0)
					{
						std::cout << "\n";
					}
				}*/

				total_attenuations[i] *= clamp(final_coeff_curr * final_coeff_curr, 0., 1.);
			}
		}
	}

	FrequencyBandResponse result(frequencies.size());
	for (int i = 0; i < frequencies.size(); i++)
	{
		result.changeAttenuation(i, total_attenuations[i]);
	}
	return result;
}


double UTD_coefficient(double n, double k, double p, double r, double thetaI, double alphaI, double alphaD)
{
	complex_double c1 = UTD_freqTerm(n, k, thetaI);
	complex_double F1 = UTD_estimateF(k * UTD_L(p, r, thetaI) * UTD_alpha(alphaD - alphaI, n, 1));
	complex_double F2 = UTD_estimateF(k * UTD_L(p, r, thetaI) * UTD_alpha(alphaD - alphaI, n, -1));
	complex_double F3 = UTD_estimateF(k * UTD_L(p, r, thetaI) * UTD_alpha(alphaD + alphaI, n, 1));
	complex_double F4 = UTD_estimateF(k * UTD_L(p, r, thetaI) * UTD_alpha(alphaD + alphaI, n, -1));

	double cot1 = UTD_cotan(Pi + (alphaD - alphaI), double(2) * n);
	double cot2 = UTD_cotan(Pi - (alphaD - alphaI), double(2) * n);
	double cot3 = UTD_cotan(Pi + (alphaD + alphaI), double(2) * n);
	double cot4 = UTD_cotan(Pi - (alphaD + alphaI), double(2) * n);

	complex_double coeff = F1 * cot1 + F2 * cot2 + F3 * cot3 + F4 * cot4;
	coeff = coeff * c1;
	coeff = UTD_euler(-k * r) * coeff;

	coeff *= std::sqrt(UTD_sphereDisKouyoumjian(r, p));

	return abs(coeff);
}

double UTD_alpha(double beta, double n, int nSign)
{
	int N = UTD_N(beta, n, nSign);
	double numer = 2. * Pi * n * N - beta;
	double denom = 2.;
	double cosine = std::cos(numer / denom);
	double alpha = 2. * cosine * cosine;

	return alpha;
}

double UTD_L(double p, double r, double thetaI)
{
	double sine = std::sin(thetaI);
	double L = UTD_sphereDis(r, p) * sine * sine;

	return L;
}

int UTD_N(double beta, double n, int nSign)
{
	if (nSign > 0)
	{
		if (beta <= Pi * (n - 1))
			return 0;
		else  // beta > Pi*(n-1)
			return 1;
	}
	else
	{
		if (beta < Pi * (1 - n))
			return -1;
		else if (beta >= Pi * (1 - n) && beta <= Pi * (1 + n)) // Pi*(1-n) <= beta <= Pi*(1+n)
			return 0;
		else // beta > Pi*(1+n)
			return 1;
	}
}

double UTD_cotan(double numer, double denom)
{
	return cotangent(numer / denom);
}

complex_double UTD_euler(double x)
{
	return complex_double(std::cos(x), std::sin(x));
}

complex_double UTD_estimateF(double X)
{
	double numer_phase = X;
	double denom_phase = X + .4;
	complex_double phase_term = UTD_euler(Pi * 0.25 * std::sqrt(numer_phase / denom_phase));
	complex_double F;

	if (X < 0.8)
	{
		double numer = std::sqrt(X);
		double denom = 0.7 * std::sqrt(X) + 1.2;
		F = phase_term * std::sqrt(Pi * X) * (1 - numer / denom);
	}
	else
	{
		double numer = 0.8;
		double denom = (X + 1.25) * (X + 1.25);
		F = (1. - numer / denom) * phase_term;
	}

	return F;
}

complex_double UTD_freqTerm(double n, double k, double thetaI)
{
	complex_double numer = UTD_euler(-Pi * 0.25);
	double denom = 2. * n * std::sqrt(2. * Pi * k) * std::sin(thetaI);

	return -numer / denom;
}

double UTD_sphereDisKouyoumjian(double r, double p)
{
	double numer = p;
	double denom = r * (p + r);

	return numer / denom;
}

double UTD_sphereDis(double r, double p)
{
	double numer = p * r;
	double denom = p + r;

	return numer / denom;
}

double cotangent(double x)
{
	if (std::abs(x) < MATH_EPSILON)
		return MATH_DOUBLE_MAX;

	return 1. / std::tan(x);
}
