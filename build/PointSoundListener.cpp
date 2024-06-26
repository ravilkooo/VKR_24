#include "PointSoundListener.h"

Vector PointSoundListener::sample_ray(bool is_2d)
{
	if (is_2d)
		return uniformSampleCircle(randomUniformNextVal());
	else
		return uniformSampleSphere(randomUniformNextVal(), randomUniformNextVal());
	//return Vector();
}

PointSoundListener::PointSoundListener()
	: position(Point(0, 0, 0))
{
	IR_cnt = std::vector<unsigned>(3);
	for (int i = 0; i < 3; i++)
	{
		IR_cnt[i] = 0;
	}

}

PointSoundListener::PointSoundListener(Point pos)
	: position(pos)
{
}

PointSoundListener::PointSoundListener(Point pos, size_t x_idx, size_t y_idx)
	: position(pos), x_idx(x_idx), y_idx(y_idx)
{
}

void PointSoundListener::archiveDiffractionIRs()
{
	if (received_diffraction_IRs.empty())
	{
		received_diffr_IRs_archived.clear();
		return;
	}
	size_t size = received_diffraction_IRs.size();
	size_t num_frequencies = received_diffraction_IRs[0].energy.attenuations.size();
	//std::cout << "num_frequencies = " << num_frequencies << "; size = " << size << "\n";
	frequencies = received_diffraction_IRs[0].energy.getFrequencies();


	received_diffr_IRs_archived.resize(num_frequencies);
	for (size_t i = 0; i < num_frequencies; i++)
	{
		//std::cout << "freq[" << i << "] = " << frequencies[i] << ";\n";
		received_diffr_IRs_archived[i].resize(size);
		for (size_t j = 0; j < size; j++)
		{
			//std::cout << "time = " << received_diffraction_IRs[j].time << "; attenuation = " << received_diffraction_IRs[j].energy.attenuations[i] << "\n";
			received_diffr_IRs_archived[i][j] = std::make_pair(
				received_diffraction_IRs[j].time,
				received_diffraction_IRs[j].energy.attenuations[i]);
		}
	}
	received_diffraction_IRs.clear();
	
	for (size_t i = 0; i < num_frequencies; i++)
	{
		std::sort(received_diffr_IRs_archived[i].begin(), received_diffr_IRs_archived[i].end(),
			[](const std::pair<double, double>& a, const std::pair<double, double>& b) {
			return a.first < b.first;
		});
	}
	return;
}
void PointSoundListener::archivePhongIRs()
{
	std::sort(received_phong_IRs_archived.begin(), received_phong_IRs_archived.end(),
		[](const std::pair<double, double>& a, const std::pair<double, double>& b) {
		return a.first < b.first;
	});
}
void PointSoundListener::saveToStream(std::ofstream& out, std::vector<std::pair<double, double>>& received_IRs_archived) const
{
	// Сохранение позиции
	out.write(reinterpret_cast<const char*>(&(position.x())), sizeof((position.x())));
	out.write(reinterpret_cast<const char*>(&(position.y())), sizeof((position.y())));
	out.write(reinterpret_cast<const char*>(&(position.z())), sizeof((position.z())));

	// Сохранение индекса
	out.write(reinterpret_cast<const char*>(&x_idx), sizeof(x_idx));
	out.write(reinterpret_cast<const char*>(&y_idx), sizeof(y_idx));

	// Сохранение размера вектора
	size_t size = received_IRs_archived.size();
	out.write(reinterpret_cast<const char*>(&size), sizeof(size));

	// Сохранение пар значений
	out.write(reinterpret_cast<const char*>(received_IRs_archived.data()), sizeof(std::pair<double, double>) * size);
}
void PointSoundListener::appendListenerToFile(const std::string& filename)
{
	std::ofstream out(filename + "_phong_IRs", std::ios::binary | std::ios::app);
	saveToStream(out, received_phong_IRs_archived);
	out.close();
	
	if (received_diffr_IRs_archived.empty())
	{
		received_diffr_IRs_archived.resize(frequencies.size());
		for (size_t i = 0; i < frequencies.size(); i++)
		{
			received_diffr_IRs_archived[i].resize(0);
		}
	}

	for (size_t i = 0; i < frequencies.size(); i++)
	{
		out = std::ofstream(filename + "_diffr_IRs_" + std::to_string(int(std::round(frequencies[i]))), std::ios::binary | std::ios::app);
		saveToStream(out, received_diffr_IRs_archived[i]);
		out.close();
	}

	out = std::ofstream(filename + "_direct_IRs", std::ios::binary | std::ios::app);
	saveToStream(out, received_direct_IRs_archived);
	out.close();
}

void PointSoundListener::allClear()
{
	received_phong_IRs_archived.clear();
	received_diffr_IRs_archived.clear();
	received_direct_IRs_archived.clear();
}