#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

// Для рисования
#include <CGAL/draw_surface_mesh.h>
#include <fstream>
#include <locale>
#include "SoundScene.h"

#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>

#define DRAW_DIFFR_EDGES false
#define OUTPUT_LOG false

int custom_geom(int argc, char* argv[])
{
	std::string scene_name;
	std::string answer = "y";
	double sx, sy, sz;
	double dx = 0.02, dy = 0.02, dz = 0.02;
	double x_low, x_high, y_low, y_high;
	unsigned list_plane = 2;
	double plane_val = 0;
	unsigned phong_flag = false;
	unsigned num_frequencies = 3;
	bool is_2d = false;
	size_t max_diffr_order;
	SoundScene scene;
	scene.getInputParams(scene_name, answer,
		sx, sy, sz,
		dx, dy,
		x_low, x_high,
		y_low, y_high,
		list_plane, plane_val, phong_flag, is_2d, max_diffr_order);
	scene.max_diffraction_order = max_diffr_order;
	scene.is_2d = is_2d;
	//scene.epsilon_h = std::max(2*std::max(dx, dy), scene.epsilon_h);

	std::string results_codename = scene.createOutputFilename(scene_name);
	scene.saveOutputFilename(results_codename);

	std::string full_filename = "Z:/course_7/blender/" + scene_name + ".obj";
	std::cout << full_filename << "\n";

	scene.setFrequencies(num_frequencies);
	scene.saveFrequencies("saved_IR/" + results_codename);
	scene.readMeshFromFile(full_filename);
	scene.mapMaterials();

	bool load_bde = false;
	std::cout << "Load latest save? (y/n) " << answer << "\n";

	if (answer != "y")
	{
		// помечаем диффр рёбра
		scene.markDiffractionEdges();

		// объединяем вместе в одно большие
		scene.combineDiffractionEdges();
		// сохраняем
		scene.saveDiffractionEdges(scene_name);

		std::cout << scene.diffr_edge_graph.edges.size() << "\n";

		// Draw!
		// CGAL::draw(scene.mesh);

		// построение графа видимости рёбер
		std::cout << "buildDiffractionEdgeGraph\n";
		scene.buildDiffractionEdgeGraph();

		// сохраняем
		scene.saveDiffractionEdgeGraph(scene_name);
	}
	else
	{
		// подгружаем
		scene.loadDiffractionEdges(scene_name);
		scene.loadDiffractionEdgeGraph(scene_name);
	}
	Point b1, b2;
	scene.getSceneBoundingBox(b1, b2);

	double min_diffraction_edge_length = 0.5;

	// просто я взял из головы
	double diffaction_flag_dir_variation_threshold = 0.95;

	// для проверки
	std::cout << scene.diffr_edge_graph.edges.size() << "; "
		<< scene.diffr_edge_graph.diffr_edge_neighbors.size() << "\n";

	double source_pos_x;
	double source_pos_y;
	double listener_pos_x;
	double listener_pos_y;

	std::cout << "start_addDiffractionFlags\n";
	scene.addDiffractionFlags();
	std::cout << "end_addDiffractionFlags\n";

	// source
	std::cout << "source position = " << sx << ", " << sy << ", " << sz << "\n";
	//(3.5, 4., 2.)
	scene.setSoundSource(PointSoundSource(Point(sx, sy, sz)));
	auto vd_source = scene.mesh.add_vertex(scene.sound_source.position);

	//scene.epsilon_h = 2*std::max(dx, dy);

	std::cout << "dx = " << dx;
	std::cout << ", dy = " << dy << "\n";
	std::cout << "x = " << x_low << ".." << x_high << "\n";
	std::cout << "y = " << y_low << ".." << y_high << "\n";
	size_t x_dots_number = std::floor((x_high - x_low) / dx);
	size_t y_dots_number = std::floor((y_high - y_low) / dy);

	if (x_low + x_dots_number * dx <= x_high)
		x_dots_number++;

	if (y_low + y_dots_number * dy <= y_high)
		y_dots_number++;

	std::cout << "x_dots_number = " << x_dots_number << "\n";
	std::cout << "y_dots_number = " << y_dots_number << "\n";

	std::cout << "max_diffr_order = " << scene.max_diffraction_order << "\n";

	size_t x_curr_idx = 0;
	size_t y_curr_idx = 0;

	double x_curr = x_low;
	double y_curr = y_low;

	double _start = x_curr; double _finish = x_high; int _full_len = 30; int _curr_len = -1; int _curr_pcnt = -1;

	//CGAL::draw(scene.mesh);

	for (x_curr_idx = 0; x_curr_idx < x_dots_number; x_curr_idx++)
	{
		x_curr = x_low + x_curr_idx * dx;
		//std::cout << "x_curr = " << x_curr << "\n";
		if ((int)std::floor((x_curr - _start) / (_finish - _start) * 100) > _curr_pcnt)
		{
			_curr_len = (int)std::floor((x_curr - _start) / (_finish - _start) * _full_len);
			_curr_pcnt = (int)std::floor((x_curr - _start) / (_finish - _start) * 100);
			std::cout << (int)std::floor((x_curr - _start) / (_finish - _start) * 100) << "% : ";
			for (int i = 0; i < _curr_len; i++)
			{
				std::cout << "#";
			}
			for (int i = _curr_len; i < _full_len; i++)
			{
				std::cout << ".";
			}
			std::cout << "\n";
		}
		y_curr = y_low;
		for (y_curr_idx = 0; y_curr_idx < y_dots_number; y_curr_idx++)
		{
			y_curr = y_low + y_curr_idx * dy;
			
			Point curr_list_pos = Point(x_curr, y_curr, plane_val);

			if (list_plane == 0) curr_list_pos = Point(plane_val, x_curr, y_curr);
			if (list_plane == 1) curr_list_pos = Point(y_curr, plane_val, x_curr);
			if (list_plane == 2) curr_list_pos = Point(x_curr, y_curr, plane_val);

			if (scene.isBackSideVisible(curr_list_pos))
			{
				continue;
			}
			if (OUTPUT_LOG)
			{
				std::cout << "(" << x_curr << ", " << y_curr << ", " << plane_val << ")\n";
			}
			// listener

			double move_listener_epsilon = std::min(dx, dy)*0.1;
			Point moved_curr_list_pos = Point(
				curr_list_pos.x() + randomUniformNextVal() * move_listener_epsilon,
				curr_list_pos.y() + randomUniformNextVal() * move_listener_epsilon,
				curr_list_pos.z() + randomUniformNextVal() * move_listener_epsilon);

			scene.setSoundListener(PointSoundListener(moved_curr_list_pos, x_curr_idx, y_curr_idx));

			auto vd_listener = scene.mesh.add_vertex(scene.sound_listener.position);

			// Draw!
			//CGAL::draw(scene.mesh);

			double frequency = scene.frequencies[0];

			bool has_direct_path = scene.addDirectPath();
			if (phong_flag < 2)
			{
				// propagate
				if (OUTPUT_LOG)
				{
					std::cout << "propagateSound\n";
				}
				for (size_t i = 0; i < scene.max_diffr_ray_number; i++)
				{
					scene.propagateSound();
				}
				if (OUTPUT_LOG)
				{
					std::cout << "END_propagateSound\n";
				}
			}
			if (phong_flag > 0)
			{
				if (OUTPUT_LOG)
				{
					std::cout << "START_addSpecularPath\n";
				}
				for (size_t i = 0; i < scene.max_ray_number; i++)
				{
					scene.addSpecularPath(frequency);
				}
				if (OUTPUT_LOG)
				{
					std::cout << "END_addSpecularPath\n";
				}
			}

			// Draw!
			if (CGAL_DRAW_MESH)
				CGAL::draw(scene.mesh);

			scene.sound_listener.archiveDiffractionIRs();
			scene.sound_listener.archivePhongIRs();

			if (OUTPUT_LOG)
			{
				double diffr_energy = 0;
				std::cout << "Num of DIFFR responses: " << scene.sound_listener.received_diffr_IRs_archived[0].size() << "\n";
				for (int i = 0; i < scene.sound_listener.received_diffr_IRs_archived[0].size(); i++)
				{
					std::cout << "Time: " << scene.sound_listener.received_diffr_IRs_archived[0][i].first
						<< " \tEnergy: " << scene.sound_listener.received_diffr_IRs_archived[0][i].second << ";\n";
					diffr_energy += scene.sound_listener.received_diffr_IRs_archived[0][i].second;
				}
				std::cout << "------------------------------------\n";


				double phong_energy = 0;
				std::cout << "Num of PHONG responses: " << scene.sound_listener.received_phong_IRs_archived.size() << "\n";
				for (int i = 0; i < scene.sound_listener.received_phong_IRs_archived.size(); i++)
				{
					std::cout << "Time: " << scene.sound_listener.received_phong_IRs_archived[i].first
						<< " \tEnergy: " << scene.sound_listener.received_phong_IRs_archived[i].second << ";\n";
					phong_energy += scene.sound_listener.received_phong_IRs_archived[i].second;
				}
				std::cout << "diffr_energy = " << diffr_energy << "\n";
				std::cout << "phong_energy = " << phong_energy << "\n";
				std::cout << "direct_energy = " << (scene.sound_listener.received_direct_IRs_archived.size() ? scene.sound_listener.received_direct_IRs_archived[0].second : 0.) << "\n";
			}

			scene.sound_listener.appendListenerToFile("saved_IR/" + results_codename);
			scene.sound_listener.allClear();

			scene.mesh.remove_vertex(vd_listener);
		}
	}

	return 0;
}



int main(int argc, char* argv[])
{
	auto start_time_measure = std::chrono::steady_clock::now();
	custom_geom(argc, argv);
	auto finish_time_measure = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::minutes>(finish_time_measure - start_time_measure);
	std::cout << "The time: " << elapsed_ms.count() << " min\n";
}