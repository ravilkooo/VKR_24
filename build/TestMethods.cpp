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

//test_flags(argc, argv);
//edge_anim(argc, argv);
//test_saveLoad("test_saveLoad", true);
//test_saveLoad("input", true, true);

int test_computePointOfClosestApproach()
{
	std::string filename = "Z:/course_7/blender/1clin.obj";

	SoundScene scene;
	scene.readMeshFromFile(filename);

	double min_diffraction_edge_length = 0.5;

	// просто я взял из головы
	double diffaction_flag_dir_variation_threshold = 0.95;

	// помечаем диффр рёбра
	scene.markDiffractionEdges();

	// объединяем вместе в одно большие
	scene.combineDiffractionEdges();

	// построение графа видимости рёбер
	std::cout << "buildDiffractionEdgeGraph\n";
	scene.buildDiffractionEdgeGraph();

	std::cout << "addDiffractionFlags\n";
	scene.addDiffractionFlags();

	BigDiffractionEdge& diffr_edge = scene.diffr_edge_graph.edges[0];

	for (int i = -4; i < 5; i++)
	{
		PointSoundSource sound_source = PointSoundSource(Point(i, 1, 1.5));
		scene.setSoundSource(sound_source);
		auto v1 = scene.mesh.add_vertex(scene.sound_source.position);

		for (int j = -4; j < 5; j++)
		{
			// Draw!
			//CGAL::draw(scene.mesh);

			// listener
			PointSoundListener sound_listener = PointSoundListener(Point(j, 1, -1.5));
			scene.setSoundListener(sound_listener);
			auto v3 = scene.mesh.add_vertex(scene.sound_listener.position);

			double t = scene.computePointOfClosestApproach(
				diffr_edge.getStart(),
				diffr_edge.getAxis(),
				sound_listener.position,
				getVectorDirection(sound_source.position - sound_listener.position));

			auto v2 = scene.mesh.add_vertex(diffr_edge.getStart()
				+ t * diffr_edge.getAxis());
			halfedge_descriptor e12 = scene.mesh.add_edge(v1, v2);
			halfedge_descriptor e23 = scene.mesh.add_edge(v2, v3);

			CGAL::draw(scene.mesh);

			scene.mesh.remove_edge(scene.mesh.edge(e12));
			scene.mesh.remove_edge(scene.mesh.edge(e23));

			scene.mesh.remove_vertex(v2);
			scene.mesh.remove_vertex(v3);


		}
		scene.mesh.remove_vertex(v1);

	}
	return 0;
}

int test_checkSourceInShadowRegion()
{
	std::string filename = "Z:/course_7/blender/1clin.obj";

	SoundScene scene;
	scene.readMeshFromFile(filename);
	//scene.mapMaterials();

	double min_diffraction_edge_length = 0.5;

	// просто я взял из головы
	double diffaction_flag_dir_variation_threshold = 0.95;

	// помечаем диффр рёбра
	scene.markDiffractionEdges();

	// объединяем вместе в одно большие
	scene.combineDiffractionEdges();

	// построение графа видимости рёбер
	std::cout << "buildDiffractionEdgeGraph\n";
	scene.buildDiffractionEdgeGraph();

	std::cout << "addDiffractionFlags\n";
	scene.addDiffractionFlags();

	BigDiffractionEdge& diffr_edge = scene.diffr_edge_graph.edges[0];
	auto v2 = scene.mesh.add_vertex(diffr_edge.getStart()
		+ diffr_edge.getLength() * 0.5 * diffr_edge.getAxis());

	for (int i = 0; i < 6; i++)
	{
		PointSoundSource sound_source = PointSoundSource(Point(0., i, 1.5));
		scene.setSoundSource(sound_source);
		auto v1 = scene.mesh.add_vertex(scene.sound_source.position);
		halfedge_descriptor e12 = scene.mesh.add_edge(v1, v2);

		std::vector<vertex_descriptor> v_arr;
		std::vector<halfedge_descriptor> e_arr;

		for (int j = 0; j < 6; j++)
		{
			// Draw!
			//CGAL::draw(scene.mesh);

			// listener
			PointSoundListener sound_listener = PointSoundListener(Point(0, j, -1.5));
			scene.setSoundListener(sound_listener);
			v_arr.push_back(scene.mesh.add_vertex(scene.sound_listener.position));

			ImpulseResponse curr_IR;

			curr_IR.diffraction_points_list.push_back(SoundDiffractionPathPoint(sound_listener.position));
			curr_IR.diffraction_points_list.push_back(SoundDiffractionPathPoint(diffr_edge.getStart()
				+ diffr_edge.getLength() * 0.5 * diffr_edge.getAxis(), 1));

			{
				int log_flag = 0;
				if (log_flag > 0) std::cout << "recursiveDiffraction_" << log_flag++ << "\n";

				const Point& source_position = sound_source.position;
				SoundDiffractionPathPoint& last_path_point = curr_IR.diffraction_points_list[curr_IR.diffraction_points_list.size() - 2];
				SoundDiffractionPathPoint& this_path_point = curr_IR.diffraction_points_list[curr_IR.diffraction_points_list.size() - 1];
				const Point& last_listener_image_position = last_path_point.position;
				const Point& listener_image_position = this_path_point.position;

				if (log_flag > 0) std::cout << "recursiveDiffraction_" << log_flag++ << "\n";

				// Determine which side of the edge the listener image position is on.
				// true == faces plane 1, false == faces plane 2.
				double plane1_distance = getSignedDistanceTo(diffr_edge.normal1, diffr_edge.getStart(), last_listener_image_position);
				double plane2_distance = getSignedDistanceTo(diffr_edge.normal2, diffr_edge.getStart(), last_listener_image_position);

				bool listener_orientation = plane1_distance > plane2_distance && plane1_distance > 0.;
				std::cout << "  listener_orientation = " << listener_orientation << "\n";

				if (log_flag > 0) std::cout << "recursiveDiffraction_" << log_flag++ << "\n";

				// Get the free triangle vertex on the same side of the edge as the image source position.
				Point triangle_free_vertex = listener_orientation ?
					diffr_edge.getFreeVertex1() : diffr_edge.getFreeVertex2();

				if (log_flag > 0) std::cout << "recursiveDiffraction_" << log_flag++ << "\n";

				// Compute the shadow region boundary plane for the edge and last source image position.
				Plane shadow_boundary(last_listener_image_position, diffr_edge.getStart(), diffr_edge.getStart() + diffr_edge.getAxis());

				Vector shadow_boundary_normal = getVectorDirection(shadow_boundary.orthogonal_vector());

				if (log_flag > 0) std::cout << "recursiveDiffraction_" << log_flag++ << "\n";

				// Make sure that the shadow boundary plane's normal points in the right direction.
				if (getSignedDistanceTo(shadow_boundary_normal, last_listener_image_position, triangle_free_vertex) < 0.)
					shadow_boundary = Plane(last_listener_image_position, diffr_edge.getStart() + diffr_edge.getAxis(), diffr_edge.getStart());

				if (log_flag > 0) std::cout << "recursiveDiffraction_" << log_flag++ << "\n";

				// Get the plane that faces the image source position.
				const Plane& listener_plane = listener_orientation ? diffr_edge.getPlane1() : diffr_edge.getPlane2();
				this_path_point.listener_plane = listener_plane;

				if (log_flag > 0) std::cout << "recursiveDiffraction_" << log_flag++ << "\n";

				// Get the plane which defines the triangle boundary on the opposite side of the edge.
				const Plane& opposite_plane = listener_orientation ? diffr_edge.getPlane2() : diffr_edge.getPlane1();
				this_path_point.source_plane = opposite_plane;

				if (log_flag > 0) std::cout << "recursiveDiffraction_" << log_flag++ << "\n";

				//SoundDiffractionPathPoint pathIDPoint(SoundDiffractionPathPoint::EDGE_DIFFRACTION,
				//	ObjectSpaceTriangle(listener_orientation ? edge.edge->triangle1 : edge.edge->triangle2, query.object),
				//	listener_orientation ? edge.edge->edge_index1 : edge.edge->edge_index2);

				//// Add a new point to the path ID for this edge diffraction.
				//path_id.addPoint(path_id_point);

				//*********************************************************************
				// Validate the path to the source, if possible.

				// Determine whether or not the source is in the shadow region by testing it against the boundary planes.
				bool source_in_shadow_region = getSignedDistanceTo(shadow_boundary_normal, last_listener_image_position, source_position) > 0 &&
					getSignedDistanceTo(opposite_plane, diffr_edge.getStart(), source_position) > 0;
				std::cout << "  source_in_shadow_region = " << source_in_shadow_region << "\n";

				if (source_in_shadow_region)
				{
					e_arr.push_back(scene.mesh.add_edge(v2, v_arr.back()));
				}

			}

		}

		CGAL::draw(scene.mesh);

		for (int j = 0; j < e_arr.size(); j++)
		{
			scene.mesh.remove_edge(scene.mesh.edge(e_arr[j]));
			scene.mesh.remove_vertex(v_arr[j]);
		}

		scene.mesh.remove_edge(scene.mesh.edge(e12));
		scene.mesh.remove_vertex(v1);

	}
	return 0;
}

int test_PhongReflection()
{
	std::string filename = "Z:/course_7/blender/test_phong.obj";

	SoundScene scene;
	scene.readMeshFromFile(filename);
	scene.mapMaterials();

	std::cout << "reflectivity = " << scene.materials_list[0].reflectivity.attenuations[3] << "\n";
	std::cout << "scattering = " << scene.materials_list[0].scattering.attenuations[3] << "\n";

	for (int i = 0; i < scene.materials_list[0].reflectivity.attenuations.size(); i++)
	{
		std::cout << scene.materials_list[0].reflectivity.attenuations[i] << "; ";
	}
	std::cout << "\n";
	for (int i = 0; i < scene.materials_list[0].scattering.attenuations.size(); i++)
	{
		std::cout << scene.materials_list[0].scattering.attenuations[i] << "; ";
	}

	std::cout << "\n" << scene.materials_list[0].reflectivity.getClosestValue(300.) << "\n";
	std::cout << scene.materials_list[0].scattering.getClosestValue(300.) << "\n";

	double min_diffraction_edge_length = 0.5;

	// просто я взял из головы
	double diffaction_flag_dir_variation_threshold = 0.95;

	// помечаем диффр рёбра
	scene.markDiffractionEdges();

	// объединяем вместе в одно большие
	scene.combineDiffractionEdges();

	std::cout << scene.diffr_edge_graph.edges.size() << "!!!\n";

	// построение графа видимости рёбер
	std::cout << "buildDiffractionEdgeGraph\n";
	scene.buildDiffractionEdgeGraph();

	std::cout << "start_addDiffractionFlags\n";
	scene.addDiffractionFlags();
	std::cout << "end_addDiffractionFlags\n";

	// Draw!
	//CGAL::draw(scene.mesh);

	// source
	scene.setSoundSource(PointSoundSource(Point(2, 0, 1)));
	scene.mesh.add_vertex(scene.sound_source.position);

	// listener
	scene.setSoundListener(PointSoundListener(Point(-2, 0, 1)));
	scene.mesh.add_vertex(scene.sound_listener.position);

	// propagate

	std::cout << "START_addSpecularPath\n";

	double frequency = 300.;
	double attenuation = 1;
	Vector ray_dir = getVectorDirection(Vector(-2, 0, -1));
	Ray ray(scene.sound_source.position, ray_dir);
	FlagSkip skip(scene, Mesh::null_face(), true);
	int refl_ray_number = 0;

	while (refl_ray_number < 1000)
	{
		Ray_intersection intersection = scene.face_tree.first_intersection(ray, skip);
		bool has_isect = ((intersection) && (boost::get<Point>(&(intersection->first))));

		if (has_isect)
		{
			Point intrersection_pos = Point(*boost::get<Point>(&(intersection->first)));

			double intersection_distance = getVectorLength(intrersection_pos - ray.source());
			face_descriptor intersection_fd = intersection->second;

			auto face_normal =
				getVectorDirection(CGAL::Polygon_mesh_processing::compute_face_normal(intersection_fd, scene.mesh));

			// Compute the dot product of the triangle's normal with the incoming ray's direction.
			double ray_dot_normal = ray_dir * face_normal;

			// Flip the normal if it points in the same direction as the ray.
			if (ray_dot_normal > 0.)
			{
				face_normal = -face_normal;
				ray_dot_normal = -ray_dot_normal;
			}

			// Get the material properties for the intersected triangle.
			SoundMaterial* material = &(scene.materials_list[scene.faces_material_map[intersection_fd]]);

			double pdf_fwd, pdf_rev, brdf_fwd, brdf_rev;
			SoundPhongPathPoint::ReflectionType refl_type;

			SoundPhongPathPoint new_point(intrersection_pos);
			new_point.material = material;
			new_point.normal = face_normal;
			new_point.frequency = frequency;
			Vector out_dir = new_point.getPhongReflection(
				-ray_dir, refl_type,
				pdf_fwd, pdf_rev, brdf_fwd, brdf_rev);
			if (refl_type == SoundPhongPathPoint::ReflectionType::NOCONTRIBUTION)
			{
				continue;
			}

			if (refl_type == SoundPhongPathPoint::ReflectionType::DIFFUSE)
			{
				scene.addArrow(intrersection_pos, intrersection_pos + out_dir);
			}
			else if (refl_type == SoundPhongPathPoint::ReflectionType::SPECULAR)
			{
				scene.addArrow(intrersection_pos, intrersection_pos + out_dir * 2);
			}
			refl_ray_number++;
		}
	}
	std::cout << "END_addSpecularPath\n";

	// Draw!
	CGAL::draw(scene.mesh);

	return 0;
}

int test_saveLoad(std::string filename, bool save_data, bool test_graph = false)
{
	// SAVE
	if (save_data)
	{
		std::string full_filename = "Z:/course_7/blender/" + filename + ".obj";

		SoundScene scene;
		scene.readMeshFromFile(full_filename);
		scene.mapMaterials();

		// помечаем диффр рёбра
		scene.markDiffractionEdges();

		// объединяем вместе в одно большие
		scene.combineDiffractionEdges();
		// сохраняем
		std::cout << "saveDiffractionEdges\n";
		scene.saveDiffractionEdges(filename);

		std::cout << scene.diffr_edge_graph.edges.size() << "\n";

		std::cout << scene.diffr_edge_graph.edges.size() << "; "
			<< scene.diffr_edge_graph.diffr_edge_neighbors.size() << "\n";

		CGAL::draw(scene.mesh);

		if (test_graph)
		{
			// построение графа видимости рёбер
			std::cout << "buildDiffractionEdgeGraph\n";
			scene.buildDiffractionEdgeGraph();
			// сохраняем
			std::cout << "saveDiffractionEdgeGraph\n";
			scene.saveDiffractionEdgeGraph(filename);

			// для проверки
			for (int i = 0; i < scene.total_num_edges; i++)
			{
				std::cout << i << ":\t";
				for (int j = 0; j < scene.diffr_edge_graph.edges[i].num_neighbors; j++)
				{
					int neigh_k = scene.diffr_edge_graph.edges[i].list_offset + j;
					int k = scene.diffr_edge_graph.diffr_edge_neighbors[neigh_k];
					std::cout << k << "{" << scene.diffr_edge_graph.edges[k].getAxis() << "}; ";
				}
				std::cout << "\n";
			}
			std::cout << scene.diffr_edge_graph.edges.size() << "; "
				<< scene.diffr_edge_graph.diffr_edge_neighbors.size() << "\n";

			CGAL::draw(scene.mesh);
		}
		else
		{
			scene.addDiffractionFlags();

			CGAL::draw(scene.mesh);
		}
	}

	// LOAD
	{
		std::string full_filename = "Z:/course_7/blender/" + filename + ".obj";

		SoundScene scene;
		scene.readMeshFromFile(full_filename);
		scene.mapMaterials();

		// подгружаем
		std::cout << "loadDiffractionEdges\n";
		scene.loadDiffractionEdges(filename);

		std::cout << scene.diffr_edge_graph.edges.size() << "\n";

		std::cout << scene.diffr_edge_graph.edges.size() << "; "
			<< scene.diffr_edge_graph.diffr_edge_neighbors.size() << "\n";

		CGAL::draw(scene.mesh);

		if (test_graph)
		{
			// подгружаем
			std::cout << "loadDiffractionEdgeGraph\n";
			scene.loadDiffractionEdgeGraph(filename);

			// для проверки
			for (int i = 0; i < scene.total_num_edges; i++)
			{
				std::cout << i << ":\t";
				for (int j = 0; j < scene.diffr_edge_graph.edges[i].num_neighbors; j++)
				{
					int neigh_k = scene.diffr_edge_graph.edges[i].list_offset + j;
					int k = scene.diffr_edge_graph.diffr_edge_neighbors[neigh_k];
					std::cout << k << "{" << scene.diffr_edge_graph.edges[k].getAxis() << "}; ";
				}
				std::cout << "\n";
			}

			CGAL::draw(scene.mesh);

			std::cout << "start_addDiffractionFlags\n";
			scene.addDiffractionFlags();
			std::cout << "end_addDiffractionFlags\n";

			// Draw!
			CGAL::draw(scene.mesh);
		}
		else
		{
			scene.addDiffractionFlags();

			CGAL::draw(scene.mesh);
		}
	}
	return 0;
}

int test_isPointInside(std::string filename)
{
	std::string full_filename = "Z:/course_7/blender/" + filename + ".obj";

	SoundScene scene;
	scene.readMeshFromFile(full_filename);
	scene.mapMaterials();

	bool load_bde = false;
	std::cout << "Load latest save? (y/n) ";
	std::string answer = "y";
	//std::cin >> answer;

	if (answer != "y")
	{
		// помечаем диффр рёбра
		scene.markDiffractionEdges();

		// объединяем вместе в одно большие
		scene.combineDiffractionEdges();
		// сохраняем
		scene.saveDiffractionEdges(filename);

		std::cout << scene.diffr_edge_graph.edges.size() << "\n";

		// Draw!
		// CGAL::draw(scene.mesh);

		// построение графа видимости рёбер
		std::cout << "buildDiffractionEdgeGraph\n";
		scene.buildDiffractionEdgeGraph();

		// сохраняем
		scene.saveDiffractionEdgeGraph(filename);
	}
	else
	{
		// подгружаем
		scene.loadDiffractionEdges(filename);
		scene.loadDiffractionEdgeGraph(filename);
	}

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



	return 0;
}

int test_flags(int argc, char* argv[])
{
	std::string scene_name;
	std::string answer = "y";
	double sx, sy, sz;
	double dx = 0.02, dy = 0.02, dz = 0.02;
	double x_low, x_high, y_low, y_high;
	unsigned list_plane;
	double plane_val;
	unsigned phong_flag;
	bool is_2d = false;
	SoundScene scene;
	size_t max_diffr_order;
	scene.getInputParams(scene_name, answer,
		sx, sy, sz,
		dx, dy,
		x_low, x_high,
		y_low, y_high,
		list_plane, plane_val, phong_flag, is_2d, max_diffr_order);
	scene.max_diffraction_order = max_diffr_order;
	scene.is_2d = is_2d;

	std::string full_filename = "Z:/course_7/blender/" + scene_name + ".obj";
	std::cout << full_filename << "\n";

	scene.setFrequencies(3);
	scene.readMeshFromFile(full_filename);
	scene.mapMaterials();

	bool load_bde = false;
	std::cout << "Load latest save? (y/n) " << answer << "\n";

	// помечаем диффр рёбра
	scene.markDiffractionEdges();

	// объединяем вместе в одно большие
	scene.combineDiffractionEdges();

	std::cout << scene.diffr_edge_graph.edges.size() << "\n";

	// Draw!
	CGAL::draw(scene.mesh);

	// построение графа видимости рёбер
	std::cout << "buildDiffractionEdgeGraph\n";

	std::cout << scene.diffr_edge_graph.edges.size() << "\n";

	// Draw!
	// CGAL::draw(scene.mesh);

	// построение графа видимости рёбер
	std::cout << "buildDiffractionEdgeGraph\n";
	scene.buildDiffractionEdgeGraph();

	// для проверки
	std::cout << scene.diffr_edge_graph.edges.size() << "; "
		<< scene.diffr_edge_graph.diffr_edge_neighbors.size() << "\n";

	scene.addDiffractionFlags();
	full_filename = "Z:/course_7/blender/" + scene_name + "_added_flags.obj";
	scene.saveMeshFromFile(full_filename);
	CGAL::draw(scene.mesh);
	return 0;
}

int edge_anim(int argc, char* argv[])
{
	std::string scene_name;
	std::string answer = "y";
	double sx, sy, sz;
	double dx = 0.02, dy = 0.02, dz = 0.02;
	double x_low, x_high, y_low, y_high;
	unsigned list_plane = 2;
	double plane_val = 0;
	unsigned phong_flag = false;
	unsigned num_frequencies = 1;
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

	std::string full_filename = "Z:/course_7/blender/" + scene_name + ".obj";
	std::cout << full_filename << "\n";

	scene.setFrequencies(num_frequencies);
	scene.readMeshFromFile(full_filename);
	scene.mapMaterials();

	bool load_bde = false;
	std::cout << "Load latest save? (y/n) " << answer << "\n";
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

	double distance = 20;

	std::cout << "start_addDiffractionFlags\n";
	scene.addDiffractionFlags();
	std::cout << "end_addDiffractionFlags\n";

	CGAL::draw(scene.mesh);

	// source
	std::cout << "source position = " << sx << ", " << sy << ", " << sz << "\n";
	scene.setSoundSource(PointSoundSource(Point(sx, sy, sz)));
	auto vd_source = scene.mesh.add_vertex(scene.sound_source.position);

	//scene.epsilon_h = 2*std::max(dx, dy);

	std::cout << "dx = " << dx;
	std::cout << ", dy = " << dy << "\n";
	std::cout << "x = " << x_low << ".." << x_high << "\n";
	std::cout << "y = " << y_low << ".." << y_high << "\n";
	size_t x_dots_number = std::floor((x_high - x_low) / dx);
	size_t y_dots_number = std::floor((y_high - y_low) / dy);

	std::cout << "x_dots_number = " << x_dots_number << "\n";
	std::cout << "y_dots_number = " << y_dots_number << "\n";

	size_t x_curr_idx = 0;
	size_t y_curr_idx = 0;

	double x_curr = x_low;
	double y_curr = y_low;

	double _start = x_curr; double _finish = x_high; int _full_len = 30; int _curr_len = -1; int _curr_pcnt = -1;


	for (x_curr_idx = 0; x_curr_idx < x_dots_number; x_curr_idx++)
	{
		x_curr = x_low + x_curr_idx * dx;
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
			std::cout << "(" << x_curr << ", " << y_curr << ", " << plane_val << ")\n";
			// listener
			double move_listener_epsilon = std::min(dx, dy) * 0.1;
			Point moved_curr_list_pos = Point(
				curr_list_pos.x() + randomUniformNextVal() * move_listener_epsilon,
				curr_list_pos.y() + randomUniformNextVal() * move_listener_epsilon,
				curr_list_pos.z() + randomUniformNextVal() * move_listener_epsilon);

			scene.setSoundListener(PointSoundListener(moved_curr_list_pos, x_curr_idx, y_curr_idx));

			auto vd_listener = scene.mesh.add_vertex(scene.sound_listener.position);

			// Draw!
			//CGAL::draw(scene.mesh);

			double frequency = scene.frequencies[0];

			//bool has_direct_path = scene.addDirectPath();
			for (size_t i = 0; i < 100; i++)
			{
				scene.propagateSound();
			}
			for (size_t i = 0; i < 20; i++)
			{
				scene.addSpecularPath(frequency);
			}

			scene.mesh.remove_vertex(vd_listener);
		}
	}

	return 0;
}
