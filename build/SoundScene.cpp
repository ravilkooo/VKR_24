#include "SoundScene.h"

// Для рисования
#include <CGAL/draw_surface_mesh.h>

SoundScene::SoundScene()
	: total_num_edges(0),
	diffr_edge_graph(DiffractionEdgeGraph())
{
}

bool SoundScene::readMeshFromFile(std::string filename)
{
    if (!CGAL::IO::read_polygon_mesh(filename, mesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return false;
    }
    face_tree = Tree(mesh.faces().first, mesh.faces().second, mesh);
    return true;
}

bool SoundScene::saveMeshFromFile(std::string filename)
{
	if (!CGAL::IO::write_polygon_mesh(filename, mesh))
	{
		std::cerr << "Invalid output file." << std::endl;
		return false;
	}
	return true;
}

//////////////////////////////////
// Материал
//////////////////////////////////

bool SoundScene::mapMaterials()
{
	bool created;
	boost::tie(faces_material_map, created) = mesh.add_property_map<face_descriptor, unsigned>("f:material", 1);
	if (!created)
		return false;

	materials_list = std::vector<SoundMaterial>();
	materials_list.push_back(SoundMaterial::SILENCE);
	materials_list.push_back(SoundMaterial::SILENCE);

	return true;
}

//////////////////////////////////
// Дифр рёбра
//////////////////////////////////

bool SoundScene::markDiffractionEdges()
{
	bool created;
	boost::tie(diffr_map, created) = mesh.add_property_map<halfedge_descriptor, bool>("h:diffr", false);
	//assert(!created);
	if (!created)
		return false;

	// 1. проходимся по треугольникам (fd - i)
	for (face_descriptor fd : mesh.faces())
	{

		// 2. проходимся по полурёбрам (hed) треугольника (fd)
		for (halfedge_descriptor hed : mesh.halfedges_around_face(mesh.halfedge(fd)))
		{
			// соседние треугольнки через данное ребро (hed/opp_hed)
			halfedge_descriptor opposite_hed = mesh.opposite(hed);
			face_descriptor neighbor_fd = mesh.face(opposite_hed);

			// нужно отсеять "бесконечную грань" (граничные рёбра)
			if (neighbor_fd == Mesh::null_face())
				continue;

			if (diffr_map[hed] || diffr_map[opposite_hed])
			{
				continue;
			}

			vertex_descriptor v1 = mesh.target(hed);
			vertex_descriptor v2 = mesh.target(opposite_hed);

			// Получаем вершины напротив рёбер
			vertex_descriptor free1 = mesh.target(mesh.next(hed));
			vertex_descriptor free2 = mesh.target(mesh.next(opposite_hed));

			Vector normal1 = getVectorDirection(CGAL::Polygon_mesh_processing::compute_face_normal(fd, mesh));
			Vector normal2 = getVectorDirection(CGAL::Polygon_mesh_processing::compute_face_normal(neighbor_fd, mesh));

			Vector vec_1to2 = mesh.point(free2) - mesh.point(free1);

			// Скипаем вогнутый угол
			if (normal1 * vec_1to2 > 0 || normal2 * vec_1to2 < 0)
				continue;

			if (normal1 * normal2 < diffraction_threshold)
			{
				diffr_map[hed] = true;
				diffr_map[opposite_hed] = true;
			}
		}
	}
	return true;
}

bool SoundScene::combineDiffractionEdges()
{
	// 0 - не окрашено
	// > 0 - окрашено
	bool created;
	boost::tie(halfedge_col_map, created) = mesh.add_property_map<halfedge_descriptor, unsigned>("h:color", 0);
	assert(created);
	if (!created)
		return false;

	unsigned curr_col = 1;
	for (auto ed : mesh.edges())
	{
		halfedge_descriptor hed = ed.halfedge();

		if (!diffr_map[hed])
			continue;

		if (halfedge_col_map[hed] > 0)
			continue;

		halfedge_col_map[hed] = curr_col;
		halfedge_col_map[mesh.opposite(hed)] = curr_col;

		vertex_descriptor vd1 = mesh.target(mesh.opposite(hed));
		vertex_descriptor vd2 = mesh.target(hed);

		std::vector<halfedge_descriptor> temp_hed_list;
		temp_hed_list.push_back(hed);

		BigDiffractionEdge curr_big_diffr_edge(mesh, vd1, vd2, hed, mesh.opposite(hed), curr_col);

		// думаю нормали нормированны
		Vector normal1_curr = getVectorDirection(CGAL::Polygon_mesh_processing::compute_face_normal(mesh.face(hed), mesh));
		Vector normal2_curr = getVectorDirection(CGAL::Polygon_mesh_processing::compute_face_normal(mesh.face(mesh.opposite(hed)), mesh));
		double wedge_curr = std::acos(-1) - std::acos(normal1_curr * normal2_curr);

		Vector normal1_sum = normal1_curr;
		Vector normal2_sum = normal2_curr;
		double wedge_sum = wedge_curr;

		// продвигаемся вперёд по рёбрам
		halfedge_descriptor hed_prev = hed;
		vertex_descriptor vd_prev = vd1;
		vertex_descriptor vd_curr = vd2;
		bool move_next = true;
		Vector last_vector = getVectorDirection(mesh.point(vd2) - mesh.point(vd1));

		// vd_prev           vd_curr             vd_next
		//   O____>hed_prev>____O____>hed_next>____O

		while (move_next)
		{
			move_next = false;

			CGAL::Halfedge_around_target_circulator<Mesh> hebegin(mesh.halfedge(vd_curr), mesh), done(hebegin);
			do
			{
				halfedge_descriptor hed_next = *hebegin;
				hebegin++;

				// скипаем не диффр рёбра
				if (!diffr_map[hed_next])
					continue;

				// скипаем ребра того же цвета (на вский случай)
				if (halfedge_col_map[hed_next] == curr_col)
					continue;

				// проверяем что vd_curr указывает на начало hed_next
				if (mesh.target(mesh.opposite(hed_next)) != vd_curr)
					hed_next = mesh.opposite(hed_next);
				assert(mesh.target(mesh.opposite(hed_next)) == vd_curr);

				vertex_descriptor vd_next = mesh.target(hed_next);

				// скипаем обратное ребро
				if (vd_next == vd_prev)
					continue;

				// скипаем новое ребро, если образует большой угол с последним
				Vector new_vector = getVectorDirection(mesh.point(vd_next) - mesh.point(vd_curr));
				if (new_vector * last_vector < big_diffraction_edge_cos_angle_threshold)
					continue;

				// проверка совпадения клинов
				Vector normal1_next = getVectorDirection(CGAL::Polygon_mesh_processing::compute_face_normal(mesh.face(hed_next), mesh));
				Vector normal2_next = getVectorDirection(CGAL::Polygon_mesh_processing::compute_face_normal(mesh.face(mesh.opposite(hed_next)), mesh));
				double wedge_next = std::acos(-1) - std::acos(normal1_next * normal2_next);
				// скипаем если угол клина отличается
				if (std::abs(wedge_next - wedge_curr) > diffaction_wedge_variation_threshold)
					continue;
				// скипаем если ориентация клина отличается
				if (normal1_next * normal1_curr < diffaction_cos_face_normal_variation_threshold)
				{
					Vector _temp = Vector(normal2_next);
					normal2_next = Vector(normal1_next);
					normal1_next = Vector(_temp);
				}
				if (normal1_next * normal1_curr < diffaction_cos_face_normal_variation_threshold ||
					normal2_next * normal2_curr < diffaction_cos_face_normal_variation_threshold)
					continue;

				temp_hed_list.push_back(hed_next);
				curr_big_diffr_edge.vd_last = vd_next;
				curr_big_diffr_edge.hed_last = mesh.opposite(hed_next);
				halfedge_col_map[hed_next] = curr_col;
				halfedge_col_map[mesh.opposite(hed_next)] = curr_col;

				vd_prev = vd_curr;
				vd_curr = vd_next;
				hed_prev = hed_next;
				last_vector = new_vector;

				normal1_sum += normal1_next;
				normal2_sum += normal2_next;
				wedge_sum += wedge_next;

				move_next = true;

				break;

			} while (hebegin != done);
		}

		// продвигаемся назад по рёбрам
		hed_prev = mesh.opposite(hed);
		vd_prev = vd2;
		vd_curr = vd1;
		move_next = true;
		last_vector = getVectorDirection(mesh.point(vd1) - mesh.point(vd2));

		// vd_prev           vd_curr             vd_next
		//   O____>hed_prev>____O____>hed_next>____O

		while (move_next)
		{
			move_next = false;

			CGAL::Halfedge_around_target_circulator<Mesh> hebegin(mesh.halfedge(vd_curr), mesh), done(hebegin);
			do
			{
				halfedge_descriptor hed_next = *hebegin;
				hebegin++;

				// скипаем не диффр рёбра
				if (!diffr_map[hed_next])
					continue;

				// скипаем ребра того же цвета (на вский случай)
				if (halfedge_col_map[hed_next] == curr_col)
					continue;

				// проверяем что vd_curr указывает на начало hed_next
				if (mesh.target(mesh.opposite(hed_next)) != vd_curr)
					hed_next = mesh.opposite(hed_next);
				assert(mesh.target(mesh.opposite(hed_next)) == vd_curr);

				vertex_descriptor vd_next = mesh.target(hed_next);

				// скипаем обратное ребро
				if (vd_next == vd_prev)
					continue;

				// скипаем новое ребро, если образует большой угол с последним
				Vector new_vector = getVectorDirection(mesh.point(vd_next) - mesh.point(vd_curr));
				if (new_vector * last_vector < big_diffraction_edge_cos_angle_threshold)
					continue;

				// проверка совпадения клинов
				Vector normal1_next = getVectorDirection(CGAL::Polygon_mesh_processing::compute_face_normal(mesh.face(hed_next), mesh));
				Vector normal2_next = getVectorDirection(CGAL::Polygon_mesh_processing::compute_face_normal(mesh.face(mesh.opposite(hed_next)), mesh));
				double wedge_next = std::acos(-1) - std::acos(normal1_next * normal2_next);
				// скипаем если угол клина отличается
				if (std::abs(wedge_next - wedge_curr) > diffaction_wedge_variation_threshold)
					continue;
				// скипаем если ориентация клина отличается
				if (normal1_next * normal1_curr < diffaction_cos_face_normal_variation_threshold)
				{
					Vector _temp = Vector(normal2_next);
					normal2_next = Vector(normal1_next);
					normal1_next = Vector(_temp);
				}
				if (normal1_next * normal1_curr < diffaction_cos_face_normal_variation_threshold ||
					normal2_next * normal2_curr < diffaction_cos_face_normal_variation_threshold)
					continue;

				temp_hed_list.push_back(hed_next);
				curr_big_diffr_edge.vd_first = vd_next;
				curr_big_diffr_edge.hed_first = mesh.opposite(hed_next);
				halfedge_col_map[hed_next] = curr_col;
				halfedge_col_map[mesh.opposite(hed_next)] = curr_col;

				vd_prev = vd_curr;
				vd_curr = vd_next;
				hed_prev = hed_next;
				last_vector = new_vector;

				normal1_sum += normal1_next;
				normal2_sum += normal2_next;
				wedge_sum += wedge_next;

				move_next = true;

				break;

			} while (hebegin != done);
		}

		curr_big_diffr_edge.wedge = wedge_sum / temp_hed_list.size();
		curr_big_diffr_edge.normal1 = getVectorDirection(normal1_sum);
		curr_big_diffr_edge.normal2 = getVectorDirection(normal2_sum);
		curr_big_diffr_edge.flag_direction = getVectorDirection(curr_big_diffr_edge.normal1 + curr_big_diffr_edge.normal2);
		diffr_edge_graph.edges.push_back(curr_big_diffr_edge);

		// получается что curr_big_diffr_edge.col совпадает с индексом в diffr_edge_graph.edges
		curr_col++;
	}
	total_num_edges = diffr_edge_graph.edges.size();

	return true;
}

void SoundScene::buildDiffractionEdgeGraph()
{
	//std::cout << "total_num_edges = " << total_num_edges << "\n";
	for (size_t e1 = 0; e1 < total_num_edges; e1++)
	{
		//std::cout << "\ne1 = " << e1 << ": ";

		BigDiffractionEdge& edge1 = diffr_edge_graph.edges[e1];

		unsigned edge1_rays = clamp(edge1.getLength() / edge_resolution, min_rays_per_edge, max_rays_per_edge);

		Vector edge1_offset = edge1.flag_direction * edge_offset;
		Vector edge1_axis = edge1.getAxis();

		edge1.list_offset = diffr_edge_graph.diffr_edge_neighbors.size();

		for (size_t e2 = 0; e2 < total_num_edges; e2++)
		{
			//std::cout << "(e2 = " << e2 << ")";

			if (e1 == e2)
				continue;

			const BigDiffractionEdge& edge2 = diffr_edge_graph.edges[e2];

			if (!edge1.testEdgeOrientation(edge2, epsilon_f /*edge_offset*/))
				continue;

			const unsigned edge2_rays = clamp(edge2.getLength() / edge_resolution, min_rays_per_edge, max_rays_per_edge);

			Vector edge2_offset = edge2.flag_direction * edge_offset;
			Vector edge2_axis = edge2.getAxis();

			bool visible = false;

			for (unsigned i = 0; i < edge1_rays; i++)
			{
				// вычисляем точку на первом ребре
				Point p1 = edge1.getStart() + edge1_axis * ((i + 1.) / (edge1_rays + 1.));
				p1 += edge1_offset;

				// если точка не находится в области дифракции второго ребра, пропустить её
				if (!edge2.testOrientation(p1, 2*epsilon_f /*edge_offset*/))
					continue;

				for (unsigned j = 0; j < edge2_rays; j++)
				{
					// вычисляем точку на втором ребре
					Point p2 = edge2.getStart() + edge2_axis * ((j + 1.) / (edge2_rays + 1.));
					p2 += edge2_offset;

					// если точка не находится в области дифракции первого ребра, пропускаем её
					if (!edge1.testOrientation(p2, 2 * epsilon_f /*edge_offset*/))
						continue;

					// вычисляем луч и расстояние вдоль луча между точками
					double distance = getVectorLength(p2 - p1);

					// пропускаем этот луч, если точки совпадают
					if (distance < MATH_EPSILON)
						continue;

					Vector ray_direction = getVectorDirection(p2 - p1);
					Ray ray(p1, ray_direction);

					// пропускаем луч, если он направлен в противоположную сторону от нормалей обоих треугольников
					// это означает, что кандидат на ребре находится вне области дифракции первого ребра
					if ((ray_direction * edge1.normal1) < -ray_direction_threshold &&
						(ray_direction * edge1.normal2) < -ray_direction_threshold)
						continue;

					// прослеживаем луч. Если он не пересекается ни с чем, ребра взаимно видимы
					Skip skip;
					Skip skip;
					Ray_intersection intersection = face_tree.first_intersection(ray, skip);

					if (!intersection)
					{
						visible = true;
						break;
					}
				}
			}

			// Если ребро видно, добавляем его в список смежности
			if (visible)
				diffr_edge_graph.diffr_edge_neighbors.push_back(e2);
		}
		// Устанавливаем количество найденных соседей для этого ребра
		edge1.num_neighbors = (unsigned)(diffr_edge_graph.diffr_edge_neighbors.size() - edge1.list_offset);
	}
}

//////////////////////////////////
// Источник/приёмник
//////////////////////////////////

void SoundScene::setSoundSource(PointSoundSource new_sound_source)
{
	sound_source = new_sound_source;
}

void SoundScene::setSoundListener(PointSoundListener new_sound_listener)
{
	sound_listener = new_sound_listener;
}

Vector SoundScene::listenerToSourceDirection()
{
	return sound_source.position - sound_listener.position;
}

double SoundScene::pdfSound(const SoundPhongPathPoint& curr, 
	const SoundPhongPathPoint& next)
{
	assert(curr.type == SoundPhongPathPoint::PointType::SOURCE);
	Vector w = next.position - curr.position;
	double inv_dist2 = 1 / w.squared_length();
	w *= std::sqrt(inv_dist2);
	
	double pdf_pos = 0, pdf_dir = uniformSpherePdf();
	double pdf = pdf_dir * inv_dist2;
	if (next.type == SoundPhongPathPoint::PointType::SURFACE)
		pdf *= std::abs(next.normal * w);
	return pdf;
}

double SoundScene::pdfListener(const SoundPhongPathPoint& curr,
	const SoundPhongPathPoint& next)
{
	// пусть будет так
	assert(curr.type == SoundPhongPathPoint::PointType::LISTENER);
	Vector w = next.position - curr.position;
	double inv_dist2 = 1 / w.squared_length();
	w *= std::sqrt(inv_dist2);

	double pdf_pos = 0, pdf_dir = uniformSpherePdf();
	double pdf = pdf_dir * inv_dist2;
	if (next.type == SoundPhongPathPoint::PointType::SURFACE)
		pdf *= std::abs(next.normal * w);
	return pdf;
}

//////////////////////////////////
// Распространение звука
//////////////////////////////////

double SoundScene::getDistanceAttenuation(double distance)
{
	double broadband_attenuation = 1. / ((4. * Pi) * (1. + distance * distance));
	return broadband_attenuation;
}

void SoundScene::propagateSound()
{
	Vector ray_direction = sound_listener.sample_ray(is_2d);

	//std::cout << "\n" << ray_direction << ": ";
	Ray ray = Ray(sound_listener.position, ray_direction);
	Skip no_skip;
	Ray_intersection intersection = face_tree.first_intersection(ray, no_skip);
	if (!intersection)
		return;
	if (!boost::get<Point>(&(intersection->first)))
		return;
	face_descriptor fd = intersection->second;
	//std::cout << "faces_flag_map[" << fd << "] = " << faces_flag_map[fd] << ";";

	SkipFew skip(fd);
	// здесь внутри блока кода реализован алгоритм нахождения флагов за уже найденными флагами
	while (faces_flag_map[fd] > 0)
	{
		addDiffractionPath(ray_direction, intersection);
		if (CGAL_DRAW_MESH)
			CGAL::draw(mesh);
		intersection = face_tree.first_intersection(ray, skip);
		if (!intersection)
			break;
		if (!boost::get<Point>(&(intersection->first)))
			break;
		fd = intersection->second;
		skip.add(fd);
	}
}

//////////////////////////////////
// Прямые пути
//////////////////////////////////

bool SoundScene::addDirectPath()
{
	if (!isVisible(sound_listener.position, sound_source.position))
		return false;

	double dist = getVectorLength(sound_listener.position - sound_source.position);
	double energy = getDistanceAttenuation(dist);

	sound_listener.received_direct_IRs_archived.push_back(std::make_pair(dist / speed_of_sound, energy));

	return true;
}

//////////////////////////////////
// Пути модели Фонга (дифузия + зерк отраж)
//////////////////////////////////

bool SoundScene::addSpecularPath(double frequency)
{

	int max_depth = 8;

	double source_I = 1;
	double source_pdf_pos = 1;
	double source_pdf_dir = uniformSpherePdf();
	double source_Le = source_I; // В моём случае

	double listener_We = 1; // всегда единица?
	double listener_pdf_pos = 1;
	double listener_pdf_dir = uniformSpherePdf();

	// bounces - количество отклонений луча.

	// listener propagation
	Ray first_l_ray = 
		Ray(sound_source.position, uniformSampleSphere(randomUniformNextVal(), randomUniformNextVal()));
	std::vector<SoundPhongPathPoint> listener_path;
	{
		listener_path.push_back(SoundPhongPathPoint(first_l_ray.source()));
		listener_path[0].distance = 0;
		listener_path[0].frequency = frequency;
		listener_path[0].pdf_fwd = listener_pdf_pos;
		listener_path[0].attenuation = 1;
		listener_path[0].type = SoundPhongPathPoint::PointType::LISTENER;
		listener_path[0].beta = listener_We;
	}

	double listener_beta = listener_path[0].beta;
	int listener_bounces = randomWalk(first_l_ray, listener_path, max_depth, listener_pdf_dir, listener_beta, frequency);

	// source propagation
	Ray first_s_ray =
		Ray(sound_listener.position, uniformSampleSphere(randomUniformNextVal(), randomUniformNextVal()));
	std::vector<SoundPhongPathPoint> source_path;
	{
		source_path.push_back(SoundPhongPathPoint(first_s_ray.source()));
		source_path[0].distance = 0;
		source_path[0].frequency = frequency;
		source_path[0].pdf_fwd = source_pdf_pos;
		source_path[0].attenuation = 1;
		source_path[0].type = SoundPhongPathPoint::PointType::SOURCE;
		source_path[0].beta = source_Le;
	}
	double source_beta = source_Le;
	int source_bounces = randomWalk(first_s_ray, source_path, max_depth, source_pdf_dir, source_beta, frequency);

	// Eric Veach. 1997. Robust Monte Carlo methods for light transport simulation
	// M. Pharr, G. Humphreys, J. Wenzel. 2016. Physically Based Rendering. From Theory To Implementation. Third Edition

	int n_l = source_bounces + 1; // количество точек в пути от источника
	int n_e = listener_bounces + 1; // количество точек в пути от приёмника

	double C_sum = 0; // итоговый вклад
	double dist_sum = 0;

	// не рассматриваем пути, где s1 == 0 ( нет точек от от источника - луч от приёмника сам попал на источник )
	int s1_low_lim = 1;
	for (int st_sum = 1 + s1_low_lim; st_sum <= n_l + n_e; st_sum++)
	{
		int s1; // s1 - количество рассматриваемых вершин в подпути от источника
		int t1; // t1 - количество рассматриваемых вершин в подпути от приёмника
		// Рассматриваем поути длиной len = st_sum = s1 + t1
		for (s1 = s1_low_lim; s1 <= n_l; s1++)
		{
			// мы рассматриваем путь:
			// (x[0], x[1], ..., x[s1 - 1], x[s1], ..., x[k - 1], x[k]) =
			// (x[0], x[1], ..., x[s1 - 1], x[s1], ..., x[s1 + t1 - 2], x[s1 + t1 - 1]) :=
			// := (y[0], y[1], ..., y[s1 - 1]) U (z[t1 - 1], ..., z[1], z[0])
			// k = s1 + t1 - 1

			// вклад данной комбинации подпутей
			double C_curr = 0;
			// длина данной комбинации подпутей
			double dist_curr = 0;

			t1 = st_sum - s1; // количество вершин в пути от приёмника
			
			// отбросить прямые пути (когда добавлю отдельный метод для их учёта)
			if (s1 == 1 && t1 == 1)
				continue;

			// отбрасываем пути, где в подпути приёмника ноль вершин
			// отбрасываем пути, где в подпути приёмника вершин больше, чем их сгенерировалось по факту
			if (t1 < 1 || t1 > n_e)
			{
				continue;
			}
			// отбрасываем пути, где конец подпути от источника не виден из конца подпути от приёмника
			if (s1 != 0 && !isVisible(source_path[s1 - 1].position, listener_path[t1 - 1].position))
			{
				continue;
			}
			// отбрасываем пути, где концы подпутей от источника и от приёмника лежат на плоскостях, смотрящих в одну сторону
			// (они не могут быть соединены)
			if (s1 != 0 && source_path[s1 - 1].type == SoundPhongPathPoint::PointType::SURFACE
				&& listener_path[t1 - 1].type == SoundPhongPathPoint::PointType::SURFACE
				&& (source_path[s1 - 1].normal * listener_path[t1 - 1].normal) > DIRECTION_EQUALITY_THRESHOLD)
			{
				continue;
			}
			// отбрасываем случай, когда подпуть от источника не содержит вершин, а приёмник не попал в источник
			// (т.к. этот случай рассматривается, только когда приёмник попал в исчточник)
			if (s1 == 0 && listener_path[t1 - 1].type != SoundPhongPathPoint::PointType::SOURCE)
			{
				continue;
			}
			// отбрасываем случай, когда приёмник попал в источник, а подпуть от источника содержит вернишы
			// (т.к. этот случай рассматривается уже при s1==0 (подпуть источника не содержит вершин) )
			if (s1 != 0 && listener_path[t1 - 1].type == SoundPhongPathPoint::PointType::SOURCE)
			{
				continue;
			}
			
			// суммарная дистанция
			dist_curr = listener_path[t1 - 1].distance;
			if (s1 != 0)
				dist_curr += getVectorLength(source_path[s1 - 1].position - listener_path[t1 - 1].position)
				+ source_path[s1 - 1].distance;

			double L = 1;
			double mis_weight = 1;

			if (OUTPUT_MIS)
			{
				std::cout << "s1 = " << s1 << "; t1 = " << t1 << "\n";
			}
			

			if (MIS_ENABLED)
			{
				if (OUTPUT_MIS)
				{
					std::cout << "Connecting...\n";
				}
				if (t1 == 1)
				{
					const SoundPhongPathPoint& qs = source_path[s1 - 1];
					{
						Vector wi = listener_path[0].position - qs.position;
						double pdf = 1;
						double source_weight = source_I / wi.squared_length();
						wi = getVectorDirection(wi);
						SoundPhongPathPoint sampled = SoundPhongPathPoint(listener_path[0].position);
						sampled.beta = source_weight / pdf;
						sampled.pdf_fwd = 0;
						// случай s1 > 1. ост отбросили выше
						L = qs.beta * qs.calcBRDF(
							source_path[s1 - 2].position,
							sampled.position) * sampled.beta;
						if (qs.type == SoundPhongPathPoint::PointType::SURFACE) L *= std::abs(wi * qs.normal);
					}
				}
				else if (s1 == 1)
				{
					const SoundPhongPathPoint& pt = listener_path[t1 - 1];
					{
						Vector wi = source_path[0].position - pt.position;
						double pdf = 1;
						double source_weight = source_I / wi.squared_length();
						wi = getVectorDirection(wi);
						SoundPhongPathPoint sampled = SoundPhongPathPoint(source_path[0].position);
						sampled.beta = source_weight / pdf;
						sampled.pdf_fwd = 0;
						// случай t1 > 1 отбросили выше
						L = pt.beta * pt.calcBRDF(
							listener_path[t1 - 2].position,
							sampled.position)
							* sampled.beta;
						if (pt.type == SoundPhongPathPoint::PointType::SURFACE) L *= std::abs(wi * pt.normal);
					}
				}
				else
				{
					const SoundPhongPathPoint& qs = source_path[s1 - 1], & pt = listener_path[t1 - 1];
					{
						L = qs.beta * qs.calcBRDF(source_path[s1 - 2], pt) * pt.calcBRDF(listener_path[t1 - 2], qs) * pt.beta;
						L *= qs.geom(pt);
					}
				}
				if (OUTPUT_MIS)
				{
					std::cout << "L = " << L << "\n";
					std::cout << "in_calcMISWeight\n";
				}
				mis_weight = calcMISWeight(source_path, listener_path, s1, t1);
				if (OUTPUT_MIS)
				{
					std::cout << "mis_weight = " << mis_weight << "\n";
				}
				assert(!std::isnan(mis_weight));
				assert(!std::isnan(L));
			}

			C_curr = mis_weight * L;

			if (C_curr < 0)
			{

				std::cout << "UJAS!!!!! C_curr < 0: ";
				std::cout << "s1 = " << s1 << " ";
				std::cout << "t1 = " << t1 << "\n";

				continue;
			}

			double t_curr = dist_curr / speed_of_sound;

			if (t_curr < max_IR_length)
			{
				if (OUTPUT_MIS)
				{
					std::cout << "C_curr =\t" << C_curr << "\n";
					std::cout << "C_curr* =\t" << C_curr * getDistanceAttenuation(dist_curr) << "\n";
					std::cout << "s1 =\t" << s1 << " ";
					std::cout << "t1 =\t" << t1 << "\n";
				}

				sound_listener.received_phong_IRs_archived.push_back(std::make_pair(dist_curr / speed_of_sound, C_curr) );

			}
		}

	}
	if (OUTPUT_MIS)
	{
		if (listener_bounces >= 1 && source_bounces >= 1)
		{
			std::cout << "listener_bounces = " << listener_bounces << "; source_bounces = " << source_bounces << "\n";
			/*if (DRAW_FINAL_PHONG_PATH)
			{
				CGAL::draw(mesh);
			}*/
		}
	}

	return true;
}

int SoundScene::randomWalk(Ray first_ray, std::vector<SoundPhongPathPoint>& path,
	int max_depth, double pdf, double beta, double frequency)
{
	assert(!path.empty());
	if (path.empty()) // проверка на дурака (по сути никогда не должен сюда заходить)
	{
		path.push_back(SoundPhongPathPoint(first_ray.source()));
		path[0].distance = 0;
		path[0].pdf_fwd = 1;
		path[0].attenuation = 1;
		path[0].type = SoundPhongPathPoint::PointType::UNDEFINED;
	}

	double max_distance = max_IR_length * speed_of_sound;

	Ray ray(first_ray);
	Vector ray_dir = getVectorDirection(first_ray.to_vector());

	int depth = 0;
	double total_distance = 0;
	double attenuation = 1;
	double pdf_fwd = pdf;

	FlagSkip skip(*this, Mesh::null_face(), true);

	while (true)
	{
		Ray_intersection intersection = face_tree.first_intersection(ray, skip);
		bool has_isect = ((intersection) && (boost::get<Point>(&(intersection->first))));

		if (has_isect)
		{
			Point intrersection_pos = Point(*boost::get<Point>(&(intersection->first)));

			double intersection_distance = getVectorLength(intrersection_pos - ray.source());
			face_descriptor intersection_fd = intersection->second;

			auto face_normal =
				getVectorDirection(CGAL::Polygon_mesh_processing::compute_face_normal(intersection_fd, mesh));

			// Вычисляем скалярное произведение нормали треугольника с направлением входящего луча
			double ray_dot_normal = ray_dir * face_normal;

			// Инвертируем нормаль, если она направлена в ту же сторону, что и луч
			if (ray_dot_normal > 0.)
			{
				face_normal = -face_normal;
				ray_dot_normal = -ray_dot_normal;
			}

			// Смещаем точку пересечения на небольшое расстояние, чтобы избежать проблем с точностью плавающей запятой
			intrersection_pos += face_normal * ray_offset;

			// Суммируем общее расстояние вдоль пути
			total_distance += intersection_distance;

			// Если общее расстояние превышает максимальное, останавливаем этот луч
			if (total_distance > max_distance)
				break;

			// Получаем свойства материала для пересеченного треугольника
			SoundMaterial* material = &materials_list[faces_material_map[intersection_fd]];

			path.push_back(SoundPhongPathPoint(intrersection_pos));
			SoundPhongPathPoint& new_point = path.back();
			SoundPhongPathPoint& prev_point = path[path.size() - 2];
			new_point.type = SoundPhongPathPoint::PointType::SURFACE;
			new_point.frequency = frequency;
			new_point.normal = face_normal;
			new_point.material = material;
			new_point.distance = prev_point.distance + intersection_distance;
			new_point.refl_type = SoundPhongPathPoint::ReflectionType::NOCONTRIBUTION;
			new_point.pdf_fwd = prev_point.convertPdf(pdf_fwd, new_point);
			new_point.attenuation = prev_point.attenuation * pdf_fwd;
			new_point.beta = beta;

			attenuation *= pdf_fwd;

			if (++depth > max_depth)
				break;

			double pdf_rev = 0, brdf_fwd = 0, brdf_rev = 0;
			SoundPhongPathPoint::ReflectionType refl_type;

			// Вычисляем новый отраженный луч с использованием BRDF материала
			Vector out_dir = new_point.getPhongReflection(-ray_dir, refl_type,
				pdf_fwd, pdf_rev, brdf_fwd, brdf_rev);

			ray = Ray(intrersection_pos, out_dir);
			ray_dir = out_dir;
			assert(!(refl_type != SoundPhongPathPoint::ReflectionType::NOCONTRIBUTION && brdf_fwd <= 0));
			assert(!(refl_type != SoundPhongPathPoint::ReflectionType::NOCONTRIBUTION && pdf_fwd <= 0));
			new_point.refl_type = refl_type;
			new_point.brdf_fwd = brdf_fwd;
			new_point.brdf_rev = brdf_rev;


			if (refl_type == SoundPhongPathPoint::ReflectionType::NOCONTRIBUTION)
			{
				break;
			}
			prev_point.pdf_rev = new_point.convertPdf(pdf_rev, prev_point);
			beta *= brdf_fwd * std::abs(out_dir * face_normal) / pdf_fwd;
		}
		else
		{
			break;
		}
	}

	return depth;

}

double SoundScene::calcMISWeight(
	std::vector<SoundPhongPathPoint> source_path,
	std::vector<SoundPhongPathPoint> listener_path,
	int s, int t)
{
	if (s + t == 2) return 1.;
	double sum_ri = 0;
	auto remap0 = [](double f) -> double { return f != 0 ? f : 1; };

	SoundPhongPathPoint
		* qs = s > 0 ? &source_path[s - 1] : nullptr,
		* pt = t > 0 ? &listener_path[t - 1] : nullptr,
		* qs_minus = s > 1 ? &source_path[s - 2] : nullptr,
		* pt_minus = t > 1 ? &listener_path[t - 2] : nullptr;

	TemporaryOverride<double> a4;
	if (pt)
	{
		if (s > 1)
			a4 = { &pt->pdf_rev, qs->convertPdf(qs->calcPhongPdf(*qs_minus, *pt), *pt) };
		else
			a4 = { &pt->pdf_rev, pdfSound(*qs, *pt) };
	}

	TemporaryOverride<double> a5;
	if (pt_minus)
		a5 = { &pt_minus->pdf_rev, pt->convertPdf(pt->calcPhongPdf(*qs, *pt_minus), *pt_minus)  };

	TemporaryOverride<double> a6;
	if (qs)
		if (t > 1)
			a6 = { &qs->pdf_rev, pt->convertPdf(pt->calcPhongPdf(*pt_minus, *qs), *qs) };
		else
			a6 = { &qs->pdf_rev, pdfListener(*pt, *qs) };
	TemporaryOverride<double> a7;
	if (qs_minus)
		a7 = { &qs_minus->pdf_rev, pt->convertPdf(qs->calcPhongPdf(*pt, *qs_minus), *qs_minus) };

	double ri = 1;
	for (int i = t - 1; i > 0; --i)
	{
		ri *= remap0(listener_path[i].pdf_rev) / remap0(listener_path[i].pdf_fwd);
		sum_ri += ri;
	}

	ri = 1;
	for (int i = s - 1; i >= 0; --i)
	{
		ri *= remap0(source_path[i].pdf_rev) / remap0(source_path[i].pdf_fwd);
		sum_ri += ri;
	}
	return 1 / (1 + sum_ri);
}

//////////////////////////////////
// Дифр пути
//////////////////////////////////

bool SoundScene::addDiffractionPath(Vector ray_direction, Ray_intersection intersection)
{

	listener_to_source_direction = getVectorDirection(sound_source.position - sound_listener.position);

	ImpulseResponse curr_IR(num_frequencies);

	face_descriptor fd = intersection->second;
	unsigned col = faces_flag_map[fd];
	const BigDiffractionEdge& diffr_edge = diffr_edge_graph.edges[col - 1];

	if (!diffr_edge.testOrientation(sound_listener.position, /*edge_offset*/ /*epsilon_h*/ epsilon_h))
	{
		return false;
	}

	const Vector edge_axis = diffr_edge.getAxis();
	const Point v_first = diffr_edge.getStart();

	double edge_t = 0.;
	if (std::abs(listenerToSourceDirection() * edge_axis) < 0.999)
		edge_t = computePointOfClosestApproach(v_first, edge_axis,
			sound_source.position, listenerToSourceDirection());
	else
		edge_t = 0.5 * diffr_edge.getLength();
	
	if (edge_t < 0. || edge_t > diffr_edge.getLength())
		return false;

	Point next_listener_image_position = v_first + edge_axis * edge_t;
	next_listener_image_position += diffr_edge.flag_direction * ray_offset;

	curr_IR.diffraction_points_list.push_back(SoundDiffractionPathPoint(sound_listener.position));
	curr_IR.diffraction_points_list.push_back(SoundDiffractionPathPoint(next_listener_image_position, diffr_edge.col, diffr_edge));
	curr_IR.diffraction_points_list[curr_IR.diffraction_points_list.size() - 1].distance = getVectorLength(next_listener_image_position - sound_listener.position);
	curr_IR.last_valid_index = 0;
	
	curr_IR.point_responses.push_back(FrequencyBandResponse(num_frequencies));
	lowest_frequency = curr_IR.point_responses[0].frequencies[0];
	highest_frequency = curr_IR.point_responses[0].frequencies[curr_IR.point_responses[0].frequencies.size() - 1];
	return recursiveDiffraction(curr_IR, diffr_edge, 1);

	return true;
}

bool SoundScene::recursiveDiffraction(ImpulseResponse& curr_IR, const BigDiffractionEdge& diffr_edge, unsigned depth)
{
	const Point& source_position = sound_source.position;
	SoundDiffractionPathPoint& last_path_point = curr_IR.diffraction_points_list[curr_IR.diffraction_points_list.size() - 2];
	SoundDiffractionPathPoint& this_path_point = curr_IR.diffraction_points_list[curr_IR.diffraction_points_list.size() - 1];
	const Point last_listener_image_position = last_path_point.position;
	const Point listener_image_position = this_path_point.position;

	// Определяем, на какой стороне ребра находится приёмник
	double plane1_distance = getSignedDistanceTo(diffr_edge.normal1, diffr_edge.getStart(), last_listener_image_position);
	double plane2_distance = getSignedDistanceTo(diffr_edge.normal2, diffr_edge.getStart(), last_listener_image_position);

	bool listener_orientation = plane1_distance > plane2_distance && plane1_distance > 0.;

	// Получаем свободную вершину треугольника на той же стороне ребра, что и приёмник
	Point triangle_free_vertex = listener_orientation ?
		diffr_edge.getFreeVertex1() : diffr_edge.getFreeVertex2();

	// Вычисляем плоскость границы теневой области для ребра и последней позиции источника изображения
	Plane shadow_boundary(last_listener_image_position, diffr_edge.getStart(), diffr_edge.getStart() + diffr_edge.getAxis());

	// Убеждаемся, что нормаль плоскости границы тени указывает в правильном направлении
	if (getSignedDistanceTo(shadow_boundary, last_listener_image_position, triangle_free_vertex) < 0.)
		shadow_boundary = Plane(last_listener_image_position, diffr_edge.getStart() + diffr_edge.getAxis(), diffr_edge.getStart());

	// Получаем плоскость, которая обращена к приёмнику
	const Plane& listener_plane = listener_orientation ? diffr_edge.getPlane1() : diffr_edge.getPlane2();
	this_path_point.listener_plane = listener_plane;

	// Получаем плоскость, которая определяет границу треугольника на противоположной стороне ребра
	const Plane& opposite_plane = listener_orientation ? diffr_edge.getPlane2() : diffr_edge.getPlane1();
	this_path_point.source_plane = opposite_plane;

	// Определяем, находится ли источник в теневой области, проверяя его относительно граничных плоскостей
	bool source_in_shadow_region = getSignedDistanceTo(shadow_boundary, last_listener_image_position, source_position) > /*-epsilon_h*/ 0 &&
		getSignedDistanceTo(opposite_plane, diffr_edge.getStart(), source_position) > /*-epsilon_h*/ 0;

	if (depth == 1 && alreadyContainsPath(diffr_edge.col))
	{
		return false;
	}

	if (source_in_shadow_region)
	{
		bool valid = true;
		size_t point_index;

		const unsigned last_point_index = curr_IR.diffraction_points_list.size() - 1;

		// Проверяем путь к приёмнику, который еще не был проверен, чтобы убедиться, что он действителен
		for (point_index = curr_IR.last_valid_index; point_index < last_point_index; point_index++)
		{
			assert(point_index < curr_IR.diffraction_points_list.size());
			SoundDiffractionPathPoint& last_point = curr_IR.diffraction_points_list[point_index];
			assert(point_index + 1 < curr_IR.diffraction_points_list.size());
			SoundDiffractionPathPoint& this_point = curr_IR.diffraction_points_list[point_index + 1];

			// Вычисляем направление от последней позиции источника изображения до текущей
			Vector direction = this_point.position - last_point.position;
			double distance = getVectorLength(direction);

			if (distance > diffraction_epsilon)
				direction /= distance;
			else
			{
				valid = false;
				break;
			}
			
			Ray ray(last_point.position + direction * diffraction_epsilon, direction);
			FlagSkip skip(*this, Mesh::null_face(), true);
			
			Ray_intersection intersection = face_tree.first_intersection(ray, skip);
			if (intersection) {
				if (boost::get<Point>(&(intersection->first)))
				{
					Point p_between = Point(*boost::get<Point>(&(intersection->first)));
					double distance_to_intersection = getVectorLength(p_between - ray.source());
					if (distance_to_intersection < distance - diffraction_epsilon * 2)
					{
						valid = false;
						break;
					}
				}
					
			}

			// Суммируем общее расстояние вдоль пути.
			this_point.distance = last_point.distance + distance;

			/*
			// Вычисляем отклик пути для последнего ребра, если оно было
			if (point_index > 0)
			{
				const SoundDiffractionPathPoint& last_last_point = curr_IR.diffraction_points_list[point_index - 1];

				FrequencyBandResponse total_attenuation = FrequencyBandResponse(frequencies.size());

				if (point_index > 1)
				{
					total_attenuation = total_attenuation * curr_IR.point_responses[point_index - 2];
				}
					

				if (curr_IR.point_responses.size() <= point_index)
					curr_IR.point_responses.push_back(total_attenuation);
				else
					curr_IR.point_responses[point_index] = total_attenuation;

			}
			*/
		}

		curr_IR.last_valid_index = point_index;

		if (valid)
		{
			// Вычисляем вектор от текущей позиции источника изображения до приёмника
			Vector source_direction = source_position - listener_image_position;
			double source_distance = getVectorLength(source_direction);

			if (source_distance > diffraction_epsilon)
			{
				source_direction /= source_distance;

				// Проверяем путь от позиции изображения источника до приёмника, чтобы убедиться, что он свободен
				bool source_visible = true;
				Ray ray_si2l(
					listener_image_position + source_direction * diffraction_epsilon,
					source_direction);
				FlagSkip skip_si2l(*this, Mesh::null_face(), true);
				Ray_intersection intersection_si2l = face_tree.first_intersection(ray_si2l, skip_si2l);
				if (intersection_si2l)
				{
					if (boost::get<Point>(&(intersection_si2l->first)))
					{
						Point p_between_si2l = Point(*boost::get<Point>(&(intersection_si2l->first)));
						double distance_to_intersection = getVectorLength(p_between_si2l - ray_si2l.source());
						if (distance_to_intersection < source_distance - diffraction_epsilon * 2)
						{
							source_visible = false;
						}
					}
				}

				if (source_visible)
				{
					FrequencyBandResponse total_attenuation;
					if (depth == 1)
					{
						total_attenuation =
							computeUTDAttenuation(source_position, this_path_point, last_listener_image_position,
								speed_of_sound, curr_IR.fbr.frequencies, epsilon_h);
					}
					else
					{

						total_attenuation =
							computeUTDAttenuation_N_order(source_position, curr_IR.diffraction_points_list, sound_listener.position,
								speed_of_sound, curr_IR.fbr.frequencies, epsilon_h);
					}


					double total_distance = curr_IR.diffraction_points_list.back().distance + source_distance;

					sound_listener.received_diffraction_IRs.push_back(
						PathData(PathType::DIFFRACTION,
							getDistanceAttenuation(total_distance)* total_attenuation,
							total_distance / speed_of_sound, curr_IR.diffraction_points_list));

					if (DRAW_ALL_CORRECT_PATHES)
					{
						for (int ss = 0; ss < curr_IR.diffraction_points_list.size() - 1; ss++)
						{
							addArrow(curr_IR.diffraction_points_list[ss].position, curr_IR.diffraction_points_list[ss + 1].position);
						}
						addArrow(curr_IR.diffraction_points_list.back().position, source_position);
						CGAL::draw(mesh);
					}
				}
			}
		}


	}

	// Возвращаемся, если достигнута максимальная глубина
	if (depth >= max_diffraction_order || diffr_edge_graph.edges.size() == 0)
	{
		return false;
	}

	// Проверяем соседние ребра на наличие путей дифракции более высокого порядка
	const DiffractionEdgeGraph* graph = &diffr_edge_graph;
	const unsigned num_neighbors = diffr_edge.num_neighbors;
	const unsigned neighbor_list_start = diffr_edge.list_offset;
	const unsigned neighbor_list_end = neighbor_list_start + num_neighbors;

	for (unsigned n = neighbor_list_start; n < neighbor_list_end; n++)
	{
		// Рассматриваем соседнее ребро
		const BigDiffractionEdge* neighbor = graph->getEdgeNeighbor(n);

		// Вычисляем точку ближайшего приближения линии от приёмника к источнику и линии ребра
		double edge_t = 0.;
		if (std::abs(listener_to_source_direction * neighbor->getAxis()) < 0.999)
			edge_t = computePointOfClosestApproach(neighbor->getStart(), neighbor->getAxis(),
				sound_listener.position, listener_to_source_direction);
		else
			edge_t = 0.5 * neighbor->getLength();

		// Убеждаемся, что ближайшая точка лежит в пределах ребра. Пропускаем недопустимые ребра
		if (edge_t < 0. || edge_t > neighbor->getLength())
			continue;

		// Вычисляем ближайшую точку на ребре. Это новая позиция изображения слушателя
		Point next_listener_image_position = neighbor->getStart() + neighbor->getAxis() * edge_t;
		next_listener_image_position += neighbor->flag_direction * diffraction_epsilon;

		// Определяем, находится ли эта точка в предыдущем клине дифракции
		bool neighbor_in_wedge = getSignedDistanceTo(shadow_boundary, diffr_edge.getStart(), next_listener_image_position) > -epsilon_h &&
			getSignedDistanceTo(opposite_plane, diffr_edge.getStart(), next_listener_image_position) > -epsilon_h;


		// Рекурсивно находим пути дифракции для допустимых ребер
		if (neighbor_in_wedge)
		{
			// Добавляем следующую точку дифракции в массив дифр точек
			curr_IR.diffraction_points_list.push_back(SoundDiffractionPathPoint(next_listener_image_position,
				neighbor->col, neighbor));

			if (CGAL_DRAW_MESH)
				CGAL::draw(mesh);

			recursiveDiffraction(curr_IR, *neighbor, depth + 1);

			// Удаляем точку дифракции из массива дифр точек
			curr_IR.diffraction_points_list.pop_back();

			if (curr_IR.last_valid_index >= depth)
				curr_IR.last_valid_index = curr_IR.last_valid_index - 1;

		}
	}

	return true;
}

double SoundScene::computePointOfClosestApproach(Point p1, Vector v1, Point p2, Vector v2)
{
	double v1DotV2 = v1 * v2;
	Vector p1ToP2 = p2 - p1;
	return ((p1ToP2 * v1) - (p1ToP2 * v2) * v1DotV2) / (1. - v1DotV2 * v1DotV2);
}

bool SoundScene::alreadyContainsPath(unsigned diffr_edge_col)
{
	for (int i = 0; i < sound_listener.received_diffraction_IRs.size(); i++)
	{
		SoundDiffractionPathPoint first_point =
			sound_listener.received_diffraction_IRs[i].diffraction_points_list[1];

		if (first_point.diffr_edge_col == diffr_edge_col)
		{
			return true;
		}
	}
	return false;
}

//////////////////////////////////
// Флаги
//////////////////////////////////

bool SoundScene::addDiffractionFlags()
{
	bool created;
	boost::tie(faces_flag_map, created) = mesh.add_property_map<face_descriptor, unsigned>("f:diffr", 0);
	if (!created)
		return false;
	for (BigDiffractionEdge bde : diffr_edge_graph.edges)
	{
		Point _v1(mesh.point(bde.vd_first));
		Point _v2(mesh.point(bde.vd_last));
		Point _v3 = _v2 + d_flag * bde.flag_direction;
		Point _v4 = _v1 + d_flag * bde.flag_direction;
		auto _v1d = mesh.add_vertex(_v1); auto _v2d = mesh.add_vertex(_v2);
		auto _v3d = mesh.add_vertex(_v3); auto _v4d = mesh.add_vertex(_v4);

		face_descriptor fd1 = mesh.add_face(_v1d, _v2d, _v3d);
		faces_flag_map[fd1] = bde.col;

		face_descriptor fd2 = mesh.add_face(_v3d, _v4d, _v1d);
		faces_flag_map[fd2] = bde.col;

		bde.flag_face_descriptor = fd1;
	}
	face_tree = Tree(mesh.faces().first, mesh.faces().second, mesh);
	return true;
}

//////////////////////////////////
// Вспомогателььные методы
//////////////////////////////////

bool SoundScene::isVisible(const Point& p1, const Point& p2)
{
	Vector dir = p2 - p1;
	dir = getVectorDirection(dir);

	FlagSkip skip(*this, Mesh::null_face(), true);
	Ray_intersection intersection = face_tree.first_intersection(Ray(p1, dir), skip);
	if (intersection)
	{
		if (boost::get<Point>(&(intersection->first)))
		{
			Point intersection_pos = Point(*boost::get<Point>(&(intersection->first)));
			return getVectorLength(intersection_pos - p1) > getVectorLength(p2 - p1);
		}
		else
		{
			return true;
		}
	}
	else
	{
		return true;
	}
}
void SoundScene::setFrequencies(unsigned new_num_frequencies)
{
	num_frequencies = new_num_frequencies;
	frequencies.clear();
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
void SoundScene::getSceneBoundingBox(Point& min_point, Point& max_point)
{
	CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);

	min_point = Point(bbox.xmin(), bbox.ymin(), bbox.zmin());
	max_point = Point(bbox.xmax(), bbox.ymax(), bbox.zmax());
}
// Функция для проверки видимости обратной стороны
bool SoundScene::isBackSideVisible(const Point& point)
{
	double back_side_visible_epsilon = 0.0001;
	double random_amplitude = 0.001;
	for (auto fd : faces(mesh))
	{
		if (faces_flag_map[fd] > 0)
		{
			continue;
		}
		auto hd = mesh.halfedge(fd);
		std::vector<Point> p_tr;
		p_tr.push_back(mesh.point(mesh.source(hd)));
		p_tr.push_back(mesh.point(mesh.target(hd)));
		p_tr.push_back(mesh.point(mesh.target(mesh.next(hd))));

		Point centroid = CGAL::centroid(p_tr[0], p_tr[1], p_tr[2]);
		unsigned intersect_cnt = 0;

		// Создание луча из точки point к центру треугольника
		for (int i = 0; i < 3; i++)
		{
			Point moved_centroid = centroid + (p_tr[i] - centroid) * (random_amplitude * randomUniformNextVal() + back_side_visible_epsilon);
			Point moved_point = point + back_side_visible_epsilon * uniformSampleSphere(randomUniformNextVal(), randomUniformNextVal());
			Vector direction = getVectorDirection(moved_centroid - moved_point);

			Vector normal = CGAL::Polygon_mesh_processing::compute_face_normal(fd, mesh);
			if (normal * direction <= 0)
			{
				continue;
			}

			Ray ray(moved_point, direction);

			// Поиск пересечений луча с поверхностью
			FlagSkip skip(*this, Mesh::null_face(), true);
			auto intersection = face_tree.first_intersection(ray, skip);
			if (!intersection)
				continue;
			if (!boost::get<Point>(&(intersection->first)))
				continue;
			{
				const face_descriptor& intersected_fd = intersection->second;
				{
					normal = CGAL::Polygon_mesh_processing::compute_face_normal(intersected_fd, mesh);
					if (normal * direction > 0)
					{
						intersect_cnt++;// Точка видит обратную сторону поверхности
					}
				}
			}
		}

		if (intersect_cnt == 3)
		{
			return true;
		}
	}

	return false;
}

//////////////////////////////////
// Дополнительно для графики
//////////////////////////////////

void SoundScene::addArrow(Point start, Point end)
{
	//std::cout << start << "::" << end << "\n";

	auto _v1 = mesh.add_vertex(start);
	auto _v2 = mesh.add_vertex(end);
	e_buffer.push_back(mesh.add_edge(_v1, _v2));

	Vector dir = getVectorDirection(end - start);
	Vector orth_dir = getVectorDirection(Vector(dir.y(), -dir.x(), 0));

	if ((std::abs(dir.x()) < 0.00001) && (std::abs(dir.y()) < 0.00001))
		orth_dir = getVectorDirection(Vector(dir.z(), 0, -dir.x()));

	auto _v3 = mesh.add_vertex(end - 0.1 * dir + 0.1 * orth_dir);
	auto _v4 = mesh.add_vertex(end - 0.1 * dir - 0.1 * orth_dir);

	e_buffer.push_back(mesh.add_edge(_v3, _v2));
	e_buffer.push_back(mesh.add_edge(_v4, _v2));

	v_buffer.push_back(_v1);
	v_buffer.push_back(_v2);
	v_buffer.push_back(_v3);
	v_buffer.push_back(_v4);
}

void SoundScene::clearArrowBuffers()
{
	for (int i = 0; i < e_buffer.size(); i++)
	{
		auto _e = e_buffer[i];
		mesh.remove_edge(mesh.edge(_e));
	}
	for (int i = 0; i < v_buffer.size(); i++)
	{
		auto _v = v_buffer[i];
		mesh.remove_vertex(_v);
	}
	v_buffer = std::vector<vertex_descriptor>();
	e_buffer = std::vector<halfedge_descriptor>();
}

void SoundScene::removeLastArrow()
{
	for (int i = 0; i < 3; i++)
	{
		if (e_buffer.size() == 0)
			break;
		mesh.remove_edge(mesh.edge( e_buffer.back() ));
		e_buffer.pop_back();
	}
	for (int i = 0; i < 4; i++)
	{
		if (v_buffer.size() == 0)
			break;
		mesh.remove_vertex(v_buffer.back());
		v_buffer.pop_back();
	}
}

//////////////////////////////////
// Сохранение и загрузка
//////////////////////////////////

std::string SoundScene::readInputFilename()
{
	std::ifstream file("input_filename.txt");

	if (!file.is_open())
	{
		std::cerr << "Incorrect input filename." << std::endl;
		return "";
	}

	std::string filename;
	file >> filename; // Считываем слово из файла

	if (file.fail())
	{
		std::cerr << "Reading error." << std::endl;
		return "";
	}

	std::cout << "Scene : " << filename << std::endl;

	file.close();

	return filename;
}

std::string SoundScene::createOutputFilename(std::string input_filename)
{
	auto now = std::chrono::system_clock::now();
	auto now_time_t = std::chrono::system_clock::to_time_t(now);

	std::tm local_tm = *std::localtime(&now_time_t);

	std::stringstream ss;
	ss << std::setfill('0') << std::setw(2) << local_tm.tm_mday << std::setw(2) << local_tm.tm_mon + 1
		<< "_" << std::setw(2) << local_tm.tm_hour << std::setw(2) << local_tm.tm_min;

	return input_filename + "_" + ss.str();
}

void SoundScene::saveOutputFilename(const std::string& data)
{
	std::string filename = "output_filename.txt";
	std::ofstream file(filename);
	if (file.is_open())
	{
		file << data;
		file.close();
		std::cout << "Results' codename: " << data << std::endl;
		std::cout << "Results' codename saved in file: " << filename << std::endl;
	}
	else
	{
		std::cerr << "Can't open fiie " << filename << " for writing." << std::endl;
	}
}

void SoundScene::getInputParams(std::string& scene_name, std::string& answer,
	double& sx, double& sy, double& sz,
	double& dx, double& dy,
	double& x_low, double& x_high,
	double& y_low, double& y_high,
	unsigned& list_plane, double& plane_val,
	unsigned& phong_flag, bool& is_2d,
	size_t& max_diffr_order)
{
	std::ifstream file("input_filename.txt");
	file.imbue(std::locale("C"));

	if (!file.is_open())
	{
		std::cerr << "Incorrect input filename." << std::endl;
		return;
	}

	file >> scene_name; // Считываем слово из файла

	if (file.fail())
	{
		std::cerr << "Reading error." << std::endl;
		return;
	}

	std::cout << "Scene : " << scene_name << std::endl;

	file >> answer;
	file >> sx >> sy >> sz;
	file >> dx >> dy;
	file >> x_low >> x_high;
	file >> y_low >> y_high;
	std::string list_plane_axis = "z";
	list_plane = 2;
	file >> list_plane_axis;
	if (list_plane_axis == "x") list_plane = 0;
	if (list_plane_axis == "y") list_plane = 1;
	file >> plane_val;
	file >> phong_flag;
	file >> is_2d;
	file >> max_diffr_order;

	file.close();

	return;
}
void SoundScene::saveFrequencies(std::string filename)
{
	// Сохранение частот
	std::ofstream freq_out(filename + "_frequencies");
	if (freq_out.is_open())
	{
		size_t num_frequencies = frequencies.size();
		for (size_t i = 0; i < num_frequencies; i++)
		{
			freq_out << frequencies[i] << " ";
			std::cout << frequencies[i] << " ";

		}
		freq_out.close();
	}
}

bool SoundScene::saveDiffractionEdges(std::string filename)
{
	std::ofstream out;

	// save diffr_map
	std::cout << "save diffr_map \n";
	out = std::ofstream(filename + "_diffr_map", std::ios::binary);
	if (!out)
	{
		throw std::runtime_error("Cannot open file for writing diffr_map");
	}
	for (halfedge_descriptor hed : mesh.halfedges())
	{
		size_t key = hed;
		bool value = diffr_map[hed];
		//std::cout << key << "; " << value << "; " << hed << "\n";
		out.write(reinterpret_cast<const char*>(&key), sizeof(key));
		out.write(reinterpret_cast<const char*>(&value), sizeof(value));
	}
	out.close();

	// save halfedge_col_map
	std::cout << "save halfedge_col_map \n";
	out = std::ofstream(filename + "_halfedge_col_map", std::ios::binary);
	if (!out)
	{
		throw std::runtime_error("Cannot open file for writing halfedge_col_map");
	}
	for (halfedge_descriptor hed : mesh.halfedges())
	{
		size_t key = hed;
		unsigned value = halfedge_col_map[hed];
		out.write(reinterpret_cast<const char*>(&key), sizeof(key));
		out.write(reinterpret_cast<const char*>(&value), sizeof(value));
	}
	out.close();

	// save diffr_edge_graph_edges
	std::cout << "save diffr_edge_graph_edges \n";
	out = std::ofstream(filename + "_diffr_edge_graph_edges", std::ios::binary);
	if (!out)
	{
		throw std::runtime_error("Cannot open file for writing diffr_edge_graph_edges");
	}
	for (size_t e1 = 0; e1 < total_num_edges; e1++)
	{
		BigDiffractionEdge& edge1 = diffr_edge_graph.edges[e1];

		size_t _vd_first	= edge1.vd_first;
		size_t _vd_last		= edge1.vd_last;
		size_t _hed_first	= edge1.hed_first;
		size_t _hed_last	= edge1.hed_last;
		unsigned col = edge1.col;
		double wedge = edge1.wedge;
		std::vector<Vector> v_arr = { edge1.flag_direction, edge1.normal1, edge1.normal2 };

		out.write(reinterpret_cast<const char*>(&_vd_first), sizeof(_vd_first));
		out.write(reinterpret_cast<const char*>(&_vd_last), sizeof(_vd_last));
		out.write(reinterpret_cast<const char*>(&_hed_first), sizeof(_hed_first));
		out.write(reinterpret_cast<const char*>(&_hed_last), sizeof(_hed_last));
		out.write(reinterpret_cast<const char*>(&col), sizeof(col));
		out.write(reinterpret_cast<const char*>(&wedge), sizeof(wedge));

		for (size_t i = 0; i < v_arr.size(); i++)
		{
			double x = v_arr[i].x();
			double y = v_arr[i].y();
			double z = v_arr[i].z();

			out.write(reinterpret_cast<const char*>(&x), sizeof(x));
			out.write(reinterpret_cast<const char*>(&y), sizeof(y));
			out.write(reinterpret_cast<const char*>(&z), sizeof(z));
		}
	}
	out.close();
	return true;
}
bool SoundScene::saveDiffractionEdgeGraph(std::string filename)
{
	std::ofstream out;

	// save diffr_edge_graph
	std::cout << "save diffr_edge_graph \n";
	out = std::ofstream(filename + "_diffr_edge_graph", std::ios::binary);
	if (!out)
	{
		throw std::runtime_error("Cannot open file for writing diffr_edge_graph");
	}
	for (size_t e1 = 0; e1 < total_num_edges; e1++)
	{
		BigDiffractionEdge& edge1 = diffr_edge_graph.edges[e1];

		unsigned list_offset = edge1.list_offset;
		unsigned num_neighbors = edge1.num_neighbors;

		out.write(reinterpret_cast<const char*>(&list_offset), sizeof(list_offset));
		out.write(reinterpret_cast<const char*>(&num_neighbors), sizeof(num_neighbors));
	}
	for (size_t i = 0; i < diffr_edge_graph.diffr_edge_neighbors.size(); i++)
	{
		unsigned col = diffr_edge_graph.diffr_edge_neighbors[i];
		out.write(reinterpret_cast<const char*>(&col), sizeof(col));
	}

	out.close();

	return true;
}

bool SoundScene::loadDiffractionEdges(std::string filename)
{
	std::ifstream in;
	bool created;

	// load diffr_map
	std::cout << "load diffr_map \n";
	in = std::ifstream(filename + "_diffr_map", std::ios::binary);
	if (!in)
	{
		throw std::runtime_error("Cannot open file for reading diffr_map");
	}

	boost::tie(diffr_map, created) = mesh.add_property_map<halfedge_descriptor, bool>("h:diffr", false);

	while (!in.eof())
	{
		size_t key;
		bool value;
		in.read(reinterpret_cast<char*>(&key), sizeof(key));
		in.read(reinterpret_cast<char*>(&value), sizeof(value));
		if (in)
		{ // Проверка на успешное чтение
			halfedge_descriptor hed(key);
			diffr_map[hed] = value;
		}
	}
	in.close();

	// load halfedge_col_map
	std::cout << "load halfedge_col_map \n";
	in = std::ifstream(filename + "_halfedge_col_map", std::ios::binary);
	if (!in)
	{
		throw std::runtime_error("Cannot open file for reading halfedge_col_map");
	}

	boost::tie(halfedge_col_map, created) = mesh.add_property_map<halfedge_descriptor, unsigned>("h:color", 0);

	while (!in.eof())
	{
		size_t key;
		unsigned value;
		in.read(reinterpret_cast<char*>(&key), sizeof(key));
		in.read(reinterpret_cast<char*>(&value), sizeof(value));
		if (in)
		{ // Проверка на успешное чтение
			halfedge_descriptor hed(key);
			halfedge_col_map[hed] = value;
		}
	}
	in.close();

	// load diffr_edge_graph_edges
	std::cout << "load diffr_edge_graph_edges \n";
	in = std::ifstream(filename + "_diffr_edge_graph_edges", std::ios::binary);
	if (!in)
	{
		throw std::runtime_error("Cannot open file for reading diffr_edge_graph_edges");
	}

	while (!in.eof())
	{
		size_t _vd_first;
		size_t _vd_last;
		size_t _hed_first;
		size_t _hed_last;
		unsigned col;
		double wedge;

		in.read(reinterpret_cast<char*>(&_vd_first), sizeof(_vd_first));
		in.read(reinterpret_cast<char*>(&_vd_last), sizeof(_vd_last));
		in.read(reinterpret_cast<char*>(&_hed_first), sizeof(_hed_first));
		in.read(reinterpret_cast<char*>(&_hed_last), sizeof(_hed_last));
		in.read(reinterpret_cast<char*>(&col), sizeof(col));
		in.read(reinterpret_cast<char*>(&wedge), sizeof(wedge));

		vertex_descriptor vd_first(_vd_first);
		vertex_descriptor vd_last(_vd_last);
		halfedge_descriptor hed_first(_hed_first);
		halfedge_descriptor hed_last(_hed_last);

		BigDiffractionEdge edge1(mesh, vd_first, vd_last,
			hed_first, hed_last, col);
		edge1.wedge = wedge;

		std::vector<Vector> v_arr(3, Vector());

		for (size_t i = 0; i < v_arr.size(); i++)
		{
			double x;
			double y;
			double z;

			in.read(reinterpret_cast<char*>(&x), sizeof(x));
			in.read(reinterpret_cast<char*>(&y), sizeof(y));
			in.read(reinterpret_cast<char*>(&z), sizeof(z));

			v_arr[i] = Vector(x, y, z);
		}

		if (in)
		{ // Проверка на успешное чтение
			edge1.flag_direction = Vector(v_arr[0]);
			edge1.normal1 = Vector(v_arr[1]);
			edge1.normal2 = Vector(v_arr[2]);
			diffr_edge_graph.edges.push_back(edge1);
		}
	}
	in.close();
	total_num_edges = diffr_edge_graph.edges.size();
	return true;
}
bool SoundScene::loadDiffractionEdgeGraph(std::string filename)
{
	std::ifstream in;

	// load diffr_edge_graph
	std::cout << "load diffr_edge_graph \n";
	in = std::ifstream(filename + "_diffr_edge_graph", std::ios::binary);
	if (!in)
	{
		throw std::runtime_error("Cannot open file for reading diffr_edge_graph");
	}

	if (!in.eof())
	{
		for (size_t e1 = 0; e1 < total_num_edges; e1++)
		{

			BigDiffractionEdge& edge1 = diffr_edge_graph.edges[e1];
			unsigned list_offset;
			unsigned num_neighbors;

			in.read(reinterpret_cast<char*>(&list_offset), sizeof(list_offset));
			in.read(reinterpret_cast<char*>(&num_neighbors), sizeof(num_neighbors));


			if (in)
			{ // Проверка на успешное чтение
				edge1.list_offset = list_offset;
				edge1.num_neighbors = num_neighbors;
			}
		}
		while (!in.eof())
		{
			unsigned col;

			in.read(reinterpret_cast<char*>(&col), sizeof(col));

			if (in)
			{ // Проверка на успешное чтение
				diffr_edge_graph.diffr_edge_neighbors.push_back(col);
			}
		}
	}

	in.close();

	return true;
}

void SoundScene::clearSavedIRs(std::string filename)
{
	
	if (std::remove((filename + "_diffr_IRs").c_str()) != 0 || std::remove((filename + "_phong_IRs").c_str()) != 0)
	{
		perror("Ошибка при удалении файлов");
	}
	else
	{
		printf("Файлы успешно удалены\n");
	}
}

//////////////////////////////////
// Анализ результатов
//////////////////////////////////

std::vector<double> SoundScene::calculateAverageEnergyPerStep(const std::vector<std::pair<double, double>>& received_phong_IRs_archived, double dt)
{
	std::vector<std::pair<double, double>> sorted_phong_IRs = received_phong_IRs_archived;
	std::sort(sorted_phong_IRs.begin(), sorted_phong_IRs.end());

	std::vector<double> average_energy_per_step;
	double current_time = sorted_phong_IRs.front().first / speed_of_sound;
	double sum_energy = 0.0;
	int count = 0;

	for (const auto& pair : sorted_phong_IRs)
	{
		if (pair.first / speed_of_sound - current_time < dt)
		{
			sum_energy += pair.second;
			count++;
		}
		else
		{
			if (count > 0)
			{
				average_energy_per_step.push_back(sum_energy / count);
			}
			else
			{
				average_energy_per_step.push_back(0.0); // если нет данных в интервале, считаем среднее как 0
			}
			sum_energy = pair.second;
			count = 1;
			current_time = pair.first / speed_of_sound;
		}
	}

	if (count > 0)
	{
		average_energy_per_step.push_back(sum_energy / count);
	}
	else
	{
		average_energy_per_step.push_back(0.0); // если нет данных в последнем интервале, считаем среднее как 0
	}

	return average_energy_per_step;
}
