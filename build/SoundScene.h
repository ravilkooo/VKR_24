#ifndef SOUNDSCENE_H
#define SOUNDSCENE_H

#include "DiffractionPropagation.h"
#include "PointSoundSource.h"
#include "PointSoundListener.h"
#include "SoundMaterial.h"
#include "SoundDiffractionPathPoint.h"
#include "SoundPhongPathPoint.h"
#include <cstdio>

// для расчёта нормалей к полигону мэша
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

// для нахождения пересечений с мэшем
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

// для нахождения границ мэша
#include <CGAL/Polygon_mesh_processing/bbox.h>

// для нахождения пересечений с мэшем
typedef Kernel::Ray_3 Ray;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

#define CGAL_DRAW_MESH false
#define CGAL_DRAW_FINAL_PATH false
#define DRAW_FINAL_PHONG_PATH false
#define DRAW_ALL_CORRECT_PATHES false
#define OUTPUT_MIS false


#define MIS_ENABLED true
#define MIS_EURISTIC 2

class SoundScene
{
public:
	// сетка сцены
	Mesh mesh;
	bool is_2d = false;

	// для быстрого поиска пересечения с объетками сцены
	Tree face_tree;

	SoundScene();

	bool readMeshFromFile(std::string filename);
	bool saveMeshFromFile(std::string filename);

	//////////////////////////////////
	// Всякие вспомогательные константы
	//////////////////////////////////

	unsigned num_frequencies = 3;
	std::vector<double> frequencies;
	double lowest_frequency = -1;
	double highest_frequency = -1;

	double min_diffraction_edge_angle = 0.5;
	double diffraction_threshold = std::abs(std::cos(degreesToRadians(mod(min_diffraction_edge_angle, 90.))));

	// просто я взял из головы
	double big_diffraction_edge_cos_angle_threshold = 0.95;
	double diffaction_wedge_variation_threshold = 0.95;
	double diffaction_cos_face_normal_variation_threshold = 0.95;

	double edge_resolution = 0.5;
	double edge_offset = 0.001;
	double min_rays_per_edge = 1;
	double max_rays_per_edge = 50;

	double ray_direction_threshold = 0.001;

	double ray_offset = 0.0001;
	// максимально допустимая длительность полученного импульса
	double max_IR_length = 0.1;

	double speed_of_sound = 343.99;

	size_t max_diffraction_order = 4;
	size_t max_ray_number = 1000;
	size_t max_diffr_ray_number = 10000;

	double diffraction_epsilon = 0.001;
	// shadow region tolerance
	double epsilon_h = edge_offset;
	// face region tolerance
	double epsilon_f = 0.01;

	//////////////////////////////////
	// Материал
	//////////////////////////////////

	std::vector<SoundMaterial> materials_list;

	// разметка материалов поверхности
	 Mesh::Property_map<face_descriptor, unsigned> faces_material_map;

	bool mapMaterials();

	//////////////////////////////////
	// Дифр рёбра
	//////////////////////////////////

	// разметка, какие ребра могут дифрагировать
	Mesh::Property_map<halfedge_descriptor, bool> diffr_map;
	// количество юольших дифр рёбр
	unsigned total_num_edges;
	// Граф видимости дифр рёбер
	DiffractionEdgeGraph diffr_edge_graph;
	// окраска смежных дифрагирующих рёбер в общий цвет
	Mesh::Property_map<halfedge_descriptor, unsigned> halfedge_col_map;
	// помечаем дифр рёбра
	bool markDiffractionEdges();
	// объединяем вместе в одно большие (окраска)
	bool combineDiffractionEdges();
	// построение графа видимости
	void buildDiffractionEdgeGraph();

	//////////////////////////////////
	// Источник/приёмник
	//////////////////////////////////
	
	// источник
	PointSoundSource sound_source;
	// приёмник
	PointSoundListener sound_listener;

	Vector listener_to_source_direction;

	// задать источник
	void setSoundSource(PointSoundSource new_sound_source);
	// задать приёмник
	void setSoundListener(PointSoundListener new_sound_listener);
	// получить направление от listener до source
	Vector listenerToSourceDirection();

	double pdfSound(const SoundPhongPathPoint& curr, const SoundPhongPathPoint& next);
	double pdfListener(const SoundPhongPathPoint& curr, const SoundPhongPathPoint& next);

	//////////////////////////////////
	// Распространение звука
	//////////////////////////////////

	// затухание от дистанции
	double getDistanceAttenuation(double distance);
	// делаем расчёты, распространяем звук от Source до Listener
	void propagateSound();

	//////////////////////////////////
	// Прямые пути
	//////////////////////////////////

	bool addDirectPath();

	//////////////////////////////////
	// Пути модели фонга (дифузия + зерк отраж)
	//////////////////////////////////

	// поиск путей модели фонга
	bool addSpecularPath(double frequency);
	// сохраняет путь (начало пути + точки в которых произошло отражение)
	// (в последней точке произошло отражение, а не затухание луча);
	// возвращает количество отражений
	int randomWalk(Ray first_ray, std::vector<SoundPhongPathPoint>& path,
		int max_depth, double pdf, double beta, double frequency);
	// подсчитать вес MIS (Multy Importance Sampling)
	double calcMISWeight(std::vector<SoundPhongPathPoint> source_path, std::vector<SoundPhongPathPoint> listener_path, int s, int t);

	//////////////////////////////////
	// Дифр пути
	//////////////////////////////////

	// поиск дифр пути
	bool addDiffractionPath(Vector ray_direction, Ray_intersection intersection);
	// рекурсивный поиск дифр пути
	bool recursiveDiffraction(ImpulseResponse& curr_IR, const BigDiffractionEdge& diffr_edge, unsigned depth);
	// находит точку на первой оси наиболее близкую ко второй оси
	double computePointOfClosestApproach(Point p1, Vector v1, Point p2, Vector v2);
	// проверка есть ли путь начин с ребра данного цвета
	bool alreadyContainsPath(unsigned diffr_edge_col);

	//////////////////////////////////
	// Флаги
	//////////////////////////////////
	
	// размер флага
	double d_flag = 1.;
	// диффр флаг или нет
	Mesh::Property_map<face_descriptor, unsigned> faces_flag_map;
	// добавить в сетки геометрию флагов
	bool addDiffractionFlags();

	class DiffractionQuery;

	//////////////////////////////////
	// Вспомогателььные методы
	//////////////////////////////////
	void setFrequencies(unsigned new_num_frequencies);
	void getSceneBoundingBox(Point& min_point, Point& max_point);
	bool isVisible(const Point& p1, const Point& p2);
	bool isBackSideVisible(const Point& point);

	//////////////////////////////////
	// Дополнительно для графики
	//////////////////////////////////

	void addArrow(Point start, Point end);
	std::vector <vertex_descriptor> v_buffer;
	std::vector <halfedge_descriptor> e_buffer;
	void clearArrowBuffers();
	void removeLastArrow();

	//////////////////////////////////
	// Сохранение и загрузка
	//////////////////////////////////
	std::string readInputFilename();
	std::string createOutputFilename(std::string input_filename);
	void saveOutputFilename(const std::string& data);
	void getInputParams(std::string& scene_name, std::string& answer,
		double& sx, double& sy, double& sz,
		double& dx, double& dy,
		double& x_low, double& x_high,
		double& y_low, double& y_high,
		unsigned& list_plane, double& plane_val,
		unsigned& phong_flag, bool& is_2d,
		size_t& max_diffr_order);
	void saveFrequencies(std::string filename);

	bool saveDiffractionEdges(std::string filename);
	bool saveDiffractionEdgeGraph(std::string filename);

	bool loadDiffractionEdges(std::string filename);
	bool loadDiffractionEdgeGraph(std::string filename);

	void clearSavedIRs(std::string filename);

	//////////////////////////////////
	// Анализ результатов
	//////////////////////////////////
	
	std::vector<double> calculateAverageEnergyPerStep(const std::vector<std::pair<double, double>>& received_phong_IRs_archived, double dt);
};

struct FlagSkip
{
	// скипать дифр флаги или нет
	bool skip_flags = false;
	// дополнительное поля для скипаемой поверхности (не флаг), если надо
	face_descriptor fd;
	const SoundScene& scene;
	FlagSkip(const SoundScene& scene)
		: scene(scene), fd(Mesh::null_face())
	{
		skip_flags = true;
	}
	FlagSkip(const SoundScene& scene, const face_descriptor fd)
		: scene(scene), fd(fd)
	{
	}
	FlagSkip(const SoundScene& scene, const face_descriptor fd, bool skip_flags)
		: scene(scene), fd(fd), skip_flags(skip_flags)
	{
	}
	bool operator()(const face_descriptor& t) const
	{
		// проверяем, является пересечёенная пов-сть флагом
		if (scene.faces_flag_map[t] > 0)
		{
			return skip_flags;
		}
		// проверяем (есть ли не флаг) нужно ли скипать эту пов-сть
		if (fd == Mesh::null_face())
		{
			// скипаемой пов-сти нет
			return false;
		}
		return (t == fd);
	}
};

#endif