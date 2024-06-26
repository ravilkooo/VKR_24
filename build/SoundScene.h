#ifndef SOUNDSCENE_H
#define SOUNDSCENE_H

#include "DiffractionPropagation.h"
#include "PointSoundSource.h"
#include "PointSoundListener.h"
#include "SoundMaterial.h"
#include "SoundDiffractionPathPoint.h"
#include "SoundPhongPathPoint.h"
#include <cstdio>

// ��� ������� �������� � �������� ����
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

// ��� ���������� ����������� � �����
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

// ��� ���������� ������ ����
#include <CGAL/Polygon_mesh_processing/bbox.h>

// ��� ���������� ����������� � �����
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
	// ����� �����
	Mesh mesh;
	bool is_2d = false;

	// ��� �������� ������ ����������� � ��������� �����
	Tree face_tree;

	SoundScene();

	bool readMeshFromFile(std::string filename);
	bool saveMeshFromFile(std::string filename);

	//////////////////////////////////
	// ������ ��������������� ���������
	//////////////////////////////////

	unsigned num_frequencies = 3;
	std::vector<double> frequencies;
	double lowest_frequency = -1;
	double highest_frequency = -1;

	double min_diffraction_edge_angle = 0.5;
	double diffraction_threshold = std::abs(std::cos(degreesToRadians(mod(min_diffraction_edge_angle, 90.))));

	// ������ � ���� �� ������
	double big_diffraction_edge_cos_angle_threshold = 0.95;
	double diffaction_wedge_variation_threshold = 0.95;
	double diffaction_cos_face_normal_variation_threshold = 0.95;

	double edge_resolution = 0.5;
	double edge_offset = 0.001;
	double min_rays_per_edge = 1;
	double max_rays_per_edge = 50;

	double ray_direction_threshold = 0.001;

	double ray_offset = 0.0001;
	// ����������� ���������� ������������ ����������� ��������
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
	// ��������
	//////////////////////////////////

	std::vector<SoundMaterial> materials_list;

	// �������� ���������� �����������
	 Mesh::Property_map<face_descriptor, unsigned> faces_material_map;

	bool mapMaterials();

	//////////////////////////////////
	// ���� ����
	//////////////////////////////////

	// ��������, ����� ����� ����� �������������
	Mesh::Property_map<halfedge_descriptor, bool> diffr_map;
	// ���������� ������� ���� ���
	unsigned total_num_edges;
	// ���� ��������� ���� ����
	DiffractionEdgeGraph diffr_edge_graph;
	// ������� ������� ������������� ���� � ����� ����
	Mesh::Property_map<halfedge_descriptor, unsigned> halfedge_col_map;
	// �������� ���� ����
	bool markDiffractionEdges();
	// ���������� ������ � ���� ������� (�������)
	bool combineDiffractionEdges();
	// ���������� ����� ���������
	void buildDiffractionEdgeGraph();

	//////////////////////////////////
	// ��������/�������
	//////////////////////////////////
	
	// ��������
	PointSoundSource sound_source;
	// �������
	PointSoundListener sound_listener;

	Vector listener_to_source_direction;

	// ������ ��������
	void setSoundSource(PointSoundSource new_sound_source);
	// ������ �������
	void setSoundListener(PointSoundListener new_sound_listener);
	// �������� ����������� �� listener �� source
	Vector listenerToSourceDirection();

	double pdfSound(const SoundPhongPathPoint& curr, const SoundPhongPathPoint& next);
	double pdfListener(const SoundPhongPathPoint& curr, const SoundPhongPathPoint& next);

	//////////////////////////////////
	// ��������������� �����
	//////////////////////////////////

	// ��������� �� ���������
	double getDistanceAttenuation(double distance);
	// ������ �������, �������������� ���� �� Source �� Listener
	void propagateSound();

	//////////////////////////////////
	// ������ ����
	//////////////////////////////////

	bool addDirectPath();

	//////////////////////////////////
	// ���� ������ ����� (������� + ���� �����)
	//////////////////////////////////

	// ����� ����� ������ �����
	bool addSpecularPath(double frequency);
	// ��������� ���� (������ ���� + ����� � ������� ��������� ���������)
	// (� ��������� ����� ��������� ���������, � �� ��������� ����);
	// ���������� ���������� ���������
	int randomWalk(Ray first_ray, std::vector<SoundPhongPathPoint>& path,
		int max_depth, double pdf, double beta, double frequency);
	// ���������� ��� MIS (Multy Importance Sampling)
	double calcMISWeight(std::vector<SoundPhongPathPoint> source_path, std::vector<SoundPhongPathPoint> listener_path, int s, int t);

	//////////////////////////////////
	// ���� ����
	//////////////////////////////////

	// ����� ���� ����
	bool addDiffractionPath(Vector ray_direction, Ray_intersection intersection);
	// ����������� ����� ���� ����
	bool recursiveDiffraction(ImpulseResponse& curr_IR, const BigDiffractionEdge& diffr_edge, unsigned depth);
	// ������� ����� �� ������ ��� �������� ������� �� ������ ���
	double computePointOfClosestApproach(Point p1, Vector v1, Point p2, Vector v2);
	// �������� ���� �� ���� ����� � ����� ������� �����
	bool alreadyContainsPath(unsigned diffr_edge_col);

	//////////////////////////////////
	// �����
	//////////////////////////////////
	
	// ������ �����
	double d_flag = 1.;
	// ����� ���� ��� ���
	Mesh::Property_map<face_descriptor, unsigned> faces_flag_map;
	// �������� � ����� ��������� ������
	bool addDiffractionFlags();

	class DiffractionQuery;

	//////////////////////////////////
	// ���������������� ������
	//////////////////////////////////
	void setFrequencies(unsigned new_num_frequencies);
	void getSceneBoundingBox(Point& min_point, Point& max_point);
	bool isVisible(const Point& p1, const Point& p2);
	bool isBackSideVisible(const Point& point);

	//////////////////////////////////
	// ������������� ��� �������
	//////////////////////////////////

	void addArrow(Point start, Point end);
	std::vector <vertex_descriptor> v_buffer;
	std::vector <halfedge_descriptor> e_buffer;
	void clearArrowBuffers();
	void removeLastArrow();

	//////////////////////////////////
	// ���������� � ��������
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
	// ������ �����������
	//////////////////////////////////
	
	std::vector<double> calculateAverageEnergyPerStep(const std::vector<std::pair<double, double>>& received_phong_IRs_archived, double dt);
};

struct FlagSkip
{
	// ������� ���� ����� ��� ���
	bool skip_flags = false;
	// �������������� ���� ��� ��������� ����������� (�� ����), ���� ����
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
		// ���������, �������� ������������� ���-��� ������
		if (scene.faces_flag_map[t] > 0)
		{
			return skip_flags;
		}
		// ��������� (���� �� �� ����) ����� �� ������� ��� ���-���
		if (fd == Mesh::null_face())
		{
			// ��������� ���-��� ���
			return false;
		}
		return (t == fd);
	}
};

#endif