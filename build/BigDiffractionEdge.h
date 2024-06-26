#ifndef BIGDIFFRACTIONEDGE_H
#define BIGDIFFRACTIONEDGE_H

#include "utils.h"

class BigDiffractionEdge
{
public:
	const Mesh& m;
	// vd_first                                  vd_last
	//   O___>hed_first>__O_..._O____<hed_last<_____O
	vertex_descriptor vd_first = Mesh::null_vertex();
	vertex_descriptor vd_last = Mesh::null_vertex();
	halfedge_descriptor hed_first = Mesh::null_halfedge();
	halfedge_descriptor hed_last = Mesh::null_halfedge();

	unsigned col = 0;

	Vector flag_direction;
	double wedge;
	Vector normal1;
	Vector normal2;

	unsigned list_offset = 0;
	unsigned num_neighbors = 0;

	face_descriptor flag_face_descriptor;

	BigDiffractionEdge(const Mesh& m) : m(m)
	{
	}

	BigDiffractionEdge(const Mesh& m, vertex_descriptor vd_first, vertex_descriptor vd_last,
		halfedge_descriptor hed_first, halfedge_descriptor hed_last, unsigned col) :
		m(m), vd_first(vd_first), vd_last(vd_last),
		hed_first(hed_first), hed_last(hed_last),
		col(col)
	{
	}

	// ����� �����
	double getLength() const;

	// �������� ����������� �� ������ ������� �� ������
	Vector getAxis() const;

	// ��������� ������ �������
	Point getStart() const;

	// �������� ��� ����� � ����
	bool testOrientation(const Point point, double offset) const;

	// �������� ��� ���� (���� �� ��������) � ����
	bool testEdgeOrientation(const BigDiffractionEdge& edge2, double offset);

	// �������� ��������� ����� �� ������� ����� 1
	Point getFreeVertex1() const;

	// �������� ��������� ����� �� ������� ����� 2
	Point getFreeVertex2() const;

	// �������� ��������� �� ������� ����� 1
	Plane getPlane1() const;

	// �������� ��������� �� ������� ����� 1
	Plane getPlane2() const;
};

#endif