#ifndef DIFFRACTIONEDGEGRAPH_H
#define DIFFRACTIONEDGEGRAPH_H

#include "utils.h"
#include "BigDiffractionEdge.h"

class DiffractionEdgeGraph
{
public:
	std::vector<BigDiffractionEdge> edges;

	std::vector<unsigned> diffr_edge_neighbors;

	DiffractionEdgeGraph()
	{
		edges = std::vector<BigDiffractionEdge>();
		diffr_edge_neighbors = std::vector<unsigned>();
	}

	const BigDiffractionEdge* getEdgeNeighbor(unsigned edgeNeighborIndex) const;
};

#endif