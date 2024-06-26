#include "DiffractionEdgeGraph.h"

const BigDiffractionEdge* DiffractionEdgeGraph::getEdgeNeighbor(unsigned edgeNeighborIndex) const
{
	assert(edgeNeighborIndex < diffr_edge_neighbors.size());
	unsigned edgeIndex = diffr_edge_neighbors[edgeNeighborIndex];

	assert(edgeIndex < edges.size());

	return &edges[edgeIndex];
}
