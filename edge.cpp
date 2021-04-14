#ifndef _EDGE_H_
#define _EDGE_H_

#include "vertex.h"
#include "edge.h"

Edge::Edge(Vertex* v, Triangle* t): length(NULL), indexA(NULL), indexB(NULL) {
	vertex = v;
	triangle = t;
	next = NULL;
	opposite = NULL;
	crease = 0;
}

Edge::~Edge() {
	if (opposite != NULL)
		opposite->opposite = NULL;
}

void Edge::extract_func(Edge* e, int& a, int& b, int& c) {
	a = e->getVertex()->getIndex();
	b = e->getNext()->getNext()->getVertex()->getIndex();
	c = 0;
}

void Edge::calculateLengthAndIndex() {
	if (this->opposite == NULL) {
		throw "Cant calculate length/index because opposite is NULL";
	}
	else {
		this->length = sqrt(pow((this->vertex->x() - this->opposite->vertex->x()), 2) +
			pow((this->vertex->y() - this->opposite->vertex->y()), 2) +
			pow((this->vertex->z() - this->opposite->vertex->z()), 2));
		
		this->indexA = this->vertex->getIndex();
		this->indexB = this->opposite->vertex->getIndex();
		
		//std::cout << "calculated length: " << length << " and index: " << Index << std::endl;
	}
}

#endif
