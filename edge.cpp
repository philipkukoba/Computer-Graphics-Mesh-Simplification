#ifndef _EDGE_H_
#define _EDGE_H_

#include "vertex.h"
#include "edge.h"

Edge::Edge(Vertex* v, Triangle* t): length(NULL) {
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

void Edge::calculateLength() {
	if (this->opposite == NULL) throw "Tried calculating length of edge without an opposite";
	length = sqrt(pow((this->vertex->x() - this->opposite->vertex->x()), 2) +
		pow((this->vertex->y() - this->opposite->vertex->z()), 2) +
		pow((this->vertex->y() - this->opposite->vertex->z()), 2));
}

#endif
