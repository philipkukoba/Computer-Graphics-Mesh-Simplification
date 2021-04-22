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
	//v_ = new Matrix();
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
#endif
