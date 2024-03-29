#ifndef EDGE_H
#define EDGE_H

#include <limits.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "vertex.h"


class Triangle;
class Matrix;

// half-edge data structure

class Edge {
public:
	Edge(Vertex* v, Triangle* t);
	~Edge();

	// here's the hash function to use for edges so they
	// can be efficiently accessed within the Bag data structure
	static void extract_func(Edge* e, int& a, int& b, int& c);

	// =========
	// ACCESSORS
	Vertex* getVertex() const { assert(vertex != NULL); return vertex; }
	Edge* getNext() const { 
		assert(next != NULL);
		return next; 
	}
	bool hasNext() const { return (next != NULL); }
	Triangle* getTriangle() const { assert(triangle != NULL); return triangle; }
	Edge* getOpposite() const {
		// warning!  the opposite edge might be NULL!
		return opposite;
	}

	float getCrease() const { return crease; }

	Vertex* operator[](int i) const {
		if (i == 0) return getVertex();
		if (i == 1) return getNext()->getNext()->getVertex();
		assert(0);
	}

	// =========
	// MODIFIERS
	void setOpposite(Edge* e) {
		assert(opposite == NULL);
		assert(e != NULL);
		assert(e->opposite == NULL);
		opposite = e;
		e->opposite = this;
	}
	void clearOpposite() {
		if (opposite == NULL) return;
		assert(opposite->opposite == this);
		opposite->opposite = NULL;
		opposite = NULL;
	}
	void setNext(Edge* e) {
		assert(next == NULL);
		assert(e != NULL);
		assert(triangle == e->triangle);
		next = e;
		length = NULL;
	}
	void setCrease(float c) { crease = c; }

	//fields
	float getLength() {
		if (length == NULL)
			length = ((*this)[0]->get() - (*this)[1]->get()).Length();
		return length;
	}
	int getIndexA() const { 
		return (*this)[1]->getIndex(); 
	}
	int getIndexB() const { 
		return (*this)[0]->getIndex(); 
	}

	//unused
	////define operator so we can store edges in sorted data structures
	//bool operator<(Edge& e) { return getLength() < e.getLength(); }

	float getError() const { return error; }
	void setError(float e) { this->error = e; }
	
	/*Matrix* getV_() const { return v_; }
	void setV_(Matrix* v__) { this->v_ = v__; }*/

	Vec4f getV_() const { return v_; }
	void setV_(Vec4f v__) { this->v_ = v__; }

private:
	Edge(const Edge&) = delete;
	Edge& operator=(const Edge&) = delete;

	// REPRESENTATION
	// in the half edge data adjacency data structure, the edge stores everything!
	Vertex* vertex;
	Triangle* triangle;
	Edge* opposite;
	Edge* next;
	float crease;

	float length;

	int indexA;
	int indexB;
	
	//calculate length (for shortest edge collapse)
	void calculateLengthAndIndex();

	// for quadric error metric
	float error;
	//Matrix* v_; //matrix representation of a vertex
	Vec4f v_;
};

#endif
