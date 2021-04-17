#ifndef EDGE_H
#define EDGE_H

#include <limits.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

class Vertex;
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
	Edge* getNext() const { assert(next != NULL); return next; }
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
		calculateLengthAndIndex();
	}
	void clearOpposite() {
		if (opposite == NULL) return;
		assert(opposite->opposite == this);
		opposite->opposite = NULL;
		opposite = NULL;
		//calculateLengthAndIndex();
	}
	void setNext(Edge* e) {
		assert(next == NULL);
		assert(e != NULL);
		assert(triangle == e->triangle);
		next = e;
	}
	void setCrease(float c) { crease = c; }

	//fields
	float getLength() {
		if (length == NULL) {
			this->calculateLengthAndIndex();
		}
		return length;  
	}
	int getIndexA() const { return indexA; }
	int getIndexB() const { return indexB; }

	//define operator so we can store edges in sorted data structures
	bool operator<(const Edge& e) { return length < e.length; }

	Matrix getError() const { return error; }
	void setError(Matrix e) { this->error = e; }

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
	Matrix error;
};

#endif
