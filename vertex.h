#ifndef _VERTEX_H
#define _VERTEX_H

#include <stdio.h>
#include <assert.h>
#include "vectors.h"
#include "matrix.h"

class Vertex; //todo is this line needed?

class Vertex {
public:
	Vertex(int i, const Vec3f& pos) : position(pos) {
		index = i;
		newFileIndex = NULL;

		Q = new Matrix(); //init Q

		//init Q
		//for (int i = 0; i < 4; i++) {
		//	for (int j = 0; j < 4; j++) {
		//		Q[i][j] = 0;
		//	}
		//}
	}
	virtual ~Vertex() { }

	// =========
	// ACCESSORS
	int getIndex() const { return index; }
	double x() const { return position.x(); }
	double y() const { return position.y(); }
	double z() const { return position.z(); }
	const Vec3f& get() const { return position; }

	// =========
	// MODIFIERS
	void set(Vec3f v) { position = v; }
	void set(double x, double y, double z) { position.Set(x, y, z); }

	//void setIndex(int i) { index = i; }
	void setNewFileIndex(int i) { newFileIndex = i; }
	int getNewFileIndex() const { return newFileIndex; }

	Matrix* getQ() const { return Q; };
	void setQ(int row, int col, float val) { Q->Set(row, col, val); };
	void addToQ(int row, int col, float val) { Q->Add(row, col, val); };
	//float* getQ() const { return (float*) Q; };
	void setQ(Matrix* m) { this->Q = m; }

	Vec3f getPosition() const { return position; };

private:
	// don't use these constructors
	Vertex() = delete;
	Vertex& operator=(const Vertex&) = delete;
	Vertex(const Vertex&) = delete;

	// ==============
	// REPRESENTATION
	Vec3f position;

	// this is the index from the original .obj file.
	// technically not part of the half-edge data structure
	int index;

	int newFileIndex;

	//float Q[4][4];
	Matrix* Q;

	// NOTE: the vertices don't know anything about adjacency.  In some
	// versions of this data structure they have a pointer to one of
	// their incoming edges.  However, this data is very complicated to
	// maintain during mesh manipulation.
};

#endif

