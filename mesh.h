#ifndef MESH_H
#define MESH_H

#include "vectors.h"
#include "array.h"
#include "bag.h"
#include "edge.h"
#include "boundingbox.h"
#include "argparser.h"
#include <queue>
#include "matrix.h"
#include <vector>
#include "my_priority_queue.h"

class Vertex;
class Triangle;
class VertexParent;

//mesh gebruikt een bag structuur voor zijn vertices edges en triangles etc

class Mesh {
public:
	// ========================
	// CONSTRUCTOR & DESTRUCTOR
	Mesh();
	virtual ~Mesh();
	void Load(const char* input_file, bool set_input_file);
	void Load(const char* input_file);

	// ========
	// VERTICES
	int numVertices() const { return vertices->Count(); }
	Vertex* addVertex(const Vec3f& pos);
	// this creates a relationship between 3 vertices (2 parents, 1 child)
	void setParentsChild(Vertex* p1, Vertex* p2, Vertex* child);
	// this accessor will find a child vertex (if it exists) when given
	// two parent vertices
	Vertex* getChildVertex(Vertex* p1, Vertex* p2) const;
	// look up vertex by index from original .obj file
	Vertex* getVertex(int i) const {
		assert(i >= 0 && i < numVertices());
		Vertex* v = (*vertices)[i];
		assert(v != NULL);
		return v;
	}

	// =====
	// EDGES
	int numEdges() const { return edges->Count(); }
	// this efficiently looks for an edge with the given vertices, using a hash table
	Edge* getEdge(Vertex* a, Vertex* b) const;

	// =========
	// TRIANGLES
	int numTriangles() const { return triangles->Count(); }
	void addTriangle(Vertex* a, Vertex* b, Vertex* c);
	void removeTriangle(Triangle* t);

	// ===============
	// OTHER ACCESSORS
	BoundingBox* getBoundingBox() const { return bbox; }

	// ===============
	// OTHER FUNCTIONS
	void Clear();

	void Paint(ArgParser* args);
	void LoopSubdivision();
	void CollapseEdge(Edge* e, float, float, float); // General collapse point, collapses all instances of AB (second level base method)
	void CollapseEdge(Edge* e); // Default value for the collapse point
	void CollapseEdge_MidPoint(Edge* e);
	void CollapseEdge_EndPoint(Edge* e);
	bool CollapseOneEdge_EndPoint(Edge* e); // Will be the base method where the others are built on top of, only collapses one instance of AB (does not delete AB or A)
	void CollapseRandomEdge();
	void CollapseShortestEdge();
	void CollapseQEM();
	void Simplification(int target_tri_count);
	
	void Save(std::string) const;
	void Save() const;

	bool LodLevelSaved(int, string&) const;
	bool ProgressiveMeshing(Vec3f); // Returns whether we have changed the mesh (hence need to redraw)

	// Quadric error metric helper functions
	void InitQuadricErrorMetric(Vertex* a, Vertex* b, Vertex* c);
	void computeContractionAndError(Edge* e, bool isMidPoint = true);

	void nextEdgeSelectionMode();
	void toggleCollapseMidPoint();

	int vertexSelectionMode = 0; // 0: not selecting, 1: select point 1, 2: select point 2
	void selectPoint(Vec3f cam_center, Vec3f cam_direction, Vec3f cam_up, int x, int y, int w, int h);
	void removeSelectedVertices();

	bool gouraud = false;
	bool progressiveMeshing = false;

	void SetLodLevel0Distance(Vec3f);

private:
	// ==============
	// REPRESENTATION
	Array<Vertex*>* vertices;
	Bag<Edge*>* edges;
	Bag<Triangle*>* triangles;
	BoundingBox* bbox;
	Bag<VertexParent*>* vertex_parents;

	int maxVertexIndex;

	const char* input_file;

	class EdgeComparer {
	public:
		int operator() (Edge* const e1, Edge* const e2) {
			return e1->getLength() > e2->getLength();
		}
	};

	class EdgeComparerQEM {
	public:
		int operator() (Edge* const e1, Edge* const e2) {
			return e1->getError() > e2->getError();
		}
	};

	my_priority_queue <Edge*, vector<Edge*>, EdgeComparer>* edgesShortestFirst;
	my_priority_queue <Edge*, vector<Edge*>, EdgeComparerQEM>* edgesQEM;

	std::vector < std::vector<int>> connectedVertices;

	Vertex* selectedPoint1 = NULL;
	Vertex* selectedPoint2 = NULL;

	int lodLevel = 0;
	float lodLevel0Distance = -1.; // not set yet

	int edgeSelectionMode = 1; // 0: random, 1: shortest, 2: QEM (can be overwritten by vertexSelectionMode)
	bool edgeCollapseToMidPoint = true; // (not applicable for QEM) 
};

#endif




