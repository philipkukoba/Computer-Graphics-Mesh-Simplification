#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <assert.h>
#include <Windows.h>
#include <GL/gl.h>

#include <iostream>
#include <fstream>
#include <regex>
#include <string>

#include "mesh.h"
#include "edge.h"
#include "vertex.h"
#include "triangle.h"
#include "vertex_parent.h"
#include "glCanvas.h"

#define INITIAL_VERTEX 10000
#define INITIAL_EDGE 10000
#define INITIAL_TRIANGLE 10000

// =======================================================================
// CONSTRUCTORS & DESTRUCTORS
// =======================================================================

Mesh::Mesh() : input_file(NULL) {
	vertices = new Array<Vertex*>(INITIAL_VERTEX);
	edges = new Bag<Edge*>(INITIAL_EDGE, Edge::extract_func);
	triangles = new Bag<Triangle*>(INITIAL_TRIANGLE, Triangle::extract_func);
	vertex_parents = new Bag<VertexParent*>(INITIAL_VERTEX, VertexParent::extract_func);
	edgesShortestFirst = new priority_queue<Edge*, std::vector<Edge*>, EdgeComparer>();
	bbox = NULL;
}

Mesh::~Mesh() {
	delete vertices;
	vertices = NULL;
	delete edges;
	edges = NULL;
	delete triangles;
	triangles = NULL;
	delete bbox;
	bbox = NULL;
}

// =======================================================================
// MODIFIERS:   ADD & REMOVE
// =======================================================================

Vertex* Mesh::addVertex(const Vec3f& position) {
	int index = numVertices();
	Vertex* v = new Vertex(index, position);
	vertices->Add(v);
	if (bbox == NULL)
		bbox = new BoundingBox(position, position);
	else
		bbox->Extend(position);
	return v;
}

void Mesh::addTriangle(Vertex* a, Vertex* b, Vertex* c) {

	// create the triangle
	Triangle* t = new Triangle();

	// create the edges
	Edge* ea = new Edge(a, t);
	Edge* eb = new Edge(b, t);
	Edge* ec = new Edge(c, t);

	// point the triangle to one of its edges
	t->setEdge(ea);

	// connect the edges to each other
	ea->setNext(ec);
	eb->setNext(ea);
	ec->setNext(eb);

	// add them to the master list
	edges->Add(ea);
	edges->Add(eb);
	edges->Add(ec);

	// connect up with opposite edges (if they exist)
	Edge* ea_op = getEdge((*ea)[1], (*ea)[0]);
	Edge* eb_op = getEdge((*eb)[1], (*eb)[0]);
	Edge* ec_op = getEdge((*ec)[1], (*ec)[0]);
	if (ea_op != NULL) {
		ea_op->setOpposite(ea);
		edgesShortestFirst->push(ea_op);
	}
	if (eb_op != NULL) {
		eb_op->setOpposite(eb);
		edgesShortestFirst->push(eb_op);
	}
	if (ec_op != NULL) {
		ec_op->setOpposite(ec);
		edgesShortestFirst->push(ec_op);
	}

	// add the triangle to the master list
	triangles->Add(t);
}

void Mesh::removeTriangle(Triangle* t) {

	Edge* ea = t->getEdge();
	Edge* eb = ea->getNext();
	Edge* ec = eb->getNext();
	assert(ec->getNext() == ea);

	// remove elements from master lists
	edges->Remove(ea);
	edges->Remove(eb);
	edges->Remove(ec);
	triangles->Remove(t);

	// clean up memory
	delete ea;
	delete eb;
	delete ec;
	delete t;
}

Edge* Mesh::getEdge(Vertex* a, Vertex* b) const {
	assert(edges != NULL);
	//std::cout << "used edges->get(a,b)" << std::endl;
	return edges->Get(a->getIndex(), b->getIndex());
}

Vertex* Mesh::getChildVertex(Vertex* p1, Vertex* p2) const {
	VertexParent* vp = vertex_parents->GetReorder(p1->getIndex(), p2->getIndex());
	if (vp == NULL) return NULL;
	return vp->get();
}

void Mesh::setParentsChild(Vertex* p1, Vertex* p2, Vertex* child) {
	vertex_parents->Add(new VertexParent(p1, p2, child));
}

// =======================================================================
// the load function parses very simple .obj files
// the basic format has been extended to allow the specification 
// of crease weights on the edges.
// =======================================================================

void Mesh::Load(const char* input_file) {

	FILE* objfile = fopen(input_file, "r");
	if (objfile == NULL) {
		printf("ERROR! CANNOT OPEN '%s'\n", input_file);
		return;
	}
	this->input_file = input_file;

	char line[200];
	char token[100];
	char atoken[100];
	char btoken[100];
	char ctoken[100];
	char dtoken[100];
	char etoken[100];
	float x, y, z;
	int a, b, c, d, e;

	int index = 0;
	int vert_count = 0;
	int vert_index = 1;

	while (fgets(line, 200, objfile)) {

		if (line[strlen(line) - 2] == '\\') {
			fgets(token, 100, objfile);
			int tmp = strlen(line) - 2;
			strncpy(&line[tmp], token, 100);
		}
		int token_count = sscanf(line, "%s\n", token);
		if (token_count == -1) continue;
		a = b = c = d = e = -1;
		if (!strcmp(token, "usemtl") ||
			!strcmp(token, "g")) {
			vert_index = 1; //vert_count + 1;
			index++;
		}
		else if (!strcmp(token, "v")) {
			vert_count++;
			sscanf(line, "%s %f %f %f\n", token, &x, &y, &z);
			addVertex(Vec3f(x, y, z));
		}
		else if (!strcmp(token, "f")) {
			int num = sscanf(line, "%s %s %s %s %s %s\n", token,
				atoken, btoken, ctoken, dtoken, etoken);
			sscanf(atoken, "%d", &a);
			sscanf(btoken, "%d", &b);
			sscanf(ctoken, "%d", &c);
			if (num > 4) sscanf(dtoken, "%d", &d);
			if (num > 5) sscanf(etoken, "%d", &e);
			a -= vert_index;
			b -= vert_index;
			c -= vert_index;
			if (d >= 0) d -= vert_index;
			if (e >= 0) e -= vert_index;
			assert(a >= 0 && a < numVertices());
			assert(b >= 0 && b < numVertices());
			assert(c >= 0 && c < numVertices());

			addTriangle(getVertex(a), getVertex(b), getVertex(c));
			if (d > -1) { assert(d < numVertices()); addTriangle(getVertex(a), getVertex(c), getVertex(d)); }
			if (e > -1) { assert(e < numVertices()); addTriangle(getVertex(a), getVertex(d), getVertex(e)); }
		}
		else if (!strcmp(token, "e")) {
			int num = sscanf(line, "%s %s %s %s\n", token, atoken, btoken, ctoken);
			assert(num == 4);
			sscanf(atoken, "%d", &a);
			sscanf(btoken, "%d", &b);
			if (!strcmp(ctoken, "inf")) x = 1000000; // this is close to infinity...
			else sscanf(ctoken, "%f", &x);
			Vertex* va = getVertex(a);
			Vertex* vb = getVertex(b);
			Edge* ab = getEdge(va, vb);
			Edge* ba = getEdge(vb, va);
			assert(ab != NULL);
			assert(ba != NULL);
			ab->setCrease(x);
			ba->setCrease(x);
		}
		else if (!strcmp(token, "vt")) {
		}
		else if (!strcmp(token, "vn")) {
		}
		else if (token[0] == '#') {
		}
		else {
			printf("LINE: '%s'", line);
		}
	}
}

// =======================================================================
// PAINT
// =======================================================================

Vec3f ComputeNormal(const Vec3f& p1, const Vec3f& p2, const Vec3f& p3) {
	Vec3f v12 = p2;
	v12 -= p1;
	Vec3f v23 = p3;
	v23 -= p2;
	Vec3f normal;
	Vec3f::Cross3(normal, v12, v23);
	normal.Normalize();
	return normal;
}

void InsertNormal(const Vec3f& p1, const Vec3f& p2, const Vec3f& p3) {
	Vec3f normal = ComputeNormal(p1, p2, p3);
	glNormal3f(normal.x(), normal.y(), normal.z());
}

void Mesh::Paint(ArgParser* args) {

	// scale it so it fits in the window
	Vec3f center; bbox->getCenter(center);
	float s = 1 / bbox->maxDim();
	glScalef(s, s, s);
	glTranslatef(-center.x(), -center.y(), -center.z());

	// this offset prevents "z-fighting" between the edges and faces
	// the edges will always win.
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.1, 4.0);

	// draw the triangles
	glColor3f(1, 1, 1);
	Iterator<Triangle*>* iter = triangles->StartIteration();
	glBegin(GL_TRIANGLES);
	while (Triangle* t = iter->GetNext()) {
		Vec3f a = (*t)[0]->get();
		Vec3f b = (*t)[1]->get();
		Vec3f c = (*t)[2]->get();
		InsertNormal(a, b, c);
		glVertex3f(a.x(), a.y(), a.z());
		glVertex3f(b.x(), b.y(), b.z());
		glVertex3f(c.x(), c.y(), c.z());
	}
	triangles->EndIteration(iter);
	glEnd();

	glDisable(GL_POLYGON_OFFSET_FILL);

	if (args->wireframe) {
		glDisable(GL_LIGHTING);

		// draw all the interior, non-crease edges
		glLineWidth(1);
		glColor3f(0, 0, 0);
		glBegin(GL_LINES);
		Iterator<Edge*>* iter = edges->StartIteration();
		while (Edge* e = iter->GetNext()) {
			if (e->getOpposite() == NULL || e->getCrease() > 0) continue;
			Vec3f a = (*e)[0]->get();
			Vec3f b = (*e)[1]->get();
			glVertex3f(a.x(), a.y(), a.z());
			glVertex3f(b.x(), b.y(), b.z());
		}
		edges->EndIteration(iter);
		glEnd();

		// draw all the interior, crease edges
		glLineWidth(3);
		glColor3f(1, 1, 0);
		glBegin(GL_LINES);
		iter = edges->StartIteration();
		while (Edge* e = iter->GetNext()) {
			if (e->getOpposite() == NULL || e->getCrease() == 0) continue;
			Vec3f a = (*e)[0]->get();
			Vec3f b = (*e)[1]->get();
			glVertex3f(a.x(), a.y(), a.z());
			glVertex3f(b.x(), b.y(), b.z());
		}
		edges->EndIteration(iter);
		glEnd();

		// draw all the boundary edges
		glLineWidth(3);
		glColor3f(1, 0, 0);
		glBegin(GL_LINES);
		iter = edges->StartIteration();
		while (Edge* e = iter->GetNext()) {
			if (e->getOpposite() != NULL) continue;
			assert(e->getCrease() == 0);
			Vec3f a = (*e)[0]->get();
			Vec3f b = (*e)[1]->get();
			glVertex3f(a.x(), a.y(), a.z());
			glVertex3f(b.x(), b.y(), b.z());
		}
		edges->EndIteration(iter);
		glEnd();

		glEnable(GL_LIGHTING);
	}

	HandleGLError();
}

// =================================================================
// SUBDIVISION
// =================================================================

void Mesh::LoopSubdivision() {
	printf("Subdivide the mesh!\n");
}

// =================================================================
// SIMPLIFICATION
// =================================================================

void Mesh::CollapseEdge_EndPoint(Edge* e, bool deleteE) {
	// Collapses the edge AB to the single point B.

	assert(e->getNext() != NULL);
	assert(e->getOpposite() != NULL);
	assert(e->getTriangle() != NULL);
	// The same should hold for all edges we encounter below, so we should really just start from a well-initialized mesh.

	Edge* AB = e;
	Edge* BC = AB->getNext();
	Edge* CA = BC->getNext();

	Vertex* A = CA->getVertex();
	Vertex* B = AB->getVertex();
	Vertex* C = BC->getVertex();
	Triangle* ABC = AB->getTriangle();

	Edge* BA = AB->getOpposite();
	Edge* AD = BA->getNext();
	Edge* DB = AD->getNext();
	Vertex* D = AD->getVertex();
	Triangle* ADB = AD->getTriangle();

	// Remove A, edges and triangles, but delay deleting them
	vertices->Remove(A);
	edges->Remove(AB);
	edges->Remove(BC); // The new BC will come from AC later.
	edges->Remove(CA);
	edges->Remove(AD);
	edges->Remove(DB); // The new DB will come from DA in the next step.
	edges->Remove(BA);
	triangles->Remove(ABC);
	triangles->Remove(ADB);


	// We'll now cycle through all other triangles containing A.
	bool cycleCompleted = false;
	Edge* eCycle = AD->getOpposite()->getNext();
	Edge* lastBP = NULL;
	while (eCycle->getVertex() != B)
	{
		Edge* AP = eCycle;
		Edge* PQ = AP->getNext();
		Edge* QA = PQ->getNext();
		Vertex* P = AP->getVertex(); Vertex* Q = PQ->getVertex();
		if (A->getIndex() == P->getIndex() || A->getIndex() == Q->getIndex())
		{
			bool b = A == B;
			int debugme2 = 0; // TODO: can happen(!): AA is an edge  (sometimes only (but not always!) indices are the same, but points are different)
		}
		Triangle* APQ = AP->getTriangle();

		if (P == B || P == Q || Q == B)
			int debugme3 = 0; // TODO: P = B can happen!

		Triangle* newBPQ = new Triangle();
		Edge* newBP = new Edge(P, newBPQ); // Will replace AP
		Edge* newPQ = new Edge(Q, newBPQ); // Will replace the old PQ
		Edge* newQB = new Edge(B, newBPQ); // Will replace QA

		newBPQ->setEdge(newBP);

		// Set up newBP
		newBP->setCrease(AP->getCrease());
		newBP->setNext(newPQ);
		// The opposite of newBP (PB) comes from the next triangle's newQB (where the next Q is our P), and will hence be set in the next step of the cycle.

		// Set up newPQ
		newPQ->setCrease(PQ->getCrease());
		newPQ->setNext(newQB);
		Edge* QP = PQ->getOpposite(); // Lies in a triangle which does not contain A, so no worries
		QP->clearOpposite();
		newPQ->setOpposite(QP); // Also sets QP's opposite to newPQ. So no need to also use QP->setOpposite(newPQ); In fact, this is not allowed as this method asserts that newPQ's opposite is NULL.

		// Set up newQB
		newQB->setCrease(QA->getCrease());
		newQB->setNext(newBP);
		// The opposite of newQB (BQ) comes from the previous triangle's newBP (where the previous P is our Q).
		if (Q == D)
		{
			// Special case if the previous triangle is deleted (hence there is no lastBP): when we are now in APD.
			Edge* BD = DB->getOpposite();
			BD->clearOpposite();
			newQB->setOpposite(BD);
		}
		else
		{
			// No need to use clearOpposite() as lastBP's opposite has not yet been set.
			newQB->setOpposite(lastBP); // Also sets the previous step's newBP's opposite.
		}

		// Remove old edges and triangles, add new ones
		edges->Remove(AP);
		edges->Remove(PQ);
		edges->Remove(QA);
		triangles->Remove(APQ);

		edges->Add(newBP);
		edges->Add(newPQ);
		edges->Add(newQB);

		/*edgesShortestFirst->push(newBP);
		edgesShortestFirst->push(newPQ);
		edgesShortestFirst->push(newQB);*/


		triangles->Add(newBPQ);

		lastBP = newBP; // Remember for the next triangle
		Edge* nextECycle = eCycle->getOpposite()->getNext();

		delete AP; AP = NULL; // Put all of the deletes after all of the Removes, as Remove looks at opposites, which might be set to NULL otherwise.
		delete PQ; PQ = NULL;
		delete QA; QA = NULL;
		delete APQ; APQ = NULL;
		// Includes deleting eCycle (AP).

		eCycle = nextECycle;


		//debugging code
		/*int missingOpposites = 0;
		Iterator<Edge*>* iter = edges->StartIteration();
		while (Edge* e = iter->GetNext()) {
			if (e->getOpposite() == NULL)
				missingOpposites++;
		}
		edges->EndIteration(iter);*/
		/*std::cout
			<< "Missing opposites (in method): " << missingOpposites << std::endl;*/
	}


	// Since we're always delaying setting the opposite of newBP, this still needs to be set for the last triangle. In that case (the last) P == C.
	if (lastBP != NULL)
		// This should always be the case unless we actually didn't loop above. Then AD->getOpposite()->getNext() == AB and AD == AB->'getPrevious()'->getOpposite() == CA->getOpposite() == AC. 
		// Thus C = D and A lies in only one triangle (excluding flipped orientation).
	{
		Edge* CB = BC->getOpposite();
		CB->clearOpposite();
		lastBP->setOpposite(CB);
	}
	else
	{
		Edge* CB = BC->getOpposite();
		CB->clearOpposite();
		Edge* BD = DB->getOpposite();
		BD->clearOpposite();
		CB->setOpposite(BD);
	}
	// When A lies in only one triangle we can simply delete our edge BC (==BD), together with their opposite CB==DB, as they all lie in the same triangle (flipped pair). We do not need to replace them.
	// In all cases these will be deleted below.

	// Finally, delete the edges from the two fully deleted triangles containing AB or BA. (We were not able to do this earlier as we still needed BC and DB).
	delete BC; BC = NULL;
	delete CA; CA = NULL;
	delete ABC; ABC = NULL;
	delete AD; AD = NULL;
	delete DB; DB = NULL;
	delete BA; BA = NULL;
	delete ADB; ADB = NULL;
	if (deleteE)
	{
		delete A; A = NULL;
		delete AB; AB = NULL;
	}
}

void Mesh::CollapseEdge_EndPoint(Edge* e) {
	// When collapsing an edge AB to B is the final goal, we should delete A and e=AB.
	CollapseEdge_EndPoint(e, true);
}

void Mesh::CollapseEdge(Edge* e, float collapse_x, float collapse_y, float collapse_z)
{
	Vertex* A = (*e)[1];
	Vertex* B = (*e)[0];
	CollapseEdge_EndPoint(e, false);
	B->set(collapse_x, collapse_y, collapse_z);
	delete A;
	delete e;
}

void Mesh::CollapseEdge_MidPoint(Edge* e)
{
	Vertex* A = (*e)[1];
	Vertex* B = (*e)[0];
	CollapseEdge(e, (A->x() + B->x()) / 2., (A->y() + B->y()) / 2., (A->z() + B->z()) / 2.);
}

void Mesh::CollapseEdge(Edge* e)
{
	// Default
	return CollapseEdge_MidPoint(e);
}


void Mesh::CollapseRandomEdge() {
	CollapseEdge(edges->ChooseRandom());
}

void Mesh::CollapseShortestEdge() {

	//if the heap is empty it needs to be refreshed
	//if (edgesShortestFirst->empty()) {

	//}

	Edge* e = edgesShortestFirst->top();

	//controleer of de edge niet al verwijderd is in de bag
	//of als de edge geen opposite heeft
	while (e->getOpposite() == NULL || 
		e->getIndexA() == e->getIndexB() || 
		edges->Get(e->getIndexA(), e->getIndexB()) == NULL)
	{
		edgesShortestFirst->pop();
		e = edgesShortestFirst->top();
	}

	CollapseEdge_EndPoint(edgesShortestFirst->top());
	edgesShortestFirst->pop();
}

void Mesh::Simplification(int target_tri_count) {
	//CollapseShortestEdge(); return;

	while (numTriangles() > target_tri_count)
	{
		//CollapseRandomEdge();
		CollapseShortestEdge();

		//debug code
		/*Iterator<Edge*>* iter = edges->StartIteration();
		int missingNexts = 0;
		int missingOpposites = 0;
		int oppositesNotInEdges = 0;
		int oppositesDifferentPoints = 0;
		int oppositesInSameTriangle = 0;
		int oppositesInFlippedTriangle = 0;
		while (Edge* e = iter->GetNext()) {
			if (e->getNext() == NULL)
				missingNexts++;
			if (e->getOpposite() == NULL)
				missingOpposites++;
			else
			{
				if (edges->Member(e->getOpposite()) == 0)
					oppositesNotInEdges++;
				if ((*e->getOpposite())[1] != (*e)[0])
					oppositesDifferentPoints++;

				Triangle* t1 = e->getTriangle();
				Triangle* t2 = e->getOpposite()->getTriangle();
				Vec3f n1 = ComputeNormal((*t1)[0]->get(), (*t1)[1]->get(), (*t1)[2]->get());
				Vec3f n2 = ComputeNormal((*t2)[0]->get(), (*t2)[1]->get(), (*t2)[2]->get());
				if (n1 == n2)
					oppositesInSameTriangle++;
				n2.Negate();
				if (n1 == n2)
					oppositesInFlippedTriangle++;
			}
		}
		edges->EndIteration(iter);*/
		/*std::cout << "Number of vertices: " << vertices->Count() << std::endl
			<< "Number of edges: " << edges->Count() << std::endl
			<< "Number of triangles: " << triangles->Count() << std::endl
			<< "Missing nexts: " << missingNexts << std::endl
			<< "Missing opposites: " << missingOpposites << std::endl
			<< "Opposites of edges not contained in edges: " << oppositesNotInEdges << std::endl
			<< "Opposites having different vertices: " << oppositesDifferentPoints << std::endl
			<< "Opposites in same triangle (orientation): " << oppositesInSameTriangle << std::endl
			<< "Opposites in flipped triangle (orientation): " << oppositesInFlippedTriangle << std::endl << std::endl;*/
	}
}

void Mesh::Save() const
{
	int triangle_count = triangles->Count();
	std::string str(input_file); // c string to c++ string
	std::string filename = std::regex_replace(str, std::regex(".obj"), std::string(""))
		+ "_simplified_to_"
		+ std::to_string(triangle_count)
		+ ".obj";

	ofstream myfile(filename);


	if (myfile.is_open())
	{
		//write all vertices
		
		int count = vertices->Count();
		for (int i = 0; i < count; i++) {

			//set new index (starts at 1)
			vertices->operator[](i)->setIndex(i + 1);

			myfile << "v "
				+ std::to_string(vertices->operator[](i)->x()) + ' '
				+ std::to_string(vertices->operator[](i)->y()) + ' '
				+ std::to_string(vertices->operator[](i)->z()) + '\n';
			//std::cout << "index: " << verticesSorted->operator[](i)->getIndex() << std::endl;
		}

		//write all faces (triangles)

		Iterator<Triangle*>* iterT = triangles->StartIteration();
		while (Triangle* t = iterT->GetNext()) {
			myfile << "f "
				+ std::to_string(t->operator[](0)->getIndex()) + ' '
				+ std::to_string(t->operator[](1)->getIndex()) + ' '
				+ std::to_string(t->operator[](2)->getIndex()) + '\n';
		}
		triangles->EndIteration(iterT);

		myfile.close();
	}
	else throw "Unable to save mesh.";
}

void Mesh::selectPoint(Vec3f cam_center, Vec3f cam_direction, Vec3f cam_up, int x, int y, int w, int h)
{
	if (vertexSelectionMode == 0)
		return;

	// Find the line defined by the camera center and the selected 2D point (x, y).
	// Step 1: Zero center in 2D: x - w/2. Similarly for y.
	// Step 2: Negate: -(x - w/2) as in 2D x points to the right, but in 3D x point to the left (in 3D y points up, z forward, hence x left). Same for y.
	// Step 3: Convert from pixels to world units. Let f be the focal length (distance to the camera plane; canonically this would be 1, but it does not matter). 
	// Assuming a (canonical) horizontal field of view of 90°, the point (0, h/2), now transformed to (w/2, 0) should correspond to (1, 0, 1) in 3D 
	//  (since (f, 0, f), (0, 0, 0) and (-f, 0, f) form an isosceles triangle with right angle at (0, 0, 0)). Thus converting from pixels to world units is done 
	//	by multiplying by 2f/w. We then end up at (-(x - w/2)*2f/w, -(y - h/2)*2f/w, f) as transformed version of (x, y).
	// Step 4: On this line we also have (by scaling) (-(x - w/2)*2/w, -(y - h/2)*2/w, 1).
	// Step 4': However, after testing and working out an example, it turns out we need to put the z-coordinate here to 12/5 = 2.4 for some reason.
	// Step 5: Transform to world coordinates using the camera matrix.

	Vec3f line_dir(-(x - w/2.)*2./w, -(y - h/2.)*2./w, 2.4); // Direction vector of line between camera center (0, 0, 0) and clicked (2D) point - camera coordinates
	
	Vec3f look = cam_direction;
	Vec3f up = cam_up - look * cam_up.Dot3(look); // up-direction of the camera, but orthogonal to the view direction
	up.Normalize();
	Vec3f hor;
	Vec3f::Cross3(hor, up, look);
	
	Matrix transf; // from the camera system to world coordinates (i.e. maps (0, 0, 0, 1) to (cam_center, 1), and (1, 0, 0, 0) to (cam_horizontal, 0) etc.)
	transf.Set(0, 0, hor.x()); transf.Set(0, 1, up.x()); transf.Set(0, 2, look.x()); transf.Set(0, 3, cam_center.x());
	transf.Set(1, 0, hor.y()); transf.Set(1, 1, up.y()); transf.Set(1, 2, look.y()); transf.Set(1, 3, cam_center.y());
	transf.Set(2, 0, hor.z()); transf.Set(2, 1, up.z()); transf.Set(2, 2, look.z()); transf.Set(2, 3, cam_center.z());
	transf.Set(3, 0, 0.);				  transf.Set(3, 1, 0.);			transf.Set(3, 2, 0.);				 transf.Set(3, 3, 1.);

	transf.TransformDirection(line_dir); // world coordinates
	line_dir.Normalize();

	// Naively iterate over all vertices, take closest (should be fast enough)
	float min_dist = FLT_MAX;
	Vertex* closestPoint = NULL;
	for (int i = 0; i < vertices->Count(); i++)
	{
		Vertex* vertex = (*vertices)[i];
		Vec3f q = vertex->get();
		float dist = ((q - cam_center) - line_dir * (q - cam_center).Dot3(line_dir)).Length(); // Remove v component from q - p, where p is any point on the line. What remains is the orthogonal component.
		if (dist < min_dist)
		{
			closestPoint = vertex;
			min_dist = dist;
		}
	}

	cout << "Selected the point at (" << closestPoint->x() << ", " << closestPoint->y() << ", " << closestPoint->z() << ") (index=" << closestPoint->getIndex() << "). ";
	if (vertexSelectionMode == 1)
	{
		selectedPoint1 = closestPoint;
		cout << "This is the first selected point. Select one more point." << endl;
		vertexSelectionMode = 2;
	}
	else
	{
		selectedPoint2 = closestPoint;
		cout << "This is the second selected point. Press 'r' to (try to) collapse the edge between them." << endl;
		vertexSelectionMode = 0;
	}
}

void Mesh::removeSelectedVertices()
{
	if (selectedPoint1 == NULL)
	{
		cout << "First select a point!" << endl;
		return;
	}
	if (selectedPoint2 == NULL)
	{
		cout << "First select a second point!" << endl;
		return;
	}
	Edge* e = getEdge(selectedPoint1, selectedPoint2); // If these do not form an edges, this will return NULL.
	if (e == NULL)
		cout << "The two selected vertices must form an edge!" << endl;
	else
		Mesh::CollapseEdge(e);
}

// =================================================================
