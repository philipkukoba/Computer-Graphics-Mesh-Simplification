#define _USE_MATH_DEFINES
#include "glCanvas.h"

// Included files for OpenGL Rendering
#include <windows.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <corecrt_math_defines.h>

// static variables of GLCanvas class

// State of the mouse cursor
int GLCanvas::mouseButton;
int GLCanvas::mouseX;
int GLCanvas::mouseY;

int GLCanvas::display_list_index;
ArgParser* GLCanvas::args;
Camera* GLCanvas::camera;
Mesh* GLCanvas::mesh;

int GLCanvas::width;
int GLCanvas::height;

long long GLCanvas::lastClickTime = 0;

void GLCanvas::InitLight() {
	// Set the last component of the position to 0 to indicate
	// a directional light source

	GLfloat position[4] = { 30,30,100, 1 };
	GLfloat diffuse[4] = { 0.75,0.75,0.75,1 };
	GLfloat specular[4] = { 0,0,0,1 };
	GLfloat ambient[4] = { 0.2, 0.2, 0.2, 1.0 };

	GLfloat zero[4] = { 0,0,0,0 };
	glLightfv(GL_LIGHT1, GL_POSITION, position);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT1, GL_AMBIENT, zero);
	glEnable(GL_LIGHT1);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glEnable(GL_COLOR_MATERIAL);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

	GLfloat spec_mat[4] = { 1,1,1,1 };
	float glexponent = 30;
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &glexponent);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec_mat);

	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	float back_color[] = { 0.0,0.0,1.0,1 };
	glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, back_color);
	glEnable(GL_LIGHT1);
}


void GLCanvas::display(void)
{
	// Clear the display buffer, set it to the background color
	glClearColor(1, 1, 1, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Set the camera parameters
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	InitLight(); // light will be a headlamp!
	camera->glPlaceCamera();

	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);

	glCallList(display_list_index);
	HandleGLError();

	// Swap the back buffer with the front buffer to display
	// the scene
	glutSwapBuffers();
}

// ========================================================
// Callback function for window resize
// ========================================================

void GLCanvas::reshape(int w, int h) {
	width = w;
	height = h;

	// Set the OpenGL viewport to fill the entire window
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);

	// Set the camera parameters to reflect the changes
	camera->glInit(w, h);

	args->width = w;
	args->height = h;
}

// ========================================================
// Callback function for mouse click or release
// ========================================================

void GLCanvas::mouse(int button, int state, int x, int y) {
	// Save the current state of the mouse.  This will be
	// used by the 'motion' function
	mouseButton = button;
	mouseX = x;
	mouseY = y;


	long long currentClickTime = chrono::time_point_cast<chrono::milliseconds>(chrono::system_clock::now()).time_since_epoch().count();
	if (mouseButton == GLUT_LEFT_BUTTON && mesh->vertexSelectionMode != 0 && currentClickTime - lastClickTime > 250) // avoid double counting of clicks
		mesh->selectPoint(camera->center, camera->direction, camera->up, x, y, width, height);
	lastClickTime = currentClickTime;
}

// ========================================================
// Callback function for mouse drag
// ========================================================

void GLCanvas::motion(int x, int y) {
	// Left button = rotation
	// (rotate camera around the up and horizontal vectors)
	if (mouseButton == GLUT_LEFT_BUTTON) {
		camera->rotateCamera(0.005 * (mouseX - x), 0.005 * (mouseY - y));
		mouseX = x;
		mouseY = y;

		glutPostRedisplay();
	}
	// Middle button = translation
	// (move camera perpendicular to the direction vector)
	else if (mouseButton == GLUT_MIDDLE_BUTTON) {
		camera->truckCamera((mouseX - x) * 0.05, (y - mouseY) * 0.05);
		mouseX = x;
		mouseY = y;

		if (mesh->ProgressiveMeshing(camera->center))
			Render(); // (includes glutPostRedisplay)
		else
			glutPostRedisplay();
	}
	// Right button = zoom
	// (move camera along the direction vector)
	else if (mouseButton == GLUT_RIGHT_BUTTON) {
		camera->dollyCamera((x - mouseX) * 0.05);
		mouseX = x;
		mouseY = y;

		if (mesh->ProgressiveMeshing(camera->center))
			Render(); // (includes glutPostRedisplay)
		else
			glutPostRedisplay();
	}

	// Redraw the scene with the new camera parameters
}

// ========================================================
// Callback function for keyboard events
// ========================================================

void GLCanvas::keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 'w':  case 'W':
		args->wireframe = !args->wireframe;
		Render();
		break;
	case 'g': case 'G':
		args->gouraud = !args->gouraud;
		mesh->gouraud = !mesh->gouraud;
		if (mesh->gouraud)
			cout << "Using Gouraud shading." << endl;
		else
			cout << "Using flat shading." << endl;
		Render();
		break;
	case 's': case 'S':
		mesh->Save();
		printf("The mesh has been saved\n");
		break;
	case 'd': case 'D':
		mesh->Simplification((int)floor(0.9 * mesh->numTriangles()));
		Render();
		cout << "Simplified the mesh." << endl;
		break;
	case 'v': case 'V':
		mesh->vertexSelectionMode = true;
		printf("Left mouse click to select a vertex.\n");
		break;
	case 'r': case 'R':
		mesh->removeSelectedVertices();
		Render();
		break;
	case 'q':  case 'Q':
		exit(0); //quit
		break;
	case 'p': case 'P':
		mesh->progressiveMeshing = !mesh->progressiveMeshing;
		if (mesh->progressiveMeshing)
			cout << "Enabled progressive meshing." << endl;
		else
			cout << "Disabled progressive meshing." << endl;
		break;
	case 'm': case 'M':
		mesh->nextEdgeSelectionMode();
		break;
	case 'l': case 'L':
		mesh->toggleCollapseMidPoint();
		break;
	default:
		printf("UNKNOWN KEYBOARD INPUT  '%c'\n", key);
	}
}

// ========================================================
// Initialize all appropriate OpenGL variables, set
// callback functions, and start the main event loop.
// This function will not return but can be terminated
// by calling 'exit(0)'
// ========================================================

void GLCanvas::initialize(ArgParser* _args, Mesh* _mesh) {

	args = _args;
	mesh = _mesh;
	camera = new PerspectiveCamera(Vec3f(0, 0, 5), Vec3f(0, 0, -1), Vec3f(0, 1, 0), 20 * M_PI / 180.0);

	// Set global lighting parameters
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	// ? Some extra stuff to make sure glut initializes properly...
	char* myargv[1];
	int myargc = 1;
	myargv[0] = _strdup("Myappname");
	glutInit(&myargc, myargv);

	// Set window parameters
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
	glEnable(GL_DEPTH_TEST);
	glutInitWindowSize(args->width, args->height);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("OpenGL Viewer");

	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glEnable(GL_NORMALIZE);

	// Ambient light
	Vec3f ambColor = Vec3f(0.2, 0.2, 0.2);
	GLfloat ambArr[] = { ambColor.x(), ambColor.y(), ambColor.z(), 1.0 };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambArr);

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glCullFace(GL_BACK);
	glDisable(GL_CULL_FACE);

	display_list_index = glGenLists(1);

	// Initialize callback functions
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);

	Render();

	mesh->SetLodLevel0Distance(camera->center);

	// Enter the main rendering loop
	glutMainLoop();
}

void GLCanvas::Render() {
	glNewList(display_list_index, GL_COMPILE_AND_EXECUTE);
	// =========================================================
	// put your GL drawing calls inside the display list for efficiency
	mesh->Paint(args);
	// =========================================================
	glEndList();
	glutPostRedisplay();
}

void GLCanvas::DrawEllipsoid(unsigned int uiStacks, unsigned int uiSlices, float fA, float fB, float fC) {
	float tStep = (M_PI) / (float)uiSlices;
	float sStep = (M_PI) / (float)uiStacks;
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (float t = -M_PI / 2; t <= (M_PI / 2) + .0001; t += tStep)
	{
		glBegin(GL_TRIANGLE_STRIP);
		for (float s = -M_PI; s <= M_PI + .0001; s += sStep)
		{
			glVertex3f(fA * cos(t) * cos(s), fB * cos(t) * sin(s), fC * sin(t));
			glVertex3f(fA * cos(t + tStep) * cos(s), fB * cos(t + tStep) * sin(s), fC * sin(t + tStep));
		}
		glEnd();
	}
}

int HandleGLError() {
	GLenum error;
	int i = 0;
	while ((error = glGetError()) != GL_NO_ERROR) {
		printf("GL ERROR(%d):  %s\n", i, gluErrorString(error));
		i++;
	}
	if (i == 0) return 1;
	return 0;
}
