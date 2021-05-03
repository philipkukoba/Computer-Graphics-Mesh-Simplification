#include <stdio.h>
#include <string>

#include "argparser.h"
#include "glCanvas.h"
#include <GL/freeglut_std.h>  // VS might suggest here to install the package even though you have it installed -> reinstall it


int main(int argc, char* argv[]) {
	//srand48(0);
	srand(time(NULL));

	ArgParser* args = new ArgParser(argc, argv);
	Mesh* mesh = new Mesh();

	string input_file;

	/* To quickly select the input mesh for the demo (or development). */
	//input_file = "C:/Users/p_kuk/Desktop/UNI/Computer Graphics/project/data/bunny_1k.obj"; // Philip
	//input_file = "../../data/bunny_1k.obj"; // Laurens

	/* General input method*/
	if (args->input_file == NULL)
	{
		cout << "Enter file location" << endl;
		std::getline(cin, input_file);
	}
	else {
		input_file = args->input_file;
	}

	mesh->Load(input_file.c_str());

	GLCanvas glcanvas;
	glcanvas.initialize(args, mesh);

	// well it never returns from the GLCanvas loop...
	delete args;
	return 0;
}
