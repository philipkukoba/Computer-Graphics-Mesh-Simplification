#include <stdio.h>
#include <string>

#include "argparser.h"
#include "glCanvas.h"
#include <GL/freeglut_std.h>

// =========================================
// =========================================

int main(int argc, char *argv[]) {
  //srand48(0);
  ArgParser *args = new ArgParser(argc, argv);

  Mesh *mesh = new Mesh();
  string input_file;
  if (args->input_file == NULL)
  {
	  cout << "Enter file location" << endl;
	  std::getline(cin, input_file);
  }
  else
	  input_file = args->input_file;
  mesh->Load(input_file.c_str());

  GLCanvas glcanvas;
  glcanvas.initialize(args,mesh); 

  // well it never returns from the GLCanvas loop...
  delete args;
  return 0;

}

// =========================================
// =========================================
