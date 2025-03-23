#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;

#include <GL/glew.h>
#include <SDL2/SDL.h>
char* file_read(const char* file_name){
    ifstream reader(file_name);
    string ret_val,line;
    getline(reader,line);
    ret_val = "";
    ret_val+=line+"\n";
    while (getline(reader,line)){
        ret_val = ret_val + line;
    }
    reader.close();
    char* output = (char*)malloc(ret_val.size()+1);
    strcpy(output,ret_val.c_str());
    return output; 
}

int main(){
    const char *fs_source =
		//"#version 100\n"  // OpenGL ES 2.0
		"#version 120\n"  // OpenGL 2.1
		"void main(void) {        "
		"  gl_FragColor[0] = 0.0; "
		"  gl_FragColor[1] = 0.0; "
		"  gl_FragColor[2] = 1.0; "
		"}";
    printf(file_read("shader.frag"));
    printf("\n");
    printf(fs_source);
    return 1;
}