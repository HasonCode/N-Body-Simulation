#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include <GL/glew.h>
#include <SDL2/SDL.h>

char* file_read(const char* file_name){
    ifstream reader(file_name, ios::binary);
    if (!reader){
        cerr<< "Bail, Bail, Bail, no reader for " << file_name << endl;
        return NULL;
    }

    stringstream buffer;
    buffer << reader.rdbuf();
    reader.close();
    string to_out = buffer.str();
    if (to_out.size()==0){
        cerr<<"FILE IS EMPTY BAIL BAIL" << endl;
    }

    char* output = (char*)malloc(to_out.size()+1);
    strcpy(output,to_out.c_str());
    output[to_out.size()]='\0';
    return output; 
}


char* file_read2(const char* filename) {
	SDL_RWops *rw = SDL_RWFromFile(filename, "rb");
	if (rw == NULL) return NULL;
	
	Sint64 res_size = SDL_RWsize(rw);
	char* res = (char*)malloc(res_size + 1);

	Sint64 nb_read_total = 0, nb_read = 1;
	char* buf = res;
	while (nb_read_total < res_size && nb_read != 0) {
		nb_read = SDL_RWread(rw, buf, 1, (res_size - nb_read_total));
		nb_read_total += nb_read;
		buf += nb_read;
	}
	SDL_RWclose(rw);
	if (nb_read_total != res_size) {
		free(res);
		return NULL;
	}
	
	res[nb_read_total] = '\0';
	return res;
}

void print_log(GLuint something){
    GLint length = 0;
    if (glIsShader(something)){
        glGetShaderiv(something, GL_INFO_LOG_LENGTH, &length);
    }
    else if (glIsProgram(something)){
        glGetProgramiv(something, GL_INFO_LOG_LENGTH, &length);
    }
    else{
        cerr << "NOT A SHADER OR PROGRAM ERROR" << endl;
    }
    char* logger = (char*)malloc(length);
    if (glIsShader(something)){
        glGetShaderInfoLog(something, length,NULL, logger);
    }
    else if (glIsProgram(something)){
        glGetProgramInfoLog(something, length,NULL, logger);
    }
    cerr << logger;
    free(logger);
}

GLuint create_shader(const char* filename, GLenum type){
    const GLchar* source = file_read(filename);
    if (source==NULL){
        cout << "Error opening "<<filename<<endl;
    }
    GLuint shad = glCreateShader(type);
    glShaderSource(shad,1,&source,NULL);
    free((void*)source);
    glCompileShader(shad);
    GLint compile_ok = GL_FALSE;
    glGetShaderiv(shad,GL_COMPILE_STATUS,&compile_ok);
    if (compile_ok==GL_FALSE){
        cerr << filename << ":";
        print_log(shad);
        glDeleteShader(shad);
        return 0;
    }
    return shad;
}