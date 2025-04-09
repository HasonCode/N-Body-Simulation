#version 120
attribute vec3 coord3d;
attribute vec3 v_color;
attribute vec3 normvec;
varying vec3 f_color;
varying vec3 normal_vector;
varying vec3 f_pos;
uniform mat4 mvp;
uniform mat4 model;
void main(void){
    gl_Position =  mvp * vec4(coord3d,1.0);
    f_color = v_color;
    f_pos = vec3(model * vec4(coord3d,1.0));
    normal_vector = normvec;
}