#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;

#include <GL/glew.h>
GLuint program;
GLint attribute_coord2d, attribute_v_color;
GLint attribute_coord3d;
GLint attribute_normvec;
GLuint vbo_triangle, vbo_triangle_colors;
GLint uniform_fade;
GLuint vbo_sphere, vbo_sphere_colors;
GLuint ibo_sphere_elements;
GLuint vbo_norm_vectors;
GLint uniform_mvp;
GLint uniform_lightpos;
GLint uniform_model;
#include <SDL2/SDL_image.h>
#include <SDL2/SDL.h>
#include "shader_utils.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
GLint uniform_m_transform;

int screen_width = 800, screen_height=600;


 
class particle_system{
    public:
        float grav_const = 6.674 * pow(10.0,-11.0);
        int layers = 50;
        int num_points = 50;
        float bounce_dist = 50.0;
        float bounce_return = 0.9;
        int num_vertices = (layers-2)*(num_points+1)+2;
        int num_triangles =(layers-3)*(2*(num_points+1))+2*(num_points+1);
        int size = 0;
        float rotation_speed = 0.00;
        float rotation_offset = 0; 
        float intensity = 10000;
        GLfloat* light_pos;
        float framerate;
        float* positions;
        float* velocities;
        float* forces;
        float* radii;
        float* masses;
        particle_system(float frames)
        {
            framerate=frames;
            light_pos = (GLfloat*)malloc(3*sizeof(GLfloat));
            light_pos[0] = 0;
            light_pos[1] = 0;
            light_pos[2] = 0;
            positions = (float*)malloc(3*sizeof(float));
            velocities = (float*)malloc(3*sizeof(float));
            forces = (float*)malloc(3*sizeof(float));
            radii = (float*)malloc(1*sizeof(float));
            masses = (float*)malloc(1*sizeof(float));

        }
        float distance_calculator(int ind1, int ind2){
            float distancex,distancey,distancez;
            distancex = pow(positions[ind1*3]-positions[ind2*3],2);
            distancey = pow(positions[ind1*3+1]-positions[ind2*3+1],2);
            distancez = pow(positions[ind1*3+2]-positions[ind2*3+2],2);
            float distance = sqrt(distancex+distancey+distancez);
            return distance;
        }
        void add_particle(float* position, float* velocity, float* force, float mass, float radius){
            float* temp_pos = (float*)malloc((size+1)*3*sizeof(float));
            float* temp_vel = (float*)malloc((size+1)*3*sizeof(float));
            float* temp_force = (float*)malloc((size+1)*3*sizeof(float));
            float* temp_mass = (float*)malloc((size+1)*sizeof(float));
            float* temp_radius = (float*)malloc((size+1)*sizeof(float));

            for (int i = 0; i<size;i++){
                temp_pos[3*i]=positions[3*i];
                temp_pos[3*i+1]=positions[3*i+1];
                temp_pos[3*i+2]=positions[3*i+2];
                temp_vel[3*i]=velocities[3*i];
                temp_vel[3*i+1]=velocities[3*i+1];
                temp_vel[3*i+2]=velocities[3*i+2];
                temp_force[3*i]=forces[3*i];
                temp_force[3*i+1]=forces[3*i+1];
                temp_force[3*i+2]=forces[3*i+2];
                temp_mass[i] = masses[i];
                temp_radius[i] = radii[i];
            }
            temp_pos[size*3]=position[0];
            temp_pos[size*3+1]=position[1];
            temp_pos[size*3+2]=position[2];

            temp_vel[size*3]=velocity[0];
            temp_vel[size*3+1]=velocity[1];
            temp_vel[size*3+2]=velocity[2];

            temp_force[size*3]=force[0];
            temp_force[size*3+1]=force[1];
            temp_force[size*3+2]=force[2];

            temp_mass[size] = mass;
            temp_radius[size] = radius;
            free(positions);
            free(velocities);
            free(forces);
            free(radii);
            free(masses);
            positions = temp_pos;
            velocities = temp_vel;
            forces = temp_force;
            radii = temp_radius;
            masses = temp_mass;
            size++;
        }
        void update_positions(){
            for (int i = 0; i<size;i++){
                positions[i*3]+=velocities[i*3]/framerate;
                positions[i*3+1]+=velocities[i*3+1]/framerate;
                positions[i*3+2]+=velocities[i*3+2]/framerate;
                if (positions[i*3]>bounce_dist){
                    positions[i*3]= 2*bounce_dist-positions[i*3];
                    velocities[i*3]= -velocities[i*3]*bounce_return;
                }
                else if (positions[i*3]<-1*bounce_dist){
                    positions[i*3]= -2*bounce_dist-positions[i*3]; 
                    velocities[i*3] = -velocities[i*3]*bounce_return;
                }
                if (positions[i*3+1]>bounce_dist){
                    positions[i*3+1]= 2*bounce_dist-positions[i*3+1];
                    velocities[i*3+1]= -velocities[i*3+1]*bounce_return;
                }
                else if (positions[i*3+1]<-1*bounce_dist){
                    positions[i*3+1]= -2*bounce_dist-positions[i*3+1];
                    velocities[i*3+1] = -velocities[i*3+1]*bounce_return;
                }
                if (positions[i*3+2]>bounce_dist){
                    positions[i*3+2]= 2*bounce_dist-positions[i*3+2];
                    velocities[i*3+2]= -velocities[i*3+2]*bounce_return;
                }
                else if (positions[i*3+2]<-1*bounce_dist){
                    positions[i*3+2]= -2*bounce_dist-positions[i*3+2];
                    velocities[i*3+2] = -velocities[i*3+2]*bounce_return;
                }
            }

        }
        void update_velocities(){
            for (int i = 0; i<size;i++){
                velocities[i*3]+=(forces[i*3]/framerate)/masses[i/3];
                velocities[i*3+1]+=(forces[i*3+1]/framerate)/masses[i/3];
                velocities[i*3+2]+=(forces[i*3+2]/framerate)/masses[i/3];
            }
        }
        float gravity_equation(float mass1,float mass2, float distance){
            float num = grav_const*mass1*mass2;
            float denom = pow(distance,2.0);
            return num/denom;
        }
        void clear_forces(){
            for (int i = 0; i<size; i++){
                forces[i*3]=0;
                forces[i*3+1]=0;
                forces[i*3+2]=0;
            }
        }
        void calc_force(int ind1, int ind2){
            float x1 = positions[ind1*3],x2 = positions[ind2*3];
            float y1 = positions[ind1*3+1],y2 = positions[ind2*3+1];
            float z1 = positions[ind1*3+2],z2 = positions[ind2*3+2];
            float x = x2-x1, y = y2-y1, z = z2-1;
            if (x<0){
                x*=-1;
            }
            if (y<0){
                y*=-1;
            }
            if (z<0){
                z*=-1;
            }
            float theta1 = asinf(y/distance_calculator(ind1,ind2));
            float theta2 = acosf((double)x/(double)distance_calculator(ind1,ind2));
            float force = gravity_equation(masses[ind1],masses[ind2],distance_calculator(ind1,ind2));
            float force_x = force*sin(theta1)*cos(theta2);
            float force_y = force*sin(theta1);
            float force_z = force*sin(theta1)*sin(theta2);
            if (x1>x2){
                forces[ind1*3]-=force_x;
                forces[ind2*3]+=force_x;
            }
            else{
                forces[ind1*3]+=force_x;
                forces[ind2*3]-=force_x;
            }

            if (y1>y2){
                forces[ind1*3+1]-=force_y;
                forces[ind2*3+1]+=force_y;
            }
            else{
                forces[ind1*3+1]+=force_y;
                forces[ind2*3+1]-=force_y;
            }

            if (z1>z2){
                forces[ind1*3+2]-=force_z;
                forces[ind2*3+2]+=force_z;
            }
            else{
                forces[ind1*3+2]+=force_z;
                forces[ind2*3+2]-=force_z;
            }
        }
        void update_forces(){
            clear_forces();
            for (int i = 0; i<size;i++){
                for (int j = i+1; j<size;j++){
                    calc_force(i,j);
                }
            }
        }
        float calc_speed(int ind1){
            float velx = velocities[ind1*3];
            float vely = velocities[ind1*3+1];
            float velz = velocities[ind1*3+2];
            return sqrt(pow(velx,2.0)+pow(vely,2.0)+pow(velz,2.0));
        }
        void merge_masses(int ind1, int ind2){
            float volume = 16.0/9.0 * M_PI*M_PI * radii[ind1]
            * radii[ind1] + radii[ind2]*radii[ind2];  
            float new_radius = pow(pow(radii[ind1]*10,3)+pow(radii[ind2]*10,3),1.0/3.0)/10;
            float new_mass = masses[ind1]+masses[ind2];
            float velx = (masses[ind1]*velocities[ind1*3]+ velocities[ind2*3]* masses[ind2])/new_mass;
            float vely = (masses[ind1]*velocities[ind1*3+1]+ velocities[ind2*3+1]* masses[ind2])/new_mass;
            float velz = (masses[ind1]*velocities[ind1*3+2]+ velocities[ind2*3+2]* masses[ind2])/new_mass;
            float forcex = forces[ind1*3]+forces[ind2*3];
            float forcey = forces[ind1*3+1]+forces[ind2*3+1];
            float forcez = forces[ind1*3+2]+forces[ind2*3+2];

            float posx,posy,posz;
            if (radii[ind1]>radii[ind2]){
                posx = positions[ind1*3];
                posy = positions[ind1*3+1];
                posz = positions[ind1*3+2];
            }
            else if (radii[ind1]<radii[ind2]){
                posx = positions[ind2*3];
                posy = positions[ind2*3+1];
                posz = positions[ind2*3+2];
            }
            else{
                if (masses[ind1]*calc_speed(ind1)>masses[ind2]*calc_speed(ind2)){
                    posx = positions[ind1*3];
                    posy = positions[ind1*3+1];
                    posz = positions[ind1*3+2];
                }
                else{
                    posx = positions[ind2*3];
                    posy = positions[ind2*3+1];
                    posz = positions[ind2*3+2];
                }
            }
            float* temp_pos = (float*)malloc((size-1)*3*sizeof(float));
            float* temp_vel = (float*)malloc((size-1)*3*sizeof(float));
            float* temp_force = (float*)malloc((size-1)*3*sizeof(float));
            float* temp_mass = (float*)malloc((size-1)*sizeof(float));
            float* temp_radius = (float*)malloc((size-1)*sizeof(float));
            int count = 0 ;
            for (int i = 0; i<size;i++){
                if (i!=ind1 && i!=ind2){
                    temp_pos[3*count]=positions[3*i];
                    temp_pos[3*count+1]=positions[3*i+1];
                    temp_pos[3*count+2]=positions[3*i+2];
                    temp_vel[3*count]=velocities[3*i];
                    temp_vel[3*count+1]=velocities[3*i+1];
                    temp_vel[3*count+2]=velocities[3*i+2];
                    temp_force[3*count]=forces[3*i];
                    temp_force[3*count+1]=forces[3*i+1];
                    temp_force[3*count+2]=forces[3*i+2];
                    temp_mass[count] = masses[i];
                    temp_radius[count] = radii[i];
                    count++;
                }
            }

            temp_pos[(size-1)*3-3] = (positions[ind1*3]+positions[ind2*3])/2.0;
            temp_pos[(size-1)*3-2] = (positions[ind1*3+1]+positions[ind2*3+1])/2.0;
            temp_pos[(size-1)*3-1] = (positions[ind1*3+2]+positions[ind2*3+2])/2.0;

            temp_vel[(size-1)*3-3] = velx;
            temp_vel[(size-1)*3-2] = vely;
            temp_vel[(size-1)*3-1] = velz;
            
            temp_force[(size-1)*3-3] = forcex;
            temp_force[(size-1)*3-2] = forcey;
            temp_force[(size-1)*3-1] = forcez;
            
            temp_radius[size-2] = new_radius;
            temp_mass[size-2] = new_mass; 
            size--;
            free(positions);
            free(velocities);
            free(forces);
            free(masses);
            free(radii);
            positions = temp_pos;
            velocities = temp_vel;
            forces = temp_force;
            masses = temp_mass;
            radii = temp_radius;
        }
    bool can_merge(int ind1, int ind2){
        float distance = distance_calculator(ind1,ind2);
        return distance<=(radii[ind1]) || distance<=(radii[ind2]);
    }
    bool update_merges(){
        bool flag = false;
        bool temp_flag = true;
        while (temp_flag){
            temp_flag = false;
            for (int i = 0; i < size; i++){
                if(temp_flag){
                    break;
                }
                for (int j = i+1; j < size; j++){
                    if (can_merge(i,j)){
                        temp_flag = true;
                        flag = true;
                        merge_masses(i,j);
                        i=0;
                        break;
                    }
                }
            }
    }
        return flag;
    }
    void updater(){
        if (!update_merges()){
            update_positions();
            update_velocities();
            update_forces();
        }
    }
    GLfloat* cross_product(GLfloat* vec1, GLfloat* vec2){
        GLfloat* ret_val = (GLfloat*)malloc(3*sizeof(GLfloat));
        ret_val[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
        ret_val[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
        ret_val[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
        return ret_val;
    }
    GLfloat* calc_normal_vector(int tnum, GLuint* triangles,GLfloat* vertices){
        GLfloat* vec1 = (GLfloat*)malloc(3*sizeof(GLfloat));
        GLfloat* vec2 = (GLfloat*)malloc(3*sizeof(GLfloat));
        vec1[0] = vertices[(triangles[tnum*3])*3]-vertices[(triangles[tnum*3+1])*3];
        vec1[1] = vertices[(triangles[tnum*3])*3+1]-vertices[(triangles[tnum*3+1])*3+1];
        vec1[2] = vertices[(triangles[tnum*3])*3+2]-vertices[(triangles[tnum*3+1])*3+2];

        vec2[0] = vertices[(triangles[tnum*3+1])*3]-vertices[(triangles[tnum*3+2])*3];
        vec2[1] = vertices[(triangles[tnum*3+1])*3+1]-vertices[(triangles[tnum*3+2])*3+1];
        vec2[2] = vertices[(triangles[tnum*3+1])*3+2]-vertices[(triangles[tnum*3+2])*3+2];

        GLfloat* ret_vec = cross_product(vec1,vec2);
        float mag = sqrt(ret_vec[0]*ret_vec[0]+ret_vec[1]*ret_vec[1]+ret_vec[2]*ret_vec[2]);
        free(vec1);
        free(vec2);
        ret_vec[0]/=mag;
        ret_vec[1]/=mag;
        ret_vec[2]/=mag;
        return ret_vec;
    }
    GLfloat* subtraction(GLfloat* vec1, GLfloat* vec2, int size){
        GLfloat* ret_vec = (GLfloat*)malloc(size*sizeof(GLfloat));
        for (int i = 0; i<size;i++){
            ret_vec[i] = vec1[i]-vec2[i];
        }
        return ret_vec;
    }
    GLfloat* get_norm_vectors(){
        GLfloat* vertices = generate_spheres();
        GLuint* triangles = get_all_indices();
        GLfloat* norm_vectors = (GLfloat*)malloc(sizeof(GLfloat)*size*num_vertices*3);
        float* count_array  = (float*)malloc(3*size*num_vertices*sizeof(float));
        for (int j = 0; j< size*num_vertices*3;j++){
            norm_vectors[j] = 0;
            count_array[j] = 0;
        }
        for (int i = 0; i < size; i++){
            for (int tri = 0; tri<num_triangles; tri++){
                GLfloat* norm_vector;
                if(true){
                    norm_vector = calc_normal_vector(i*num_triangles+tri,triangles,vertices);
                }
                for (int j = 0; j < 3; j++){
                    if (norm_vector[0]){
                        norm_vectors[triangles[(i*num_triangles+tri)*3+j]*3] += norm_vector[0];
                        count_array[triangles[(i*num_triangles+tri)*3+j]*3]++;
                    }
                    if (norm_vector[1]){
                        norm_vectors[triangles[(i*num_triangles+tri)*3+j]*3+1] += norm_vector[1];
                        count_array[triangles[(i*num_triangles+tri)*3+j]*3+1]++;
                    }
                    if (norm_vector[2]){
                        norm_vectors[triangles[(i*num_triangles+tri)*3+j]*3+2] += norm_vector[2];
                        count_array[triangles[(i*num_triangles+tri)*3+j]*3+2]++;
                    }
                }
                free(norm_vector);
            }
        }
        for (int k = 0; k<size*num_vertices; k++){
            norm_vectors[k]/=count_array[k];
        }
        for (int rep = 0; rep<size; rep++){
            norm_vectors[rep*num_vertices*3]=0.0;
            norm_vectors[rep*num_vertices*3+1]=1.0;
            norm_vectors[rep*num_vertices*3+2]=0.0;
            if (rep>0){
                norm_vectors[rep*num_vertices*3-3]=0.0;
                norm_vectors[rep*num_vertices*3-2]=-1.0;
                norm_vectors[rep*num_vertices*3-1]=0.0;
            }
        }
        free(count_array);
        free(triangles);
        free(vertices);
        return norm_vectors;
    }   
    GLfloat* generate_sphere(int ind){
        GLfloat* ret_arr = (GLfloat*)malloc(num_vertices*3*sizeof(GLfloat));
        int p = 0;
        float init_y = -radii[ind];
        float radius = radii[ind];
        for (int layer = 0; layer<layers; layer++){
            for (int point = 0; point<=num_points;point++){
                float rot = (2*M_PI)/num_points*point;
                float calc_val = radius*radius - abs(init_y)*abs(init_y);
                if (calc_val<0){
                    calc_val=0.0;
                }
                float x = sqrt(calc_val)*cos(rot);
                float z = sqrt(calc_val)*sin(rot);
                ret_arr[p*3] = (GLfloat)(x+positions[ind*3]);
                ret_arr[p*3+1] = (GLfloat)(init_y+positions[ind*3+1]);
                ret_arr[p*3+2] = (GLfloat)(z+positions[ind*3+2]);
                if (layer==0||layer==layers-1){
                    ret_arr[p*3]=positions[ind*3];
                    ret_arr[p*3+2]=positions[ind*3+2];
                    p++;
                    break;
                }
                else{
                    p++;
                }
            }
            init_y+=(radii[ind]*2)/(layers-1);
        } 
        return ret_arr;
    }
    GLfloat* generate_spheres(){
        GLfloat* ret_arr = (GLfloat*)malloc(size*num_vertices*3*sizeof(GLfloat));
        int p=0;
        for (int sphere = 0; sphere<size;sphere++){
            float init_y = -radii[sphere];
            float radius = radii[sphere];
            for (int layer = 0; layer<layers; layer++){
                for (int point = 0; point<=num_points;point++){
                    float rot = (2*M_PI)/num_points*point;
                    float x = sqrt(radius*radius - abs(init_y)*abs(init_y))*cos(rot+rotation_offset);
                    float z = sqrt(radius*radius - abs(init_y)*abs(init_y))*sin(rot+rotation_offset);
                    ret_arr[p*3] = (GLfloat)(x+positions[sphere*3]);
                    ret_arr[p*3+1] = (GLfloat)(init_y+positions[sphere*3+1]);
                    ret_arr[p*3+2] = (GLfloat)(z+positions[sphere*3+2]);
                    if (layer==0||layer==layers-1){
                        ret_arr[p*3]=positions[sphere*3];
                        ret_arr[p*3+2]=positions[sphere*3+2];
                        p++;
                        break;
                    }
                    else{
                        p++;
                    }
                }
                init_y+=(radii[sphere]*2)/(layers-1);
            } 
        }
        rotation_offset+=rotation_speed;
        return ret_arr;
    }
    GLushort* get_indexes(){
        GLushort* fun_arr = (GLushort*)malloc(3*num_triangles*sizeof(GLshort));
        int point = 0;
        for (int layer = 0; layer < layers; layer++){
            if (layer>0 && layer<layers-2){
                for (int i = (layer-1)*num_points+1; i < (layer)*num_points+1; i++){
                    if (i%num_points!=0){
                        fun_arr[point*3] = i;
                        fun_arr[point*3+1] = i+1;
                        fun_arr[point*3+2] = i+num_points+1;
                        fun_arr[point*3+3] = i;
                        fun_arr[point*3+4] = i+num_points;
                        fun_arr[point*3+5] = i+num_points+1;
                        point+=2;
                        }
                    else{
                        fun_arr[point*3] = i;
                        fun_arr[point*3+1] = (layer-1)*num_points+1;
                        fun_arr[point*3+2] = i+num_points;
                        fun_arr[point*3+3] = (layer-1)*num_points+1;
                        fun_arr[point*3+4] = i+num_points;
                        fun_arr[point*3+5] = (layer)*num_points+1;
                        point+=2;
                    }
                }
            }
            else if(layer==0){
                for (int i = 0; i<num_points;i++){
                    if (i!=num_points-1){
                        fun_arr[point*3] = 0;
                        fun_arr[point*3+1] = i+1;
                        fun_arr[point*3+2] = i+2;
                    }
                    else{
                        fun_arr[point*3]=0;
                        fun_arr[point*3+1]=i+1;
                        fun_arr[point*3+2] = 1;
                    }
                    point++;
                }
            }
            else if (layer==layers-2){
                for (int i = num_points*(layer-1)+1; i<num_points*(layer)+1;i++){
                    if (i%num_points!=0){
                        fun_arr[point*3]= num_points*(layer)+1;
                        fun_arr[point*3+1] = i;
                        fun_arr[point*3+2] = 1+i;
                    }
                    else{
                        fun_arr[point*3] = num_points*(layer)+1;
                        fun_arr[point*3+1] = i;
                        fun_arr[point*3+2] = num_points*(layer-1)+1;
                    }
                    point++;
                }
            }
        }
        return fun_arr;
    }
    GLuint* get_all_indices(){
        GLuint* fun_arr = (GLuint*)malloc(size*num_triangles*3*sizeof(GLuint));
        int point = 0;
        int tnumpoints = num_points+1;
        for (int sphere = 0; sphere<size;sphere++){
            for (int layer = 0; layer < layers; layer++){
                if (layer>0 && layer<layers-2){
                    for (int i = (layer-1)*tnumpoints+1; i < (layer)*tnumpoints+1; i++){
                        if (i%tnumpoints!=0){
                            fun_arr[point*3] = i+sphere*num_vertices;
                            fun_arr[point*3+1] = i+1+sphere*num_vertices;
                            fun_arr[point*3+2] = i+tnumpoints+1+sphere*num_vertices;
                            fun_arr[point*3+3] = i+sphere*num_vertices;
                            fun_arr[point*3+4] = i+tnumpoints+sphere*num_vertices;
                            fun_arr[point*3+5] = i+tnumpoints+1+sphere*num_vertices;
                            point+=2;
                            }
                        else{
                            fun_arr[point*3] = i+sphere*num_vertices;
                            fun_arr[point*3+1] = i+1+sphere*num_vertices;
                            fun_arr[point*3+2] = i+tnumpoints+1+sphere*num_vertices;
                            fun_arr[point*3+3] = i+sphere*num_vertices;
                            fun_arr[point*3+4] = i+tnumpoints+sphere*num_vertices;
                            fun_arr[point*3+5] = i+tnumpoints+1+sphere*num_vertices;
                            point+=2;
                        }
                    }
                }
                else if(layer==0){
                    for (int i = 0; i<tnumpoints;i++){
                        if (i!=tnumpoints-1){
                            fun_arr[point*3] = 0+sphere*num_vertices;
                            fun_arr[point*3+1] = i+1+sphere*num_vertices;
                            fun_arr[point*3+2] = i+2+sphere*num_vertices;
                        }
                        else{
                            fun_arr[point*3] = 0+sphere*num_vertices;
                            fun_arr[point*3+1] = i+1+sphere*num_vertices;
                            fun_arr[point*3+2] = i+2+sphere*num_vertices;
                        }
                        point++;
                    }
                }
                else if (layer==layers-2){
                    for (int i = tnumpoints*(layer-1)+1; i<tnumpoints*(layer)+1;i++){
                        if (i%tnumpoints!=0){
                            fun_arr[point*3]= tnumpoints*(layer)+1+sphere*num_vertices;
                            fun_arr[point*3+1] = i+sphere*num_vertices;
                            fun_arr[point*3+2] = 1+i+sphere*num_vertices;
                        }
                        else{
                            fun_arr[point*3]= tnumpoints*(layer)+1+sphere*num_vertices;
                            fun_arr[point*3+1] = i+sphere*num_vertices;
                            fun_arr[point*3+2] = 1+i+sphere*num_vertices;
                        }
                        point++;
                    }
                }
            }
        }
        void* ret = malloc(3);
        return fun_arr;
    }  
    GLfloat* get_sphere_colors(){
        GLfloat* colors = (GLfloat*)malloc(3*num_vertices*sizeof(GLfloat));
        int p = 0;
        for (int i = 0; i < layers; i++){
            for (int j = 0; j<num_points+1; j++ ){

                GLfloat color_val1 = 1.0;
                GLfloat color_val2 = 0.75;
                GLfloat color_val3 = 0.33;
                if (j%3==0){
                    colors[p*3] = 0.2;
                    colors[p*3+1] = 0.2;
                    colors[p*3+2] = 0.8;
                }
                if (j%3==1){
                    colors[p*3] = 0.2;
                    colors[p*3+1] = 0.2;
                    colors[p*3+2] = 0.8;
                }
                if (j%3==2){
                    colors[p*3] = 0.2;
                    colors[p*3+1] = 0.2;
                    colors[p*3+2] = 0.8;
                }
                p++;
                if (i==0 || i==layers-1){
                    break;
                }
            }
        }
        return colors;
    }
    GLfloat* get_all_colors(){
        int p = 0;
        GLfloat* colors = (GLfloat*)malloc(size*num_vertices*sizeof(GLfloat)*3);
        for (int sphere = 0; sphere<size;sphere++){
            for (int i = 0; i < layers; i++){
                for (int j = 0; j<num_points+1; j++ ){
                    GLfloat color_val1 = 1.0;
                    GLfloat color_val2 = 0.75;
                    GLfloat color_val3 = 0.33;
                    if (j%3==0){
                        colors[p*3] = 0.2;
                        colors[p*3+1] = 0.2;
                        colors[p*3+2] = 0.8;
                    }
                    if (j%3==1){
                        colors[p*3] = 0.2;
                        colors[p*3+1] = 0.2;
                        colors[p*3+2] = 0.8;
                    }
                    if (j%3==2){
                        colors[p*3] = 0.2;
                        colors[p*3+1] = 0.2;
                        colors[p*3+2] = 0.8;
                    }
                    p++;
                    if (i==0 || i==layers-1){
                        break;
                    }
                }
            }
        }
        return colors;
    } 
};

particle_system p(250);
GLfloat* sphere_verticies;
GLuint* sphere_elements;
GLfloat* sphere_colors;
GLfloat* norm_vectors;
float aspectaxis(){
    float outputzoom = 1.0f;
    float aspectorigin = 16.0f / 9.0f;
    int aspectconstraint = 1;
    switch(aspectconstraint){
        case 1:
            if ((screen_width/screen_height)<aspectorigin){
                outputzoom *= (((float)screen_width/screen_height)/aspectorigin);
            }
            else{
                outputzoom*= ((float)aspectorigin/aspectorigin);
            }
            break;
        case 2: 
            outputzoom *= (((float)screen_width/screen_height)/aspectorigin);
            break;
        default:
            outputzoom *= ((float)aspectorigin/aspectorigin);
    }
    return outputzoom;
}


float recalculatefov(){
    return 2.0f * glm::atan(glm::tan(glm::radians(45.0f/2.0f))/aspectaxis());
}

void on_resize(int width, int height){
    screen_width = width;
    screen_height = height;
    glViewport(0,0,screen_width,screen_height);
}


glm::mat4 projection = glm::perspective(glm::radians(45.0f), 1.0f*screen_width/screen_height,0.05f,100.0f);
bool init_resources(){

    sphere_elements = p.get_all_indices();
    sphere_colors = p.get_all_colors();
    sphere_verticies = p.generate_spheres();
    norm_vectors = p.get_norm_vectors();


    glGenBuffers(1, &ibo_sphere_elements);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_sphere_elements);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,p.num_triangles*sizeof(GLuint)*3*p.size,sphere_elements,GL_STATIC_DRAW);



    glGenBuffers(1, &vbo_sphere_colors);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_sphere_colors);
    glBufferData(GL_ARRAY_BUFFER, p.num_vertices*3*sizeof(GLfloat)*p.size, sphere_colors, GL_STATIC_DRAW);



    glGenBuffers(1, &vbo_sphere);
    glBindBuffer(GL_ARRAY_BUFFER,vbo_sphere);
    glBufferData(GL_ARRAY_BUFFER,p.num_vertices*3*sizeof(GLfloat)*p.size,sphere_verticies,
                GL_STATIC_DRAW);

    glGenBuffers(1, &vbo_norm_vectors);
    glBindBuffer(GL_ARRAY_BUFFER,vbo_norm_vectors);
    glBufferData(GL_ARRAY_BUFFER,p.num_vertices*3*sizeof(GLfloat)*p.size, norm_vectors,
    GL_STATIC_DRAW);

    GLint compile_ok = GL_FALSE, link_ok = GL_FALSE;
    GLuint vs,fs;
    vs = create_shader("shader.vert",GL_VERTEX_SHADER);
    fs = create_shader("shader.frag",GL_FRAGMENT_SHADER);
    if (vs == 0 || fs == 0){
        return false;
    }
    program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &link_ok);
    glUseProgram(program);

    const GLchar* attribute_name = "coord3d";
    attribute_coord3d = glGetAttribLocation(program,attribute_name);
    if (attribute_coord3d == -1){
        cerr<<"Could not bind the attribute "<<attribute_coord3d<<endl;
        return -1;
    }
    const char* uniform_name;
    uniform_name = "mvp";
    uniform_mvp = glGetUniformLocation(program,uniform_name);
    if (uniform_mvp == -1){
        cerr << "Could not bind uniform "<< uniform_name << endl;
        return 0;
    }
    if (!link_ok){
        cerr << "Error in glLinkProgram" << endl;
        return false;
    }
    uniform_name = "lightpos";
    uniform_lightpos = glGetUniformLocation(program,uniform_name);
    if (uniform_lightpos == -1 ){
        cerr << "Could not bind uniform " << uniform_name <<endl;
        return false;
    }
    uniform_name = "model";
    uniform_model = glGetUniformLocation(program,uniform_name);
    if (uniform_model == -1){
        cerr << "Could not bind uniform " << uniform_name << endl;
        return false;
    }
    attribute_name = "v_color";
    attribute_v_color = glGetAttribLocation(program,attribute_name);
    if (attribute_v_color == -1){
        cerr << "Could not bind the attribute " << attribute_name << endl;
        return false;
    }
    if (attribute_coord3d == -1){
        cerr << "Could not bind the attribute " << attribute_name << endl;
        return false;
    }
    attribute_name = "normvec";
    attribute_normvec = glGetAttribLocation(program,attribute_name);
    if (attribute_normvec == -1){
        cerr << "Could not bind the attribute " << attribute_name << endl;
        return false;
    }
    return true;
}

float facex = -2.0, facey = 0.0, facez = 20.0;
float lookx = 0.0;
float looky = 0.0;
float lookz = 0.0;

GLfloat* create_circle(int num_points,int layers){
    int num_vertices = (layers-2)*num_points+2;
    GLfloat* ret_arr = (GLfloat*)malloc(sizeof(GLfloat)*3*(num_points*(layers-2)+2));
    float radius = 0.02, rotation = 0,step = 2*M_PI/num_points;
    float init_y = -radius;
    int p = 0;
    for (int layer = 0; layer<layers; layer++){
        for (int point = 0; point<num_points;point++){
                float rot = (2*M_PI)/num_points*point;
                float x = sqrt(radius*radius - abs(init_y)*abs(init_y))*cos(rot);
                float z = sqrt(radius*radius - abs(init_y)*abs(init_y))*sin(rot);
                ret_arr[p*3] = (GLfloat)(x+lookx);
                ret_arr[p*3+1] = (GLfloat)(init_y+looky);
                ret_arr[p*3+2] = (GLfloat)(z+lookz);
                if (layer==0||layer==layers-1){
                    ret_arr[p*3]=lookx;
                    ret_arr[p*3+2]=lookz;
                    p++;
                    break;
                }
                else{
                    p++;
                }
            }
            init_y+=(radius*2)/(layers-1);
        } 
    return ret_arr;
}

GLfloat* get_color_array(int num_points, int layers){
    GLfloat* ret_arr = (GLfloat*)malloc(sizeof(GLfloat)*3*(num_points*(layers-2)+2));
    for (int i = 0; i < (((layers-2)*num_points)+2);i++){
        ret_arr[i*3] = 0.0;
        ret_arr[i*3+1] = 0.0;
        ret_arr[i*3+2] = 0.0;
    }
    return ret_arr;
}
int count =0;
GLuint* get_cursor(int num_points, int layers){
    int num_triangles =(layers-3)*(2*num_points)+2*num_points;  
    GLuint* fun_arr = (GLuint*)malloc(3*num_triangles*sizeof(GLuint));
    int point = 0;
    for (int layer = 0; layer < layers; layer++){
        if (layer>0 && layer<layers-2){
            for (int i = (layer-1)*num_points+1; i < (layer)*num_points+1; i++){
                if (i%num_points!=0){
                    fun_arr[point*3] = i;
                    fun_arr[point*3+1] = i+1;
                    fun_arr[point*3+2] = i+num_points+1;
                    fun_arr[point*3+3] = i;
                    fun_arr[point*3+4] = i+num_points;
                    fun_arr[point*3+5] = i+num_points+1;
                    point+=2;
                    }
                else{
                    fun_arr[point*3] = i;
                    fun_arr[point*3+1] = (layer-1)*num_points+1;
                    fun_arr[point*3+2] = i+num_points;
                    fun_arr[point*3+3] = (layer-1)*num_points+1;
                    fun_arr[point*3+4] = i+num_points;
                    fun_arr[point*3+5] = (layer)*num_points+1;
                    point+=2;
                }
            }
        }
        else if(layer==0){
            for (int i = 0; i<num_points;i++){
                if (i!=num_points-1){
                    fun_arr[point*3] = 0;
                    fun_arr[point*3+1] = i+1;
                    fun_arr[point*3+2] = i+2;
                }
                else{
                    fun_arr[point*3]=0;
                    fun_arr[point*3+1]=i+1;
                    fun_arr[point*3+2] = 1;
                }
                point++;
            }
        }
        else if (layer==layers-2){
            for (int i = num_points*(layer-1)+1; i<num_points*(layer)+1;i++){
                if (i%num_points!=0){
                    fun_arr[point*3]= num_points*(layer)+1;
                    fun_arr[point*3+1] = i;
                    fun_arr[point*3+2] = 1+i;
                }
                else{
                    fun_arr[point*3] = num_points*(layer)+1;
                    fun_arr[point*3+1] = i;
                    fun_arr[point*3+2] = num_points*(layer-1)+1;
                }
                point++;
            }
        }
    }
        return fun_arr;
}
void render(SDL_Window* window){
    glClearColor(1.0,1.0,1.0,1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glUseProgram(program);
    glBindBuffer(GL_ARRAY_BUFFER,vbo_sphere);
    glEnableVertexAttribArray(attribute_coord3d);
    glVertexAttribPointer(
        attribute_coord3d,
        3,
        GL_FLOAT,
        GL_FALSE,
        3*sizeof(GLfloat),
        0);
    
    glEnableVertexAttribArray(attribute_v_color);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_sphere_colors);
    glVertexAttribPointer(
        attribute_v_color,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(GLfloat)*3,
        0
    );


    glEnableVertexAttribArray(attribute_normvec);
    glBindBuffer(GL_ARRAY_BUFFER,vbo_norm_vectors);
    glVertexAttribPointer(attribute_normvec,3,GL_FLOAT,GL_FALSE,sizeof(GLfloat)*3,0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_sphere_elements);
    int size; glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER,GL_BUFFER_SIZE, &size);
    glDrawElements(GL_TRIANGLES,p.num_triangles*p.size*3,GL_UNSIGNED_INT,0);


    glBindBuffer(GL_ARRAY_BUFFER,vbo_sphere);
    int num_points = 10,layers = 10;
    GLfloat* cursor_points = create_circle(num_points,layers);


    glBufferData(GL_ARRAY_BUFFER,sizeof(GLfloat)*3*(num_points*(layers-2)+2),cursor_points,GL_STATIC_DRAW);
    free(cursor_points);

    glBindBuffer(GL_ARRAY_BUFFER,vbo_sphere_colors);
    GLfloat* cursor_colors = get_color_array(num_points,layers);
    glBufferData(GL_ARRAY_BUFFER,sizeof(GLfloat)*3*(num_points*(layers-2)+2),cursor_colors,GL_STATIC_DRAW);
    free(cursor_colors);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,ibo_sphere_elements);

    GLuint* cursor_elements = get_cursor(num_points,layers);

    glBufferData(GL_ELEMENT_ARRAY_BUFFER,3*sizeof(GLuint)*((layers-3)*(2*num_points)+2*num_points),cursor_elements,GL_STATIC_DRAW);
    free(cursor_elements);

    glUseProgram(program);
    glBindBuffer(GL_ARRAY_BUFFER,vbo_sphere);
    glEnableVertexAttribArray(attribute_coord3d);
    glVertexAttribPointer(
        attribute_coord3d,
        3,
        GL_FLOAT,
        GL_FALSE,
        3*sizeof(GLfloat),
        0);
    
    glEnableVertexAttribArray(attribute_v_color);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_sphere_colors);
    glVertexAttribPointer(
        attribute_v_color,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(GLfloat)*3,
        0
    );
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_sphere_elements);
    glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER,GL_BUFFER_SIZE, &size);
    glDrawElements(GL_TRIANGLES,((layers-3)*(2*num_points)+2*num_points)*3,GL_UNSIGNED_INT,0);
    // glDrawArrays(GL_TRIANGLES,0,6);

    glDisableVertexAttribArray(attribute_coord3d);
    glDisableVertexAttribArray(attribute_v_color);
    
    SDL_GL_SwapWindow(window);
}

void free_resources(){
    glDeleteProgram(program);
    glDeleteBuffers(1,&vbo_sphere);
    glDeleteBuffers(1,&vbo_sphere_colors);
    glDeleteBuffers(1,&ibo_sphere_elements);
    glDeleteBuffers(1,&vbo_norm_vectors);

}
GLfloat* lightpos = (GLfloat*)malloc(3*sizeof(GLfloat));
void logic(){
    lightpos[0] = 0.0;
    lightpos[1] = 0.0;
    lightpos[2] = 0.0; 
    glm::mat4 model = glm::translate(glm::mat4(1.0f),glm::vec3(0.0,0.0,0.0));
    glm::mat4 view = glm::lookAt(glm::vec3(facex,facey,facez),glm::vec3(lookx,looky,lookz),glm::vec3(0.0,1.0,0.0));
    float angle = SDL_GetTicks() / 1000.0 * 45;
    glm::vec3 axis_y(0,1,0);
    glm::mat4 anim = glm::rotate(glm::mat4(1.0f), glm::radians(angle),axis_y);
    glm::mat4 mvp = projection * view * model;
    float move  = sinf(SDL_GetTicks() / 1000.0 * (2*3.14) / 5) / 2 + 0.5;
    glUniform3fv(uniform_lightpos,1,lightpos);
    glUniformMatrix4fv(uniform_mvp,1,GL_FALSE,glm::value_ptr(mvp));
    glUniformMatrix4fv(uniform_model,1,GL_FALSE,glm::value_ptr(model));
    
    glm::vec3 axis_z(0,0,1);
    glUseProgram(program);
    

}
int mouse_state = 0;
float camera_speed = 100.0;
float divisor = 250.0;
float turn_speed = 0.05;
bool flag = false;
int mousex=0;
int mousey=0;
int prev_x,prev_y;

float look_dist = 2.0, rotation = 0.0,rotation2 = 0.0;
float* current = (float*)malloc(3*sizeof(float));
float* look = (float*)malloc(3*sizeof(float));

float* movement_vector(float* cam_pos, float* look_at,float move_amount){
    float* ret_arr = (float*)malloc(3*sizeof(float));
    float x = abs(cam_pos[0]-look_at[0]);
    float y = abs(cam_pos[1]-look_at[2]);
    float z = abs(cam_pos[2]-look_at[2]);
    float dist = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0));
    ret_arr[0] = move_amount*cos(rotation)*cos(rotation2);
    ret_arr[1] = move_amount * sin(rotation2);
    ret_arr[2] = move_amount * sin(rotation)*cos(rotation2);
    if (cam_pos[0]>look_at[0]){
        cam_pos[0]*=-1;
    }
    if (cam_pos[1]>look_at[1]){
        cam_pos[1]*=-1;
    }
    if (cam_pos[1]>look_at[1]){
        cam_pos[1]*=-1;
    }
    return ret_arr;
}

void main_loop(SDL_Window* window){
    while (true){
        SDL_Event ev;
        while (SDL_PollEvent(&ev)){
            if (ev.type == SDL_QUIT){
                return;
            }
            if (ev.type == SDL_WINDOWEVENT && ev.window.event == SDL_WINDOWEVENT_SIZE_CHANGED){
                on_resize(ev.window.data1,ev.window.data2);
            }
            if (ev.type == SDL_KEYDOWN){
                float* move_arr;
                current[0]=facex,current[1]=facey,current[2]=facez;
                look[0] = lookx,look[1]=looky,look[2]=lookz;
                move_arr = movement_vector(current,look,camera_speed/divisor);
                if (ev.key.keysym.scancode == SDL_SCANCODE_S){
                    facex-=move_arr[0];
                    facey-=move_arr[1];
                    facez-=move_arr[2];
                    lookx-=move_arr[0];
                    looky-=move_arr[1];
                    lookz-=move_arr[2];
                }
                else if (ev.key.keysym.scancode == SDL_SCANCODE_W){
                    facex+=move_arr[0];
                    facey+=move_arr[1];
                    facez+=move_arr[2];
                    lookx+=move_arr[0];
                    looky+=move_arr[1];
                    lookz+=move_arr[2];
                }
                if (ev.key.keysym.scancode == SDL_SCANCODE_RETURN){
                    if (camera_speed<1000){
                        camera_speed+=10;
                    }
                    
                }
                if (ev.key.keysym.scancode == SDL_SCANCODE_RSHIFT){
                    if (camera_speed>0.1){
                        camera_speed-=10;
                    }
                }
                if (ev.key.keysym.scancode == SDL_SCANCODE_Q){
                    if (turn_speed>0.006){
                        turn_speed-=0.005;
                    }
                }
                if (ev.key.keysym.scancode == SDL_SCANCODE_E){
                    if (turn_speed<0.1){
                        turn_speed+=0.005;
                    }
                }
                if (ev.key.keysym.scancode == SDL_SCANCODE_UP){
                    if (p.framerate>6){
                        p.framerate-=5;
                    }
                }
                if (ev.key.keysym.scancode == SDL_SCANCODE_DOWN){
                    p.framerate+=5;
                    if (p.framerate>1000){
                        p.framerate+=20;
                    }
                }
            }
            if (ev.type = SDL_MOUSEBUTTONDOWN){
                float dist = sqrt(pow(facex-lookx,2)+pow(facey-looky,2)+pow(facez-lookz,2));
                if(!flag){
                    mouse_state = SDL_GetMouseState(&mousex,&mousey);
                    rotation = mousex;
                    rotation2 = mousey;
                    lookx = look_dist*cos(rotation)*cos(rotation2)+facex;
                    lookz = look_dist*cos(rotation2)*sin(rotation)+facez;
                    looky = look_dist*sin(rotation2)+facey;
                    flag = true;
                }
                else{
                    if (mouse_state==1){ 
                        rotation+=(mousex-prev_x)*turn_speed;
                        rotation2+=(mousey-prev_y)*turn_speed;
                        lookx = look_dist*cos(rotation)*cos(rotation2)+facex;
                        lookz = look_dist*sin(rotation)*cos(rotation2)+facez;
                        looky = look_dist*sin(rotation2)+facey;
                    }
                    prev_x = mousex;
                    prev_y = mousey;
                }
                mouse_state = SDL_GetMouseState(&mousex,&mousey);
            }
        }
        p.updater();
        free(sphere_verticies);
        sphere_verticies = p.generate_spheres();
        glBindBuffer(GL_ARRAY_BUFFER,vbo_sphere);
        glBufferData(GL_ARRAY_BUFFER,p.num_vertices*3*sizeof(GLfloat)*p.size,sphere_verticies,
                GL_STATIC_DRAW);
        free(sphere_elements);
        sphere_elements = p.get_all_indices();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,ibo_sphere_elements);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,p.num_triangles*3*p.size*sizeof(GLuint),sphere_elements,GL_STATIC_DRAW);
        free(sphere_colors);
        sphere_colors = p.get_all_colors();
        glBindBuffer(GL_ARRAY_BUFFER,vbo_sphere_colors);
        glBufferData(GL_ARRAY_BUFFER,p.num_vertices*3*sizeof(GLfloat)*p.size,sphere_colors,GL_STATIC_DRAW);
        free(norm_vectors);
        norm_vectors = p.get_norm_vectors();
        glBindBuffer(GL_ARRAY_BUFFER,vbo_norm_vectors);
        glBufferData(GL_ARRAY_BUFFER,p.num_vertices*3*sizeof(GLfloat)*p.size, norm_vectors,
        GL_STATIC_DRAW);
        logic();
        render(window);
    }
}


int main(int argc, char* argv[]){
    float* pos = (float*)malloc(3*sizeof(float));
    pos[0]=0.0;pos[1]=0.0;pos[2]=0.0;
    float* vel = (float*)malloc(3*sizeof(float));
    vel[0]=0.0;vel[1]=0.00;vel[2]=0.0;
    float* force = (float*)malloc(3*sizeof(float));
    force[0]=0.0;force[1]=0.0;force[2]=0.0;
    p.add_particle(pos,vel,force,2500000000.0,1);
    vel[0]=0.0;
    vel[1]=0.0;
    pos[0]=1.0;
    pos[1]=1.0;
    pos[2]=1.0;

    p.add_particle(pos,vel,force,45000000000000.0,5.0);
    
    pos[0]=-1.0;
    pos[1]=-1.0;
    pos[2]=-1.0;

    p.add_particle(pos,vel,force,250000000.0,1.0);
    pos[0]=5.0;
    pos[1]=5.0;
    pos[2] = 5.0;
    p.add_particle(pos,vel,force,250000000.0,1.0);

    pos[0]=-5.0;
    pos[1]=-5.0;
    pos[2] =-5.0;
    p.add_particle(pos,vel,force,250000000.0,1.0);

    pos[0]=-10.0;
    pos[1]=-10.0;
    pos[2] =-10.0;
    p.add_particle(pos,vel,force,250000000.0,1.0);

    pos[0]=10.0;
    pos[1]=10.0;
    pos[2] =10.0;
    p.add_particle(pos,vel,force,250000000.0,1.0);
    for (int i = 0;i<20;i++){
        pos[0]=-i*3;
        pos[1]=i*3;
        pos[2]=50;
        p.add_particle(pos,vel,force,250000000.0,1.0);
    }
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("My Fun N-Body Simulation",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, screen_width, screen_height,
        SDL_WINDOW_RESIZABLE | SDL_WINDOW_OPENGL);
    SDL_GL_CreateContext(window);

    GLenum glew_status = glewInit();
    if (glew_status != GLEW_OK){
        cerr << "Error: glewInit" << glewGetErrorString(glew_status)
        << endl;
        return EXIT_FAILURE;
    }
    if (!init_resources()){
        return EXIT_FAILURE;
    }
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE,1);
    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    main_loop(window);
    free_resources();
    free(p.velocities);
    free(p.forces);
    free(p.positions);
    free(p.radii);
    free(p.masses);
    free(p.light_pos);
    free(sphere_colors);
    free(sphere_elements);
    free(sphere_verticies);

    return EXIT_SUCCESS;
}

