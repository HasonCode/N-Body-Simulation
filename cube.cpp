#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;

#include <GL/glew.h>
GLuint program;
GLint attribute_coord2d, attribute_v_color;
GLint attribute_coord3d;
GLuint vbo_triangle, vbo_triangle_colors;
GLint uniform_fade;
GLuint vbo_cube, vbo_cube_colors;
GLuint ibo_cube_elements;
GLint uniform_mvp;
#include <SDL2/SDL_image.h>
#include <SDL2/SDL.h>
#include "shader_utils.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
// #include "SDL_image.h"
GLint uniform_m_transform;

int screen_width = 800, screen_height=600;
struct attributes{
    GLfloat coord3d[3];
    GLfloat v_color[3];
};


class particle_system{
    public:
        float grav_const = 6.674 * pow(10.0,-11.0);
        int layers = 50;
        int num_points = 50;
        float bounce_dist = 50.0;
        float bounce_return = 0.9;
        int num_vertices = (layers-2)*num_points+2;
        int num_triangles =(layers-3)*(2*num_points)+2*num_points;
        int size = 0;
        float framerate;
        float* positions;
        float* velocities;
        float* forces;
        float* radii;
        float* masses;
        particle_system(float frames)
        {
            // layers = num_layers;
            // n_points = n_points;
            // num_vertices = (layers-2)*num_points+2;
            // num_triangles =(layers-3)*(2*num_points)+2*num_points;
            framerate=frames;
            positions = (float*)malloc(0*sizeof(float));
            velocities = (float*)malloc(0*sizeof(float));
            forces = (float*)malloc(0*sizeof(float));
            radii = (float*)malloc(0*sizeof(float));
            masses = (float*)malloc(0*sizeof(float));

        }
        float distance_calculator(int ind1, int ind2){
            float distancex,distancey,distancez;
            distancex = pow(positions[ind1*3]-positions[ind2*3],2);
            distancey = pow(positions[ind1*3+1]-positions[ind2*3+2],2);
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
            // cout<<radii[size]<<endl;
            // if(size==1){
            //     cout<<radii[size-1]<<endl;
            // }
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
                    // cout<<"Bounce Pos: "<<positions[i*3] << " Bounce Vel: "<< velocities[i*3]<<endl;
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
            //theta1 = atan(z/x)
            //theta2 = atan(y/x)
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
            float theta1 = atan2(distance_calculator(ind1,ind2),y);
            float theta2 = acosf((double)x/(double)distance_calculator(ind1,ind2));
            float force = gravity_equation(masses[ind1],masses[ind2],distance_calculator(ind1,ind2));
            float force_x = force*sin(theta1)*cos(theta2);
            float force_y = force*sin(theta1);
            float force_z = force*sin(theta1)*sin(theta2);
            // cout<<"Force X: "<<force_x<< " Force Y " << force_y << " Force Z "<<force_z<<endl;
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
            //m1v1 + m2v2 = mv

            float posx,posy,posz;
            if (masses[ind1]>masses[ind2]){
                posx = positions[ind1*3];
                posy = positions[ind1*3+1];
                posz = positions[ind1*3+2];
            }
            else if (masses[ind1]<masses[ind2]){
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
            // cout<<(size-1)<<endl;
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
                    // cout<<"interesting..."<<endl;
                    count++;
                }
            }
            // cout<<"PRE MERGER:"<<"X: " << positions[ind1]<< " Y " << positions[ind1+1] << " Z " << positions[ind1+2]<<endl;
            // cout<<"PRE MERGER:"<<"X: " << positions[ind2]<< " Y " << positions[ind2+1] << " Z " << positions[ind2+2]<<endl;
            // cout<<"PRE MERGER:"<<"X: " << (positions[ind1]+positions[ind2])/2.0<< " Y " << (positions[ind1+1]+positions[ind2+1])/2.0 << " Z " << (positions[ind1+2]+positions[ind2+2])/2.0<<endl;


            temp_pos[(size-1)*3-3] = (positions[ind1*3]+positions[ind2*3])/2.0;
            temp_pos[(size-1)*3-2] = (positions[ind1*3+1]+positions[ind2*3+1])/2.0;
            temp_pos[(size-1)*3-1] = (positions[ind1*3+2]+positions[ind2*3+2])/2.0;
            // cout<<size*3-3<<endl;

            temp_vel[(size-1)*3-3] = velx;
            temp_vel[(size-1)*3-2] = vely;
            temp_vel[(size-1)*3-1] = velz;
            
            temp_force[(size-1)*3-3] = forcex;
            temp_force[(size-1)*3-2] = forcey;
            temp_force[(size-1)*3-1] = forcez;
            
            temp_radius[size-2] = new_radius;
            temp_mass[size-2] = new_mass; 
            // cout<<"POST MERGER:"<<"X: " << temp_pos[0]<< " Y " << temp_pos[1] << " Z " << temp_pos[2]<<endl;
            size--;
            free(positions);
            free(velocities);
            free(forces);
            free(masses);
            free(radii);
            positions = temp_pos;
            // cout<<"POST MERGER:"<<"X: " << positions[0]<< " Y " << positions[1] << " Z " << positions[2]<<endl;
            velocities = temp_vel;
            forces = temp_force;
            masses = temp_mass;
            radii = temp_radius;
        }
    bool can_merge(int ind1, int ind2){
        float distance = distance_calculator(ind1,ind2);
        return distance<=(radii[ind1]);
    }
    void update_merges(){
        for (int i = 0; i < size; i++){
            for (int j = i+1; j < size; j++){
                if (can_merge(i,j)){
                    merge_masses(i,j);
                    i=0;
                    break;
                }
            }
        }
    }
    void updater(){
        update_positions();
        update_velocities();
        update_forces();
        update_merges();
    }

    GLfloat* generate_sphere(int ind){
        GLfloat* ret_arr = (GLfloat*)malloc(num_vertices*3*sizeof(GLfloat));
        int p = 0;
        float init_y = -radii[ind];
        float radius = radii[ind];
        for (int layer = 0; layer<layers; layer++){
            for (int point = 0; point<num_points;point++){
                // cout<<p<<endl;
                // float rot = num_points/2*M_PI*point;
                float rot = (2*M_PI)/num_points*point;
                // cout << "Rot: " << rot << endl;
                float x = sqrt(radius*radius - abs(init_y)*abs(init_y))*cos(rot);
                float z = sqrt(radius*radius - abs(init_y)*abs(init_y))*sin(rot);
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
                for (int point = 0; point<num_points;point++){
                    // cout<<p<<endl;
                    // float rot = num_points/2*M_PI*point;
                    float rot = (2*M_PI)/num_points*point;
                    // cout << "Rot: " << rot << endl;
                    float x = sqrt(radius*radius - abs(init_y)*abs(init_y))*cos(rot);
                    float z = sqrt(radius*radius - abs(init_y)*abs(init_y))*sin(rot);
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
        return ret_arr;
    }
    //Should be good now, plug into GPT for some feedback though
    GLushort* get_indexes(){
        //10*3 for 8, 
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
    GLushort* get_all_indices(){
        GLushort* fun_arr = (GLushort*)malloc(size*num_triangles*3*sizeof(GLushort));
        int point = 0;
        for (int sphere = 0; sphere<size;sphere++){
            for (int layer = 0; layer < layers; layer++){
                if (layer>0 && layer<layers-2){
                    for (int i = (layer-1)*num_points+1; i < (layer)*num_points+1; i++){
                        if (i%num_points!=0){
                            fun_arr[point*3] = i+sphere*num_vertices;
                            fun_arr[point*3+1] = i+1+sphere*num_vertices;
                            fun_arr[point*3+2] = i+num_points+1+sphere*num_vertices;
                            fun_arr[point*3+3] = i+sphere*num_vertices;
                            fun_arr[point*3+4] = i+num_points+sphere*num_vertices;
                            fun_arr[point*3+5] = i+num_points+1+sphere*num_vertices;
                            point+=2;
                            }
                        else{
                            fun_arr[point*3] = i+sphere*num_vertices;
                            fun_arr[point*3+1] = (layer-1)*num_points+1+sphere*num_vertices;
                            fun_arr[point*3+2] = i+num_points+sphere*num_vertices;
                            fun_arr[point*3+3] = (layer-1)*num_points+1+sphere*num_vertices;
                            fun_arr[point*3+4] = i+num_points+sphere*num_vertices;
                            fun_arr[point*3+5] = (layer)*num_points+1+sphere*num_vertices;
                            point+=2;
                        }
                    }
                }
                else if(layer==0){
                    for (int i = 0; i<num_points;i++){
                        if (i!=num_points-1){
                            fun_arr[point*3] = 0+sphere*num_vertices;
                            fun_arr[point*3+1] = i+1+sphere*num_vertices;
                            fun_arr[point*3+2] = i+2+sphere*num_vertices;
                        }
                        else{
                            fun_arr[point*3]=0+sphere*num_vertices;
                            fun_arr[point*3+1]=i+1+sphere*num_vertices;
                            fun_arr[point*3+2] = 1+sphere*num_vertices;
                        }
                        point++;
                    }
                }
                else if (layer==layers-2){
                    for (int i = num_points*(layer-1)+1; i<num_points*(layer)+1;i++){
                        if (i%num_points!=0){
                            fun_arr[point*3]= num_points*(layer)+1+sphere*num_vertices;
                            fun_arr[point*3+1] = i+sphere*num_vertices;
                            fun_arr[point*3+2] = 1+i+sphere*num_vertices;
                        }
                        else{
                            fun_arr[point*3] = num_points*(layer)+1+sphere*num_vertices;
                            fun_arr[point*3+1] = i+sphere*num_vertices;
                            fun_arr[point*3+2] = num_points*(layer-1)+1+sphere*num_vertices;
                        }
                        point++;
                    }
                }
            }
        }
        return fun_arr;
    }  
    GLfloat* get_sphere_colors(){
        // cout<<"Vertices: "<<num_vertices<<endl;
        GLfloat* colors = (GLfloat*)malloc(3*num_vertices*sizeof(GLfloat));
        int p = 0;
        for (int i = 0; i < layers; i++){
            for (int j = 0; j<num_points; j++ ){

                GLfloat color_val1 = 1.0;
                GLfloat color_val2 = 0.75;
                GLfloat color_val3 = 0.33;
                if (j%3==0){
                    colors[p*3] = 1.0;
                    colors[p*3+1] = 0.0;
                    colors[p*3+2] = 0.0;
                }
                if (j%3==1){
                    colors[p*3] = 0.0;
                    colors[p*3+1] = 1.0;
                    colors[p*3+2] = 0.0;
                }
                if (j%3==2){
                    colors[p*3] = 0.0;
                    colors[p*3+1] = 0.0;
                    colors[p*3+2] = 1.0;
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
                for (int j = 0; j<num_points; j++ ){
                    GLfloat color_val1 = 1.0;
                    GLfloat color_val2 = 0.75;
                    GLfloat color_val3 = 0.33;
                    if (j%3==0){
                        colors[p*3] = 1.0;
                        colors[p*3+1] = 0.0;
                        colors[p*3+2] = 0.0;
                    }
                    if (j%3==1){
                        colors[p*3] = 0.0;
                        colors[p*3+1] = 1.0;
                        colors[p*3+2] = 0.0;
                    }
                    if (j%3==2){
                        colors[p*3] = 0.0;
                        colors[p*3+1] = 0.0;
                        colors[p*3+2] = 1.0;
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

particle_system p(10);
GLfloat* sphere_verticies;
GLushort* sphere_elements;


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

    GLfloat triangle_verticies[] = {
        0.0, 0.0,
        0.5, 0.0,
        0.5, 0.3,
 
        0.0, 0.0,
        0.0, 0.3,
        0.5, 0.3

    };

    GLfloat triangle_colors[] = {
        1.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 0.0, 0.0,
    };

    GLfloat cube_verticies[] = {
        -1.0, -1.0, 1.0,
        1.0, -1.0, 1.0, 
        1.0, 1.0, 1.0, 
        -1.0, 1.0, 1.0,

        -1.0, -1.0, -1.0,
        1.0, -1.0, -1.0,
        1.0, 1.0, -1.0, 
        -1.0, 1.0, -1.0
    };

    GLfloat cube_colors[] = {
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0, 
        1.0, 1.0, 1.0,

        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 1.0, 1.0
    };
    struct attributes triangle_attributes[] = {
        {{0.0,  0.8,0.0},   {1.0, 1.0, 0.0,}},
        {{-0.8, -0.8,0.0},   {0.0, 0.0, 1.0,}},
         {{0.8, -0.8,0.0 },  {1.0, 0.0, 0.0,}}
    };

    GLushort cube_elements[] = {
        // front
        0,1,2,
        2,3,0,
        // right 
        1, 5, 6,
        6, 2, 1,
        // back
        7, 6, 5,
        5, 4, 7,
        // left
        4, 0, 3, 
        3, 7, 4,
        //bottom 
        4, 5, 1,
        1, 0, 4,
        // top 
        3, 2, 6, 
        6, 7, 3

    };
    GLfloat temp_sphere_colors[] = {
        0.0,1.0,0.0,
        1.0,0.0,0.0,
        0.0,1.0,0.0,
        0.0,0.0,1.0,
        1.0,0.0,0.0,
        1.0,0.0,0.0,
        0.0,1.0,0.0,
        0.0,0.0,1.0,
        1.0,0.0,0.0,
        1.0,0.0,0.0,
        0.0,1.0,0.0,
        0.0,0.0,1.0,
        1.0,0.0,0.0,
        1.0,0.0,0.0,
        0.0,1.0,0.0
    };
    GLshort temp_sphere_elements[] = {
        0,1,2,
        0,2,3,
        0,3,4,
        0,4,1,
        1,2,6,
        1,5,6,
        2,3,7,
        2,6,7,
        3,4,8,
        3,7,8,
        4,1,8,
        1,8,5,
        5,6,10,
        5,9,10,
        6,7,11,
        6,10,11,
        7,8,12,
        7,11,12,
        8,5,12,
        5,12,9,
        13,9,10,
        13,10,11,
        13,11,12,
        13,12,9,
    };

    GLfloat temp_sphere_verticies[] = {
        0,-1,0,
        0.866, -0.5, 0,
        0, -0.5,0.866,
        -0.866,-0.5,0,
        0,-0.5,-0.866,
        1,0,0,
        0,0,1,
        -1,0,0,
        0,0,-1,
        0.866,0.5,0,
        0,0.5,0.866,
        -0.866,0.5,0,
        0,0.5,-0.866,
        0,1,0,
    };

    sphere_elements = p.get_all_indices();
    GLfloat* sphere_colors = p.get_all_colors();
    sphere_verticies = p.generate_spheres();


    // GLfloat* sphere_colors = temp_sphere_colors;
    // GLfloat* sphere_verticies = temp_sphere_verticies;
    // for (int i = 0; i < p.num_vertices*p.size; i++){
    //     cout << "X: " << sphere_verticies[i*3] << " Y: " << sphere_verticies[i*3+1] << " Z: " << sphere_verticies[i*3+2]<<endl;
    // }

    // for (int i = 0; i < p.num_vertices; i++){
    //     cout << "R: " << sphere_colors[i*3] << " G: " << sphere_colors[i*3+1] << " B: " << sphere_colors[i*3+2]<<endl;
    // }


    // for (int i = 0; i < p.num_triangles*2; i++){
    //     cout << "A: " << sphere_elements[i*3] << " B: " << sphere_elements[i*3+1] << " C: " << sphere_elements[i*3+2]<<endl;
    // }
    glGenBuffers(1, &ibo_cube_elements);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_cube_elements);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,p.num_triangles*sizeof(GLshort)*3*p.size,sphere_elements,GL_STATIC_DRAW);



    glGenBuffers(1, &vbo_cube_colors);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_cube_colors);
    glBufferData(GL_ARRAY_BUFFER, p.num_vertices*3*sizeof(GLfloat)*p.size, sphere_colors, GL_STATIC_DRAW);



    glGenBuffers(1, &vbo_cube);
    glBindBuffer(GL_ARRAY_BUFFER,vbo_cube);
    glBufferData(GL_ARRAY_BUFFER,p.num_vertices*3*sizeof(GLfloat)*p.size,sphere_verticies,
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
    // uniform_name = "fade";
    // uniform_fade = glGetUniformLocation(program,uniform_name);
    // if (uniform_fade == -1){
    //     cerr << "Could not bind uniform " << uniform_name << endl;
    //     return false;
    // }
    uniform_name = "mvp";
    uniform_mvp = glGetUniformLocation(program,uniform_name);
    if (uniform_mvp == -1){
        cerr << "Could not bind uniform "<< uniform_name << endl;
        return 0;
    }
    const char* uniform_name2;
    // uniform_name2 = "m_transform";
    // uniform_m_transform = glGetUniformLocation(program,uniform_name2);
    // if (uniform_m_transform == -1){
    //     cerr << "Could not bind uniform " << uniform_name2 << endl;
    //     return false;
    // }
    if (!link_ok){
        cerr << "Error in glLinkProgram" << endl;
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
    return true;
}
void render(SDL_Window* window){
    glClearColor(1.0,1.0,1.0,1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glUseProgram(program);
    glBindBuffer(GL_ARRAY_BUFFER,vbo_cube);
    glEnableVertexAttribArray(attribute_coord3d);
    glVertexAttribPointer(
        attribute_coord3d,
        3,
        GL_FLOAT,
        GL_FALSE,
        3*sizeof(GLfloat),
        0);
    
    glEnableVertexAttribArray(attribute_v_color);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_cube_colors);
    glVertexAttribPointer(
        attribute_v_color,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(GLfloat)*3,
        0
    );
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_cube_elements);
    int size; glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER,GL_BUFFER_SIZE, &size);
    glDrawElements(GL_TRIANGLES,p.num_triangles*p.size*3,GL_UNSIGNED_SHORT,0);



    // glDrawArrays(GL_TRIANGLES,0,6);

    glDisableVertexAttribArray(attribute_coord3d);
    glDisableVertexAttribArray(attribute_v_color);
    
    SDL_GL_SwapWindow(window);
}

void free_resources(){
    glDeleteProgram(program);
    glDeleteBuffers(1,&vbo_triangle);
}
float facex = -2.0, facey = 0.0, facez = 20.0;
void logic(){
    // float current_fade  = sinf(SDL_GetTicks() / 1000.0 * (2*3.14) / 5) / 2 + 0.5;
    // glm::mat4 projection = glm::mat4(1.0f); 
    glm::mat4 model = glm::translate(glm::mat4(1.0f),glm::vec3(0.0,0.0,-4.0));
    glm::mat4 view = glm::lookAt(glm::vec3(facex,facey,facez),glm::vec3(0.0,0.0,0.0),glm::vec3(0.0,1.0,0.0));
    float angle = SDL_GetTicks() / 1000.0 * 45;
    glm::vec3 axis_y(0,1,0);
    glm::mat4 anim = glm::rotate(glm::mat4(1.0f), glm::radians(angle),axis_y);
    glm::mat4 mvp = projection * view * model;
    float move  = sinf(SDL_GetTicks() / 1000.0 * (2*3.14) / 5) / 2 + 0.5;
    glUniformMatrix4fv(uniform_mvp,1,GL_FALSE,glm::value_ptr(mvp));
    glm::vec3 axis_z(0,0,1);
    // glm::mat4 m_transform = glm::translate(glm::mat4(1.0f),glm::vec3(move,0.0,0.0))* glm::rotate(glm::mat4(1.0f),glm::radians(angle),axis_z);
    glUseProgram(program);
    // glUniform1f(uniform_fade,current_fade);
    // glUniformMatrix4fv(uniform_m_transform,1,GL_FALSE,glm::value_ptr(m_transform));
    // glUniformMatrix4fv(uniform_mvp,1,GL_FALSE,glm::value_ptr(mvp));

}
int mouse_state = 0;
float camera_speed = 100.0;
float divisor = 250.0;
float turn_speed = 0.05;
bool flag = false;
int mousex=0;
int mousey=0;
int prev_x,prev_y;

float* movement_vector(float* cam_pos, float* look_at,float move_amount){
    float* ret_arr = (float*)malloc(3*sizeof(float));
    float x = abs(cam_pos[0]-look_at[0]);
    float y = abs(cam_pos[1]-look_at[2]);
    float z = abs(cam_pos[2]-look_at[2]);
    float dist = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0));
    ret_arr[0] = move_amount*sin(acos(y/dist))*cos(atan2(z,x));
    ret_arr[1] = move_amount * sin(acos(y/dist))*cos(atan2(z,x));
    ret_arr[2] = move_amount * cos(acos(y/dist));
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
float lookx = 0.0;
float looky = 0.0;
float lookz = -4.0;
float* current = (float*)malloc(3*sizeof(float));
float* look = (float*)malloc(3*sizeof(float));
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
                if (ev.key.keysym.scancode == SDL_SCANCODE_W){
                    facex-=move_arr[0];
                    facey-=move_arr[1];
                    facez-=move_arr[2];
                }
                // else if (ev.key.keysym.scancode == SDL_SCANCODE_A){
                //     facex -= camera_speed/divisor;
                // }
                else if (ev.key.keysym.scancode == SDL_SCANCODE_S){
                    facex+=move_arr[0];
                    facey+=move_arr[1];
                    facez+=move_arr[2];
                }
                
                if (ev.key.keysym.scancode == SDL_SCANCODE_RETURN){
                    if (camera_speed<1000){
                        camera_speed+=10;
                    }
                }
                if (ev.key.keysym.scancode == SDL_SCANCODE_RSHIFT){
                    if (camera_speed>10){
                        camera_speed-=10;
                    }
                }
                // else if (ev.key.keysym.scancode == SDL_SCANCODE_D){
                //     facez-=camera_speed/divisor;
                //     cout<<"D is down and facez has a value of: "<<facez<<endl;
                // }
            }
            if (ev.type = SDL_MOUSEBUTTONDOWN){
                if(!flag){
                    mouse_state = SDL_GetMouseState(&mousex,&mousey);
                    prev_x = mousex;
                    prev_y = mousey;
                    flag = true;
                }
                else{
                    if (mouse_state==1){ 
                        facex+=(mousex-prev_x)*turn_speed;
                        facey+=(mousey-prev_y)*turn_speed;
                    }
                    prev_x = mousex;
                    prev_y = mousey;
                }
                mouse_state = SDL_GetMouseState(&mousex,&mousey);
            }
            // if (ev.type == SDL_MOUSEMOTION){
            //     facex-=ev.motion.xrel/camera_speed;
            //     facey-=ev.motion.yrel/camera_speed;
            // }

        }
        // p.update_forces();
        cout<<"X position: " << p.positions[0] << " Y position: " << p.positions[1] << " Z position: " << p.positions[2]<< endl;
        
        // cout<<"X force: " << p.forces[0] << " Y force: " << p.forces[1] << " Z force: " << p.forces[2]<< endl;
        // cout<<"X velocity: " << p.velocities[3] << " Y velocity: " << p.velocities[4] << " Z velocity: " << p.velocities[5]<< endl;
        p.updater();
        free(sphere_verticies);
        sphere_verticies = p.generate_spheres();
        glBindBuffer(GL_ARRAY_BUFFER,vbo_cube);
        glBufferData(GL_ARRAY_BUFFER,p.num_vertices*3*sizeof(GLfloat)*p.size,sphere_verticies,
                GL_STATIC_DRAW);
        sphere_elements = p.get_all_indices();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,ibo_cube_elements);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,p.num_triangles*3*p.size*sizeof(GLushort),sphere_elements,GL_STATIC_DRAW);
        logic();
        // cout << "Pointer to window:" << window<<endl;
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

    p.add_particle(pos,vel,force,45000000000.0,1.3);
    
    pos[0]=-1.0;
    pos[1]=-1.0;
    pos[2]=-1.0;
    p.add_particle(pos,vel,force,250000000.0,1.0);
    // cout<< "Velocity: "<<p.velocities[8]<<" Vel[2]"<<vel[2]<<endl;
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("My Textured Cube",
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
    // glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    main_loop(window);
    free_resources();
    return EXIT_SUCCESS;
}

