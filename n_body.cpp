// F: Represents the force of gravity between two objects. 
// G: Is the universal gravitational constant, approximately 6.674 × 10^-11 N(m²/kg²). 
// m1 and m2: Are the masses of the two objects involved in the gravitational interaction. 
// r: Is the distance between the centers of the two objects. 
// r²: The distance is squared, meaning that the force of gravity decreases rapidly as the objects move farther apart. 
// G(m1m2)/r^2
#include <cstdlib>
#include <iostream>
using namespace std;
#include <fstream>
#include <math.h>

#include <SDL2/SDL_image.h>
#include <SDL2/SDL.h>
#include "shader_utils.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


class particle_system{
    public:
        float grav_const = 6.674 * pow(10.0,11.0);
        int layers = 25;
        int num_points = 99;
        int size = 0;
        float framerate;
        float* positions;
        float* velocities;
        float* forces;
        float* radii;
        float* masses;
        particle_system(float frames){
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
            size++;
            float* temp_pos = (float*)malloc(sizeof(positions)+3*sizeof(float));
            float* temp_vel = (float*)malloc(sizeof(velocities)+3*sizeof(float));
            float* temp_force = (float*)malloc(sizeof(forces)+3*sizeof(float));
            float* temp_mass = (float*)malloc(sizeof(masses)+sizeof(float));
            float* temp_radius = (float*)malloc(sizeof(radii)+sizeof(float));

            for (int i = 0; i<sizeof(positions)/sizeof(float)/3;i++){
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
            temp_pos[sizeof(positions)/sizeof(float)]=position[0];
            temp_pos[sizeof(positions)/sizeof(float)+1]=position[1];
            temp_pos[sizeof(positions)/sizeof(float)+2]=position[2];

            temp_vel[sizeof(velocities)/sizeof(float)]=velocity[0];
            temp_vel[sizeof(velocities)/sizeof(float)+1]=velocity[1];
            temp_vel[sizeof(velocities)/sizeof(float)+2]=velocity[2];
            temp_force[sizeof(forces)/sizeof(float)]=force[0];
            temp_force[sizeof(forces)/sizeof(float)+1]=force[1];
            temp_force[sizeof(forces)/sizeof(float)+2]=force[2];

            temp_mass[sizeof(masses)/sizeof(float)] = mass;
            temp_radius[sizeof(radii)/sizeof(float)] = radius;
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
        }
        void update_positions(){
            for (int i = 0; i<sizeof(velocities)/sizeof(float);i++){
                positions[i]+=velocities[i]/framerate;
            }
        }
        void update_velocities(){
            for (int i = 0; i<sizeof(forces)/sizeof(float);i++){
                velocities[i]+=(forces[i]/framerate)/masses[i/3];
            }
        }
        float gravity_equation(float mass1,float mass2, float distance){
            float num = grav_const*mass1*mass2;
            float denom = pow(distance,2.0);
            return num/denom;
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
            float theta1 = atan2((double)z,(double)x);
            float theta2 = atan2((double)y,(double)x);
            float force = gravity_equation(masses[ind1],masses[ind2],distance_calculator(ind1,ind2));
            float force_x = pow(tan(theta1),3.0) + sqrt(pow(tan(theta1),6.0)-4*pow(force,2.0)*(pow(tan(theta1),2.0)-pow(cos(theta1),2.0)));
            force_x/=2.0;
            float force_y = force_x*tan(theta1);
            float force_z = force_x*tan(theta2);
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
            for (int i = 0; i<size;i++){
                for (int j = i+1; j<size;j++){
                    calc_force(i,j);
                }
            }
        }
        void merge_masses(int ind1, int ind2){
            float volume = 16.0/9.0 * M_PI*M_PI * radii[ind1]
            * radii[ind1] + radii[ind2]*radii[ind2];  
            float new_radius = pow(volume*3.0/4.0/M_PI,1/3);
            float new_mass = masses[ind1]+masses[ind2];
            float velx = (masses[ind1]*velocities[ind1*3]+ velocities[ind2*3]* masses[ind2])/new_mass;
            float vely = (masses[ind1]*velocities[ind1*3+1]+ velocities[ind2*3+1]* masses[ind2])/new_mass;
            float velz = (masses[ind1]*velocities[ind1*3+2]+ velocities[ind2*3+2]* masses[ind2])/new_mass;
            float forcex = forces[ind1*3]+forces[ind2*3];
            float forcey = forces[ind1*3+1]+forces[ind2*3+1];
            float forcez = forces[ind1*3+2]+forces[ind2*3+2];

            //m1v1 + m2v2 = mv




            float* temp_pos = (float*)malloc(sizeof(positions)-3*sizeof(float));
            float* temp_vel = (float*)malloc(sizeof(velocities)-3*sizeof(float));
            float* temp_force = (float*)malloc(sizeof(forces)-3*sizeof(float));
            float* temp_mass = (float*)malloc(sizeof(masses)-sizeof(float));
            float* temp_radius = (float*)malloc(sizeof(radii)-sizeof(float));

            for (int i = 0; i<size;i++){
                if (i!=ind1 && i!=ind2){
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
            }

            temp_pos[size*sizeof(float)*3-1] = positions[ind1*3];
            temp_pos[size*sizeof(float)*3+1-1] = positions[ind1*3+1];
            temp_pos[size*sizeof(float)*3+2-1] = positions[ind1*3+2];
            
            temp_vel[size*sizeof(float)*3-1] = velx;
            temp_vel[size*sizeof(float)*3+1-1] = vely;
            temp_vel[size*sizeof(float)*3+2-1] = velz;
            
            temp_force[size*sizeof(float)*3-1] = forcex;
            temp_force[size*sizeof(float)*3+1-1] = forcey;
            temp_force[size*sizeof(float)*3+2-1] = forcez;
            
            temp_radius[size-1] = new_radius;
            temp_mass[size-1] = new_mass;  
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
        return distance<=(radii[ind1]+radii[ind2])/3.0;
    }
    void update_merges(){
        for (int i = 0; i < size; i++){
            for (int j = 1; j < size; j++){
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
        GLfloat* ret_arr = (GLfloat*)malloc(num_points*layers*3*sizeof(GLfloat));
        int p = 0;
        for (int layer = 0; layer<layers; layer++){
            for (int point = 0; point<num_points;point++){
                float rot = num_points/2*M_PI*point;
                float x = radii[ind]*cos(rot);
                float z = radii[ind]*sin(rot);
                float y = radii[ind]*2/layers*layer;
                ret_arr[p] = (GLfloat)x;
                ret_arr[p+1] = (GLfloat)y;
                ret_arr[p+2] = (GLfloat)z;
            }
        } 
        return ret_arr;
    }
    GLint* get_indexes(){
        //10*3 for 8, 
        GLint* fun_arr = (GLint*)malloc(sizeof(GLint)*(layers-1)*(num_points+1)*6);
        for (int layer = 0; layer < layers; layer++){
            if (layer!=layers-1){
                for (int i = layer*num_points; i < layer*num_points+num_points; i++){
                    if (i!=num_points-1){
                        fun_arr[i*3] = i;
                        fun_arr[i*3+1] = i+1;
                        fun_arr[i*3+2] = i+num_points+1;
                        fun_arr[i*3+3] = i;
                        fun_arr[i*3+4] = i+num_points;
                        fun_arr[i*3+5] = i+num_points+1;
                        }
                    else{
                        fun_arr[i*3] = i;
                        fun_arr[i*3+1] = 0;
                        fun_arr[i*3+2] = i+1;
                    }
                }
            }
        }
        return fun_arr;
    }   
};


int main(){
    cout << "Hello world" << endl;
    return 1;
}
//F = ma
//a = F/m
//v = v0 + at
//v = v0 + Ft/m
//V = 4/3pir^3