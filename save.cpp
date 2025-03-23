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

particle_system p(60);




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


glm::mat4 projection = glm::perspective(recalculatefov(), 1.0f*screen_width/screen_height,0.1f,10.0f);
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

    glGenBuffers(1, &ibo_cube_elements);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_cube_elements);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(cube_elements),cube_elements,GL_STATIC_DRAW);



    glGenBuffers(1, &vbo_cube_colors);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_cube_colors);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cube_colors), cube_colors, GL_STATIC_DRAW);



    glGenBuffers(1, &vbo_cube);
    glBindBuffer(GL_ARRAY_BUFFER,vbo_cube);
    glBufferData(GL_ARRAY_BUFFER,sizeof(cube_verticies),cube_verticies,
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
    glDrawElements(GL_TRIANGLES,size/sizeof(GLushort),GL_UNSIGNED_SHORT,0);



    // glDrawArrays(GL_TRIANGLES,0,6);

    glDisableVertexAttribArray(attribute_coord3d);
    glDisableVertexAttribArray(attribute_v_color);
    
    SDL_GL_SwapWindow(window);
}

void free_resources(){
    glDeleteProgram(program);
    glDeleteBuffers(1,&vbo_triangle);
}

void logic(){
    // float current_fade  = sinf(SDL_GetTicks() / 1000.0 * (2*3.14) / 5) / 2 + 0.5;
    // glm::mat4 projection = glm::mat4(1.0f); 
    glm::mat4 model = glm::translate(glm::mat4(1.0f),glm::vec3(0.0,0.0,-4.0));
    glm::mat4 view = glm::lookAt(glm::vec3(0.0,2.0,0.0),glm::vec3(0.0,0.0,-4.0),glm::vec3(0.0,1.0,0.0));
    float angle = SDL_GetTicks() / 1000.0 * 45;
    glm::vec3 axis_y(0,1,0);
    glm::mat4 anim = glm::rotate(glm::mat4(1.0f), glm::radians(angle),axis_y);
    glm::mat4 mvp = projection * view * model * anim;
    float move  = sinf(SDL_GetTicks() / 1000.0 * (2*3.14) / 5) / 2 + 0.5;
    glUniformMatrix4fv(uniform_mvp,1,GL_FALSE,glm::value_ptr(mvp));
    glm::vec3 axis_z(0,0,1);
    // glm::mat4 m_transform = glm::translate(glm::mat4(1.0f),glm::vec3(move,0.0,0.0))* glm::rotate(glm::mat4(1.0f),glm::radians(angle),axis_z);
    glUseProgram(program);
    // glUniform1f(uniform_fade,current_fade);
    // glUniformMatrix4fv(uniform_m_transform,1,GL_FALSE,glm::value_ptr(m_transform));
    // glUniformMatrix4fv(uniform_mvp,1,GL_FALSE,glm::value_ptr(mvp));

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
        }
        logic();
        // cout << "Pointer to window:" << window<<endl;
        render(window);
    }
}

int main(int argc, char* argv[]){
    float* pos = (float*)malloc(3*sizeof(float));
    pos[0]=0.0;pos[1]=0.0;pos[2]=0.0;
    float* vel = (float*)malloc(3*sizeof(float));
    vel[0]=0.0;vel[1]=0.0;vel[2]=0.0;
    float* force = (float*)malloc(3*sizeof(float));
    force[0]=0.0;force[1]=0.0;force[2]=0.0;
    p.add_particle(pos,vel,force,25.0,2.0);
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
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    main_loop(window);
    free_resources();
    return EXIT_SUCCESS;
}

