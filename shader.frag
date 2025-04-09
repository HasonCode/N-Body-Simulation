#version 120
varying vec3 f_color;
varying vec3 normal_vector;
varying vec3 f_pos;
uniform float fade;
uniform vec3 lightpos;
void main(void){
    float ambStrength = 0.4;
    vec3 ambient = ambStrength * vec3(1.0,1.0,1.0);
    vec3 lightdir = normalize(lightpos-f_pos);
    float diff = max(dot(normal_vector,lightdir),0.0);
    diff = diff;
    vec3 res = (ambient+diff) * f_color;
    gl_FragColor = vec4(res,fade);

}