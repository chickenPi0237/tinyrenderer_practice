#include <vector>
#include <iostream>
#include <limits>

#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h"

Model *model     = NULL;
float *shadow_buffer;
const int width  = 800;
const int height = 800;

//light_dir here is start from surface, previous leesson i use start from light point. btw light_dir = - light_dir_from_light_point
Vec3f light_dir(1,1,1);
Vec3f       eye(1,1,3);
Vec3f    center(0,0,0);
Vec3f        up(0,1,0);

struct DepthShader : public IShader {
    mat<3,3,float> varying_tri;

    DepthShader() : varying_tri(){}

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Viewport*Projection*ModelView*gl_Vertex;     // transform it to screen coordinates
        varying_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        Vec3f p = varying_tri*bar;
        color = TGAColor(255,255,255,255)*(p.z/depth);
        return false;
    }
};

struct PhongShader : public IShader {
    mat<2,3,float> varying_uv;  // triangle uv coordinates, written by the vertex shader, read by the fragment shader
    mat<3,3,float> varying_nrm; // normal per vertex to be interpolated by FS
    mat<4,3,float> varying_tri;
    mat<3,3,float> view_tri; //triagnle in view coordinates
    mat<3,3,float> ndc_tri;
    mat<4,4,float> uniform_m;
    mat<4,4,float> uniform_mti;
    

    virtual Vec4f vertex(int iface, int nthvert) {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        varying_nrm.set_col(nthvert, proj<3>(uniform_mti*embed<4>(model->normal(iface, nthvert), 0.f)));
        Vec4f gl_Vertex = Viewport*Projection*ModelView*embed<4>(model->vert(iface, nthvert));
        varying_tri.set_col(nthvert, gl_Vertex);
        // ndc_tri, varying_tri, view_tri, I have tried. cause slighty changed showing as suffix _2 _3 .tga
        //ndc_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        ndc_tri.set_col(nthvert, model->vert(iface, nthvert));
        view_tri.set_col(nthvert, proj<3>(ModelView*embed<4>(model->vert(iface, nthvert))));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        Vec3f bn = (varying_nrm*bar).normalize();
        
        Vec2f uv = varying_uv*bar;
        //care this embed<4>(light_dir, 0.f) use embed<4>(light_dir, 1.f(default)) will be different. looks like light coming from another direction.
        Vec3f l = proj<3>(uniform_m*embed<4>(light_dir, 0.f)).normalize();

        mat<3,3,float> A; //to transform to tangent basis vector
        //suffix _4.
        // A[1] = proj<3>(varying_tri.col(2)/varying_tri.col(2)[3]-varying_tri.col(0)/varying_tri.col(0)[3]);
        // A[0] = proj<3>(varying_tri.col(1)/varying_tri.col(1)[3]-varying_tri.col(0)/varying_tri.col(0)[3]);
        //sufix _6
        // A[1] = proj<3>(varying_tri.col(2)-varying_tri.col(0));
        // A[0] = proj<3>(varying_tri.col(1)-varying_tri.col(0));
        // suffix _5
        // A[0] = view_tri.col(1)-view_tri.col(0);
        // A[1] = view_tri.col(2)-view_tri.col(0);
        // suffix _7 with ndc_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        // suffix _8 with ndc_tri.set_col(nthvert, model->vert(iface, nthvert));
        A[0] = ndc_tri.col(1)-ndc_tri.col(0);
        A[1] = ndc_tri.col(2)-ndc_tri.col(0);
        A[2] = bn;
        mat<3,3,float> A_it = A.invert();
        Vec3f i = A_it * Vec3f(varying_uv[0][1]-varying_uv[0][0], varying_uv[0][2]-varying_uv[0][0], 0);
        Vec3f j = A_it * Vec3f(varying_uv[1][1]-varying_uv[1][0], varying_uv[1][2]-varying_uv[1][0], 0);
        mat<3,3,float> B;
        B.set_col(0, i.normalize());
        B.set_col(1, j.normalize());
        B.set_col(2, bn);
        // B[0] = i.normalize();
        // B[1] = j.normalize();
        // B[2] = bn;
        //B = B.transpose();
        Vec3f n = (B*model->normal(uv)).normalize();

        float diff = std::min(1.f, std::max(0.f, n*l));
        color = model->diffuse(uv)*diff;
        //color = TGAColor(255,255,255,255)*diff;
        return false;
    }
};

struct GouraudShader : public IShader {
    //Vec3f varying_intensity; // written by vertex shader, read by fragment shader
    mat<2,3,float> varying_uv;
    mat<4,4,float> uniform_m;
    mat<4,4,float> uniform_mti;

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Viewport*Projection*ModelView*gl_Vertex;     // transform it to screen coordinates
        //after affine mapping, normal vector should mapped by inverse(transpose(map)).
        //varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity 
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        //float intensity = varying_intensity*bar;   // interpolate intensity for the current pixel
        Vec2f uv = varying_uv*bar;
        //normal(Vec2f) read from normal map, which stored normal vector's xyz as rgb.
        Vec3f n = proj<3>(uniform_mti*embed<4>(model->normal(uv))).normalize();
        //why light_dir need transform, isn't light_dir stationary? if we don't transform light_dir, it would be a light come from the screen coordinate. 
        //check comparison of no_light_transform/light_transform pictures. especially *l-100_e300
        //to my code, uniform_mti*light_dir is work, author use uniform_m.
        Vec3f i = proj<3>(uniform_mti*embed<4>(light_dir)).normalize();

        //author code phong-reflect model
        // Vec3f r = (n*(n*i*2.f) - i).normalize();   // reflected light
        // float spec = pow(std::max(r.z, 0.0f), model->specular(uv));
        // //std::cout << model->specular(uv) << " ";
        // float diff = std::max(0.f, n*i);
        // TGAColor c = model->diffuse(uv);
        // color = c;
        // for (int i=0; i<3; i++) color[i] = std::min<float>(5 + c[i]*(0.8*diff + .6*spec), 255);

        //blinn-phong reflection
        Vec3f eye_transformed = proj<3>(uniform_mti*embed<4>(eye)).normalize();
        int glossy_level = 50;
        Vec3f h = (i+eye_transformed).normalize();
        float spec = std::pow(std::max(n*h, 0.0f), glossy_level);
        float diffuse = std::max(0.f, n*i);
        float abiment = 5/255;
        //color = TGAColor(255, 255, 255)*intensity; // well duh
        TGAColor c = model->diffuse(uv);
        color = c;
        //normally sum of scalar coefficient must be equal to 1
        for(int i=0; i<3; ++i) { color[i]=std::min<float>( c[i]*(0.1*abiment + 0.6*diffuse + 0.3*spec) ,255); }
        return false;                              // no, we do not discard this pixel
    }
};

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    light_dir.normalize();
    //TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);
    float* zbuffer_f = new float[width*height];
    shadow_buffer = new float[width*height];
    for(int i=0; i<width*height; ++i){
        zbuffer_f[i] = shadow_buffer[i] = -std::numeric_limits<float>::max();
    }
    //first pass, find the shadow map.
    lookat(light_dir, center, up);
    viewport(0, 0, width, height);
    projection(-1.f/(light_dir-center).norm());
    TGAImage depth(width, height, TGAImage::RGB);
    DepthShader depthshader;
    for (int i=0; i<model->nfaces(); i++) {
        Vec4f screen_coords[3];
        for (int j=0; j<3; j++) {
            screen_coords[j] = depthshader.vertex(i, j);
        }
        
        triangle_my(screen_coords, depthshader, depth, shadow_buffer);
    }
    depth.flip_vertically();
    depth.write_tga_file("depth.tga");

    //second pass, render all the thing.
    // TGAImage image  (width, height, TGAImage::RGB);
    
    // lookat(eye, center, up);
    // viewport(0, 0, width, height);
    // projection(-1.f/(eye-center).norm());
    // GouraudShader shader;
    // shader.uniform_m = Projection*ModelView;
    // shader.uniform_mti = (Projection*ModelView).invert_transpose();
    // for (int i=0; i<model->nfaces(); i++) {
    //     Vec4f screen_coords[3];
    //     for (int j=0; j<3; j++) {
    //         screen_coords[j] = shader.vertex(i, j);
    //     }
    //     //triangle(screen_coords, shader, image, zbuffer);
    //     triangle_my(screen_coords, shader, image, zbuffer_f);
    // }

    // image.  flip_vertically(); // to place the origin in the bottom left corner of the image
    // zbuffer.flip_vertically();
    // image.  write_tga_file("output_my_phongshader_tangent_t.tga");
    // zbuffer.write_tga_file("zbuffer_my_phongshader_tangent_t.tga");

    // { // dump z-buffer (debugging purposes only)
    //     TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
    //     for (int i=0; i<width; i++) {
    //         for (int j=0; j<height; j++) {
    //             //unsigned char color  = ((zbuffer[i+j*width]+1)/2)*255;
    //             zbimage.set(i, j, TGAColor(zbuffer_f[i+j*width]));
    //         }
    //     }
    //     zbimage.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    //     zbimage.write_tga_file("zbimage_phongshader.tga");
    // }
    delete [] zbuffer_f;
    delete [] shadow_buffer;
    delete model;
    return 0;
}
