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
TGAImage frame(width, height, TGAImage::RGB);

//light_dir here is start from surface, previous leesson i use start from light point. btw light_dir = - light_dir_from_light_point
Vec3f light_dir(1,1,1);
Vec3f       eye(0,0,3);
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

    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color) {
        std::cout << "do nothing";
        return false;
    }
};

struct ZShader : public IShader {
    mat<4,3,float> varying_tri;

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Projection*ModelView*gl_Vertex;     // transform it to screen coordinates
        varying_tri.set_col(nthvert, gl_Vertex);
        return gl_Vertex;
    }
    virtual bool fragment(Vec3f bar, TGAColor &color){
        color = TGAColor(0,0,0);
        return false;
    }
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color) {
        color = TGAColor(0,0,0);
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
    mat<4,4,float> uniform_shadow;
    mat<3,3,float> varying_tri;

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Viewport*Projection*ModelView*gl_Vertex;     // transform it to screen coordinates
        //after affine mapping, normal vector should mapped by inverse(transpose(map)).
        //varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity 
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        varying_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        //float intensity = varying_intensity*bar;   // interpolate intensity for the current pixel
        Vec2f uv = varying_uv*bar;
        //add shadow
        //transform tri vertex to shadow buffer coordinates
        Vec4f sp_b = (uniform_shadow*embed<4>(varying_tri*bar));
        sp_b = sp_b/sp_b[3];
        //43.34 is just magic number to avoid z-fighting.
        float shadow = 0.3+0.7*(sp_b[2]+43.34>shadow_buffer[int(sp_b[0])+int(sp_b[1])*width]);

        //ambiment from SSAO
        Vec2f ssao_p = proj<2>(varying_tri*bar);
        float ambient_SSAO = frame.get(ssao_p[0], ssao_p[1])[0] * (25/255.f);

        //add glowmap
        TGAColor glow = model->glow(uv);

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
        for(int i=0; i<3; ++i) { color[i]=std::min<float>(ambient_SSAO + c[i]*shadow*( + 0.8*diffuse + 1.5*spec) + glow[i]*3.f ,255.f); }
        //color = glow + color;
        return false;                              // no, we do not discard this pixel
    }
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color){return true;}

};

float max_elevation_angle(float *zbuffer, Vec2f p, Vec2f dir) {
    float maxangle = 0;
    //go 1000 step in each direction.
    for (float t=0.; t<1000.; t+=1.) {
        Vec2f cur = p + dir*t;
        if (cur.x>=width || cur.y>=height || cur.x<0 || cur.y<0) return maxangle;

        float distance = (p-cur).norm();
        if (distance < 1.f) continue;
        float elevation = zbuffer[int(cur.x)+int(cur.y)*width]-zbuffer[int(p.x)+int(p.y)*width];
        maxangle = std::max(maxangle, atanf(elevation/distance));
    }
    //result will in range 0~90 degree.
    return maxangle;
}

int main(int argc, char** argv) {
    // if (2==argc) {
    //     model = new Model(argv[1]);
    // } else {
    //     model = new Model("obj/african_head.obj");
    // }

    light_dir.normalize();
    TGAImage image  (width, height, TGAImage::RGB);
    for(int i=1; i<argc; ++i){
        model = new Model(argv[i]);

        //TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);
        float* zbuffer_f = new float[width*height];
        shadow_buffer = new float[width*height];
        float* abiment_buffer = new float[width*height];
        for(int i=0; i<width*height; ++i){
            zbuffer_f[i] = shadow_buffer[i] = abiment_buffer[i] = -std::numeric_limits<float>::max();
        }
        // //first pass, find the shadow map. basically watch obj from light source so we can determine where we can't see aka where shadow should occure.
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
        depth.write_tga_file("depth_african_head_more.tga");

        mat<4,4,float> shadow_m = Viewport*Projection*ModelView;
        // ambient occulison pass
        // screen space ambient occlusion.
        // algorithm : 
        // emit rays outward from it in UV space.
        //     ex: up/down/right/left
        // record max slope at each direction ray go through
        // calculate how "open" or "occluded" is this point. if all max slope at every direction is 0, it is open.

        lookat(eye, center, up);
        viewport(0, 0, width, height);
        projection(-1.f/(eye-center).norm());
        ZShader zshader;
        for (int i=0; i<model->nfaces(); i++) {
            Vec4f screen_coords[3];
            for (int j=0; j<3; j++) {
                //didn't use screen_coords, but code can't run successfully without a variable to assign.
                screen_coords[j] = zshader.vertex(i, j);
            }
            //triangle(screen_coords, shader, image, zbuffer);
            //std::cout << zshader.varying_tri[2];
            triangle_my(zshader.varying_tri, zshader, frame, abiment_buffer);
        }
        for (int x=0; x<width; x++) {
            for (int y=0; y<height; y++) {
                if (abiment_buffer[x+y*width] < -1e5) continue;
                float total = 0;
                //use 8 direction to approximate result. should use solid angle.
                for (float a=0; a<M_PI*2-1e-4; a += M_PI/4) {
                    // M_PI/2 - maxnangle, because, if maxangle bigger, the point should be darker.
                    total += M_PI/2 - max_elevation_angle(abiment_buffer, Vec2f(x, y), Vec2f(cos(a), sin(a)));
                }
                // convert degree to 0~1.
                total /= (M_PI/2)*8;
                // pow 100 to result reasonable ambient light feeling. the higher the more parts is dark(lower value).
                //https://www.symbolab.com/solver/functions-calculator/f%5Cleft(x%5Cright)%3Dx%5E%7B100%7D?or=input to visulize function.
                total = pow(total, 100.f);
                frame.set(x, y, TGAColor(total*255, total*255, total*255));
            }
        }
        


        //second pass, render all the thing.
        lookat(eye, center, up);
        viewport(0, 0, width, height);
        projection(-1.f/(eye-center).norm());
        GouraudShader shader;
        shader.uniform_m = Projection*ModelView;
        shader.uniform_mti = (Projection*ModelView).invert_transpose();
        shader.uniform_shadow = shadow_m*((Viewport*Projection*ModelView).invert());
        for (int i=0; i<model->nfaces(); i++) {
            Vec4f screen_coords[3];
            for (int j=0; j<3; j++) {
                screen_coords[j] = shader.vertex(i, j);
            }
            //triangle(screen_coords, shader, image, zbuffer);
            triangle_my(screen_coords, shader, image, zbuffer_f);
        }
        frame.flip_vertically();
        frame.write_tga_file("SSAO_african_head_more.tga");
        image.  flip_vertically(); // to place the origin in the bottom left corner of the image
        //zbuffer.flip_vertically();
        image.  write_tga_file("output_my_african_head_shadow_SSAO_glow.tga");
        //zbuffer.write_tga_file("zbuffer_my_shadow.tga");

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
        delete [] abiment_buffer;
        delete model;
    }
    
    return 0;
}
