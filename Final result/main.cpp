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
TGAImage SSAO_frame(width, height, TGAImage::RGB);

//light_dir here is start from surface, previous leesson i use start from light point. btw light_dir = - light_dir_from_light_point
Vec3f light_dir(1,1,1);
Vec3f       eye(0,0,2);
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
        gl_Vertex = Projection*ModelView*gl_Vertex;     // not transform to screen coordinates yet.
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
        float ambient_SSAO = SSAO_frame.get(ssao_p[0], ssao_p[1])[0] * (25/255.f);

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

struct faceContourShader : public IShader{
    mat<4,3,float> varying_tri;
    virtual Vec4f vertex(int iface, int nthvert){
        Vec4f tri = Viewport*Projection*ModelView*embed<4>(model->vert(iface, nthvert));
        tri = tri/tri[3];
        varying_tri.set_col(nthvert, tri);
        return tri;
    }
    virtual bool fragment(Vec3f gl_FragCoord  ,Vec3f bar, TGAColor &color){
        //nothing
        return true;
    }
    virtual bool fragment(Vec3f bar, TGAColor &color){
        //try to use barycentric coordinated to check if on edge of triangle, but result dots on edge.
        if((1-bar[0]-bar[1])<0.05 || (1-bar[1]-bar[2])<0.05 || (1-bar[0]-bar[2])<0.05){
            color = TGAColor(255,255,255);
            return false;
        }
        return true;
    }
};

struct FlatShader : public IShader{
    mat<2,3,float> varying_uv;
    mat<3,3,float> varying_tri;
    mat<4,4,float> uniform_m;
    mat<4,4,float> uniform_mti;
    mat<4,3,float> dc_tri;
    virtual Vec4f vertex(int iface, int nthvert){
        Vec4f tri = embed<4>(model->vert(iface, nthvert));
        dc_tri.set_col(nthvert, Projection*ModelView*tri);
        tri = Viewport*Projection*ModelView*tri;
        varying_tri.set_col(nthvert, proj<3>(tri/tri[3]));
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        return tri;
    }
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color){return true;}
    virtual bool fragment(Vec3f bar, TGAColor &color){
        //flat_normal no need to transform by uniform_mti becasue varying_tri is transformed.
        Vec3f flat_normal = cross((varying_tri.col(1)-varying_tri.col(0)), (varying_tri.col(2)-varying_tri.col(0)));
        flat_normal = flat_normal.normalize();
        Vec2f uv = varying_uv*bar;
        //transform light_dir, so light's coordinates will be in world coordinates. otherwise will be in screen coordinates.
        Vec3f l = proj<3>(uniform_mti*embed<4>(light_dir, 0.f)).normalize();
        float intensity = flat_normal*l;
        TGAColor c = model->diffuse(uv);
        color = c*intensity;
        return false;
    }
};

struct GouraudShader_wo_ : public IShader {
    //Vec3f varying_intensity; // written by vertex shader, read by fragment shader
    mat<2,3,float> varying_uv;
    mat<4,4,float> uniform_m;
    mat<4,4,float> uniform_mti;
    mat<3,3,float> varying_tri;
    mat<3,3,float> varying_noraml;
    mat<4,3,float> dc_tri; //device tri

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        dc_tri.set_col(nthvert, Projection*ModelView*gl_Vertex);
        gl_Vertex = Viewport*Projection*ModelView*gl_Vertex;     // transform it to screen coordinates
        //after affine mapping, normal vector should mapped by inverse(transpose(map)).
        //varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity 
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        varying_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        varying_noraml.set_col(nthvert, model->normal(iface, nthvert));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        //get interploated coordinates of texture.
        Vec2f uv = varying_uv*bar;
        //get interploated coordinates of normal.
        Vec3f n = varying_noraml*bar;
        n = n.normalize();
        //why light_dir need transform, isn't light_dir stationary? if we don't transform light_dir, it would be a light come from the screen coordinate. 
        //check comparison of no_light_transform/light_transform pictures. especially *l-100_e300
        //to my code, uniform_mti*light_dir is work, author use uniform_m.
        Vec3f i = proj<3>(uniform_mti*embed<4>(light_dir, 0.f)).normalize();

        float diffuse = std::max(0.f, n*i);
        //color = TGAColor(255, 255, 255)*intensity; // well duh
        TGAColor c = model->diffuse(uv);
        color = c*diffuse;
        //normally sum of scalar coefficient must be equal to 1
        //for(int i=0; i<3; ++i) { color[i]=std::min<float>(ambient_SSAO + c[i]*shadow*( + 0.8*diffuse + 1.5*spec) + glow[i]*3.f ,255.f); }
        //color = glow + color;
        return false;                              // no, we do not discard this pixel
    }
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color){return true;}
};

struct GouraudShader_add_normalmap : public IShader {
    //Vec3f varying_intensity; // written by vertex shader, read by fragment shader
    mat<2,3,float> varying_uv;
    mat<4,4,float> uniform_m;
    mat<4,4,float> uniform_mti;
    mat<3,3,float> varying_tri;
    mat<4,3,float> dc_tri; //device tri

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        dc_tri.set_col(nthvert, Projection*ModelView*gl_Vertex);
        gl_Vertex = Viewport*Projection*ModelView*gl_Vertex;     // transform it to screen coordinates
        //after affine mapping, normal vector should mapped by inverse(transpose(map)).
        //varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity 
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        varying_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        //get interploated coordinates of texture.
        Vec2f uv = varying_uv*bar;
        //get interploated coordinates of normal.
        // Vec4f n_tmp = uniform_mti*embed<4>(model->normal(uv), 0.f);
        // Vec3f n = proj<3>(n_tmp/n_tmp[3]).normalize();
        //or ?
        Vec3f n = proj<3>(uniform_mti*embed<4>(model->normal(uv), 0.f)).normalize();

        //why light_dir need transform, isn't light_dir stationary? if we don't transform light_dir, it would be a light come from the screen coordinate. 
        //check comparison of no_light_transform/light_transform pictures. especially *l-100_e300
        //to my code, uniform_mti*light_dir is work, author use uniform_m.
        Vec3f i = proj<3>(uniform_mti*embed<4>(light_dir, 0.f)).normalize();

        float diffuse = std::max(0.f, n*i);
        //color = TGAColor(255, 255, 255)*intensity; // well duh
        TGAColor c = model->diffuse(uv);
        color = c*diffuse;
        //normally sum of scalar coefficient must be equal to 1
        //for(int i=0; i<3; ++i) { color[i]=std::min<float>(ambient_SSAO + c[i]*shadow*( + 0.8*diffuse + 1.5*spec) + glow[i]*3.f ,255.f); }
        //color = glow + color;
        return false;                              // no, we do not discard this pixel
    }
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color){return true;}
};

struct GouraudShader_add_normalmap_tangent : public IShader {
    //Vec3f varying_intensity; // written by vertex shader, read by fragment shader
    mat<2,3,float> varying_uv;
    mat<4,4,float> uniform_m;
    mat<4,4,float> uniform_mti;
    mat<3,3,float> ndc_tri;     // triangle in normalized device coordinates
    mat<3,3,float> varying_normal;
    mat<4,3,float> dc_tri;

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Projection*ModelView*gl_Vertex;
        dc_tri.set_col(nthvert, gl_Vertex);     // note here didn't transform to screen coordinates
        //after affine mapping, normal vector should mapped by inverse(transpose(map)).
        //varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity 
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        ndc_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        varying_normal.set_col(nthvert, proj<3>(uniform_mti*embed<4>(model->normal(iface, nthvert), 1.f)));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        //get interploated coordinates of texture.
        Vec2f uv = varying_uv*bar;
        //get interploated coordinates of normal.
        Vec3f n = (varying_normal*bar).normalize();

        mat<3,3,float> A;
        A[0] = ndc_tri.col(1)-ndc_tri.col(0);
        A[1] = ndc_tri.col(2)-ndc_tri.col(0);
        A[2] = n;
        A = A.invert();

        mat<3,3,float> Darboux;
        // Darboux[0] = (A*Vec3f(varying_uv[0][1]-varying_uv[0][0], varying_uv[0][2]-varying_uv[0][0], 0)).normalize();
        // Darboux[1] = (A*Vec3f(varying_uv[1][1]-varying_uv[1][0], varying_uv[1][2]-varying_uv[1][0], 0)).normalize();
        // Darboux[2] = n;
        // Darboux = Darboux.transpose();
        Darboux.set_col(0, (A*Vec3f(varying_uv[0][1]-varying_uv[0][0], varying_uv[0][2]-varying_uv[0][0], 0)).normalize());
        Darboux.set_col(1, (A*Vec3f(varying_uv[1][1]-varying_uv[1][0], varying_uv[1][2]-varying_uv[1][0], 0)).normalize());
        Darboux.set_col(2, n);

        n = (Darboux*model->normal(uv)).normalize();
        
        //why light_dir need transform, isn't light_dir stationary? if we don't transform light_dir, it would be a light come from the screen coordinate. 
        //check comparison of no_light_transform/light_transform pictures. especially *l-100_e300
        //to my code, uniform_mti*light_dir is work, author use uniform_m.
        Vec3f i = proj<3>(uniform_mti*embed<4>(light_dir, 0.f)).normalize();

        float diffuse = std::max(0.f, n*i);
        //color = TGAColor(255, 255, 255)*intensity; // well duh
        TGAColor c = model->diffuse(uv);
        color = c*diffuse;
        //normally sum of scalar coefficient must be equal to 1
        //for(int i=0; i<3; ++i) { color[i]=std::min<float>(ambient_SSAO + c[i]*shadow*( + 0.8*diffuse + 1.5*spec) + glow[i]*3.f ,255.f); }
        //color = glow + color;
        return false;                              // no, we do not discard this pixel
    }
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color){return true;}
};

//normalmap_tangent + specular
struct GouraudShader_add_spec : public IShader {
    //Vec3f varying_intensity; // written by vertex shader, read by fragment shader
    mat<2,3,float> varying_uv;
    mat<4,4,float> uniform_m;
    mat<4,4,float> uniform_mti;
    mat<3,3,float> ndc_tri;     // triangle in normalized device coordinates
    mat<3,3,float> varying_normal;
    mat<4,3,float> dc_tri;

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Projection*ModelView*gl_Vertex;
        dc_tri.set_col(nthvert, gl_Vertex);     // note here didn't transform to screen coordinates
        //after affine mapping, normal vector should mapped by inverse(transpose(map)).
        //varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity 
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        ndc_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        varying_normal.set_col(nthvert, proj<3>(uniform_mti*embed<4>(model->normal(iface, nthvert), 1.f)));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        //get interploated coordinates of texture.
        Vec2f uv = varying_uv*bar;
        //get interploated coordinates of normal.
        Vec3f n = (varying_normal*bar).normalize();

        mat<3,3,float> A;
        A[0] = ndc_tri.col(1)-ndc_tri.col(0);
        A[1] = ndc_tri.col(2)-ndc_tri.col(0);
        A[2] = n;
        A = A.invert();

        mat<3,3,float> Darboux;
        Darboux.set_col(0, (A*Vec3f(varying_uv[0][1]-varying_uv[0][0], varying_uv[0][2]-varying_uv[0][0], 0)).normalize());
        Darboux.set_col(1, (A*Vec3f(varying_uv[1][1]-varying_uv[1][0], varying_uv[1][2]-varying_uv[1][0], 0)).normalize());
        Darboux.set_col(2, n);

        n = (Darboux*model->normal(uv)).normalize();
        
        //why light_dir need transform, isn't light_dir stationary? if we don't transform light_dir, it would be a light come from the screen coordinate. 
        //check comparison of no_light_transform/light_transform pictures. especially *l-100_e300
        //to my code, uniform_mti*light_dir is work, author use uniform_m.
        Vec3f i = proj<3>(uniform_mti*embed<4>(light_dir, 0.f)).normalize();

        float diffuse = std::max(0.f, n*i);
        //float glossy_level = model->specular(uv);
        // Vec3f r = (n*(n*i*2.f) - i).normalize();;
        // float spec = pow(std::max<float>(r*Vec3f(0,0,1), 0.f), 20+model->specular(uv)); //after Projection and ModelView, camera is lying on z-axis now, so eye is simply (0,0,1)
        Vec3f eye_transformed = proj<3>(uniform_mti*embed<4>(eye)).normalize();
        //std::cout << eye_transformed << std::endl;
        int glossy_level = 50;
        Vec3f h = (i+eye_transformed).normalize(); //same as (i+Vec3f(0,0,1)).normalize();
        float spec = std::pow(std::max(n*h, 0.0f), glossy_level); 
        
        //color = TGAColor(255, 255, 255)*intensity; // well duh
        TGAColor c = model->diffuse(uv);
        color = c;
        //normally sum of scalar coefficient must be equal to 1
        for(int i=0; i<3; ++i) { color[i]=std::min<float>(5 + c[i]*(0.8*diffuse + 0.8*spec),255.f); }
        return false;                              // no, we do not discard this pixel
    }
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color){return true;}
};
//normalmap_tangent + specular + shadow
struct GouraudShader_add_shadow : public IShader {
    //Vec3f varying_intensity; // written by vertex shader, read by fragment shader
    mat<2,3,float> varying_uv;
    mat<4,4,float> uniform_m;
    mat<4,4,float> uniform_mti;
    mat<3,3,float> ndc_tri;     // triangle in normalized device coordinates
    mat<3,3,float> varying_normal;
    mat<4,3,float> dc_tri;
    mat<4,4,float> uniform_shadow;

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Projection*ModelView*gl_Vertex;
        dc_tri.set_col(nthvert, gl_Vertex);     // note here didn't transform to screen coordinates
        //after affine mapping, normal vector should mapped by inverse(transpose(map)).
        //varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity 
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        ndc_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        varying_normal.set_col(nthvert, proj<3>(uniform_mti*embed<4>(model->normal(iface, nthvert), 1.f)));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        //add shadow
        //transform tri vertex to shadow buffer coordinates
        Vec4f sp_b = (uniform_shadow*embed<4>(Viewport*dc_tri*bar));
        sp_b = sp_b/sp_b[3];
        //43.34 is just magic number to avoid z-fighting.
        float shadow = 0.3+0.7*(sp_b[2]+43.34>shadow_buffer[int(sp_b[0])+int(sp_b[1])*width]);
        
        //get interploated coordinates of texture.
        Vec2f uv = varying_uv*bar;
        //get interploated coordinates of normal.
        Vec3f n = (varying_normal*bar).normalize();

        mat<3,3,float> A;
        A[0] = ndc_tri.col(1)-ndc_tri.col(0);
        A[1] = ndc_tri.col(2)-ndc_tri.col(0);
        A[2] = n;
        A = A.invert();

        mat<3,3,float> Darboux;
        Darboux.set_col(0, (A*Vec3f(varying_uv[0][1]-varying_uv[0][0], varying_uv[0][2]-varying_uv[0][0], 0)).normalize());
        Darboux.set_col(1, (A*Vec3f(varying_uv[1][1]-varying_uv[1][0], varying_uv[1][2]-varying_uv[1][0], 0)).normalize());
        Darboux.set_col(2, n);

        n = (Darboux*model->normal(uv)).normalize();
        
        //why light_dir need transform, isn't light_dir stationary? if we don't transform light_dir, it would be a light come from the screen coordinate. 
        //check comparison of no_light_transform/light_transform pictures. especially *l-100_e300
        //to my code, uniform_mti*light_dir is work, author use uniform_m.
        Vec3f i = proj<3>(uniform_mti*embed<4>(light_dir, 0.f)).normalize();

        float diffuse = std::max(0.f, n*i);
        //float glossy_level = model->specular(uv);
        // Vec3f r = (n*(n*i*2.f) - i).normalize();;
        // float spec = pow(std::max<float>(r*Vec3f(0,0,1), 0.f), 20+model->specular(uv)); //after Projection and ModelView, camera is lying on z-axis now, so eye is simply (0,0,1)
        Vec3f eye_transformed = proj<3>(uniform_mti*embed<4>(eye)).normalize();
        //std::cout << eye_transformed << std::endl;
        int glossy_level = 50;
        Vec3f h = (i+eye_transformed).normalize(); //same as (i+Vec3f(0,0,1)).normalize();
        float spec = std::pow(std::max(n*h, 0.0f), glossy_level); 
        
        //color = TGAColor(255, 255, 255)*intensity; // well duh
        TGAColor c = model->diffuse(uv);
        color = c;
        //normally sum of scalar coefficient must be equal to 1
        for(int i=0; i<3; ++i) { color[i]=std::min<float>(5 + c[i]*shadow*(0.8*diffuse + 0.8*spec),255.f); }
        return false;                              // no, we do not discard this pixel
    }
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color){return true;}
};
//normalmap_tangent + specular + shadow + SSAO
struct GouraudShader_add_SSAO : public IShader {
    //Vec3f varying_intensity; // written by vertex shader, read by fragment shader
    mat<2,3,float> varying_uv;
    mat<4,4,float> uniform_m;
    mat<4,4,float> uniform_mti;
    mat<3,3,float> ndc_tri;     // triangle in normalized device coordinates
    mat<3,3,float> varying_normal;
    mat<4,3,float> dc_tri;
    mat<4,4,float> uniform_shadow;

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Projection*ModelView*gl_Vertex;
        dc_tri.set_col(nthvert, gl_Vertex);     // note here didn't transform to screen coordinates
        //after affine mapping, normal vector should mapped by inverse(transpose(map)).
        //varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity 
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        ndc_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        varying_normal.set_col(nthvert, proj<3>(uniform_mti*embed<4>(model->normal(iface, nthvert), 1.f)));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        //add shadow
        //transform tri vertex to shadow buffer coordinates,
        //notice uniform_shadow is shadow_m*(Viewport*Projection*ModelView).invert() which invert the transform of dc_tri then transform to shadow
        Vec4f sp_b = (uniform_shadow*embed<4>(Viewport*dc_tri*bar));
        sp_b = sp_b/sp_b[3];
        //43.34 is just magic number to avoid z-fighting.
        float shadow = 0.3+0.7*(sp_b[2]+43.34>shadow_buffer[int(sp_b[0])+int(sp_b[1])*width]);

        //ambiment from SSAO
        Vec2f ssao_p = proj<2>(Viewport*dc_tri*bar);
        float ambient_SSAO = SSAO_frame.get(ssao_p[0], ssao_p[1])[0] * (5/255.f);
        
        //get interploated coordinates of texture.
        Vec2f uv = varying_uv*bar;
        //get interploated coordinates of normal.
        Vec3f n = (varying_normal*bar).normalize();

        mat<3,3,float> A;
        A[0] = ndc_tri.col(1)-ndc_tri.col(0);
        A[1] = ndc_tri.col(2)-ndc_tri.col(0);
        A[2] = n;
        A = A.invert();

        mat<3,3,float> Darboux;
        Darboux.set_col(0, (A*Vec3f(varying_uv[0][1]-varying_uv[0][0], varying_uv[0][2]-varying_uv[0][0], 0)).normalize());
        Darboux.set_col(1, (A*Vec3f(varying_uv[1][1]-varying_uv[1][0], varying_uv[1][2]-varying_uv[1][0], 0)).normalize());
        Darboux.set_col(2, n);

        n = (Darboux*model->normal(uv)).normalize();
        
        //why light_dir need transform, isn't light_dir stationary? if we don't transform light_dir, it would be a light come from the screen coordinate. 
        //check comparison of no_light_transform/light_transform pictures. especially *l-100_e300
        //to my code, uniform_mti*light_dir is work, author use uniform_m.
        Vec3f i = proj<3>(uniform_mti*embed<4>(light_dir, 0.f)).normalize();

        float diffuse = std::max(0.f, n*i);
        //float glossy_level = model->specular(uv);
        // Vec3f r = (n*(n*i*2.f) - i).normalize();;
        // float spec = pow(std::max<float>(r*Vec3f(0,0,1), 0.f), 20+model->specular(uv)); //after Projection and ModelView, camera is lying on z-axis now, so eye is simply (0,0,1)
        Vec3f eye_transformed = proj<3>(uniform_mti*embed<4>(eye)).normalize();
        //std::cout << eye_transformed << std::endl;
        int glossy_level = 50;
        Vec3f h = (i+eye_transformed).normalize(); //same as (i+Vec3f(0,0,1)).normalize();
        float spec = std::pow(std::max(n*h, 0.0f), glossy_level); 
        
        //color = TGAColor(255, 255, 255)*intensity; // well duh
        TGAColor c = model->diffuse(uv);
        color = c;
        //normally sum of scalar coefficient must be equal to 1
        for(int i=0; i<3; ++i) { color[i]=std::min<float>(ambient_SSAO + c[i]*shadow*(0.8*diffuse + 0.8*spec),255.f); }
        return false;                              // no, we do not discard this pixel
    }
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color){return true;}
};
//normalmap_tangent + specular + shadow + SSAO + glow
struct GouraudShader_add_glow : public IShader {
    //Vec3f varying_intensity; // written by vertex shader, read by fragment shader
    mat<2,3,float> varying_uv;
    mat<4,4,float> uniform_m;
    mat<4,4,float> uniform_mti;
    mat<3,3,float> ndc_tri;     // triangle in normalized device coordinates
    mat<3,3,float> varying_normal;
    mat<4,3,float> dc_tri;
    mat<4,4,float> uniform_shadow;

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Projection*ModelView*gl_Vertex;
        dc_tri.set_col(nthvert, gl_Vertex);     // note here didn't transform to screen coordinates
        //after affine mapping, normal vector should mapped by inverse(transpose(map)).
        //varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity 
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        ndc_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        varying_normal.set_col(nthvert, proj<3>(uniform_mti*embed<4>(model->normal(iface, nthvert), 1.f)));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        //add shadow
        //transform tri vertex to shadow buffer coordinates,
        //notice uniform_shadow is shadow_m*(Viewport*Projection*ModelView).invert() which invert the transform of dc_tri then transform to shadow
        Vec4f sp_b = (uniform_shadow*embed<4>(Viewport*dc_tri*bar));
        sp_b = sp_b/sp_b[3];
        //43.34 is just magic number to avoid z-fighting.
        float shadow = 0.3+0.7*(sp_b[2]+43.34>shadow_buffer[int(sp_b[0])+int(sp_b[1])*width]);

        //ambiment from SSAO
        Vec2f ssao_p = proj<2>(Viewport*dc_tri*bar);
        float ambient_SSAO = SSAO_frame.get(ssao_p[0], ssao_p[1])[0] * (5/255.f);
        
        //get interploated coordinates of texture.
        Vec2f uv = varying_uv*bar;
        //get interploated coordinates of normal.
        Vec3f n = (varying_normal*bar).normalize();

        mat<3,3,float> A;
        A[0] = ndc_tri.col(1)-ndc_tri.col(0);
        A[1] = ndc_tri.col(2)-ndc_tri.col(0);
        A[2] = n;
        A = A.invert();

        mat<3,3,float> Darboux;
        Darboux.set_col(0, (A*Vec3f(varying_uv[0][1]-varying_uv[0][0], varying_uv[0][2]-varying_uv[0][0], 0)).normalize());
        Darboux.set_col(1, (A*Vec3f(varying_uv[1][1]-varying_uv[1][0], varying_uv[1][2]-varying_uv[1][0], 0)).normalize());
        Darboux.set_col(2, n);

        n = (Darboux*model->normal(uv)).normalize();
        
        //why light_dir need transform, isn't light_dir stationary? if we don't transform light_dir, it would be a light come from the screen coordinate. 
        //check comparison of no_light_transform/light_transform pictures. especially *l-100_e300
        //to my code, uniform_mti*light_dir is work, author use uniform_m.
        Vec3f i = proj<3>(uniform_mti*embed<4>(light_dir, 0.f)).normalize();

        float diffuse = std::max(0.f, n*i);
        //float glossy_level = model->specular(uv);
        // Vec3f r = (n*(n*i*2.f) - i).normalize();;
        // float spec = pow(std::max<float>(r*Vec3f(0,0,1), 0.f), 20+model->specular(uv)); //after Projection and ModelView, camera is lying on z-axis now, so eye is simply (0,0,1)
        Vec3f eye_transformed = proj<3>(uniform_mti*embed<4>(eye)).normalize();
        //std::cout << eye_transformed << std::endl;
        int glossy_level = 50;
        Vec3f h = (i+eye_transformed).normalize(); //same as (i+Vec3f(0,0,1)).normalize();
        float spec = std::pow(std::max(n*h, 0.0f), glossy_level); 
        
        //glow texture
        TGAColor glow = model->glow(uv);

        TGAColor c = model->diffuse(uv);
        color = c;
        //normally sum of scalar coefficient must be equal to 1
        for(int i=0; i<3; ++i) { color[i]=std::min<float>(ambient_SSAO + c[i]*shadow*(0.8*diffuse + 0.8*spec) + 1.6f*glow[i],255.f); }
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
        depth.write_tga_file("diablo3_pose_shadow.tga");

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
            triangle_my(zshader.varying_tri, zshader, SSAO_frame, abiment_buffer);
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
                SSAO_frame.set(x, y, TGAColor(total*255, total*255, total*255));
            }
        }
        


        //second pass, render all the thing.
        lookat(eye, center, up);
        viewport(0, 0, width, height);
        projection(-1.f/(eye-center).norm());
        //faceContourShader shader;
        //FlatShader shader;
        //GouraudShader_wo_ shader;
        //GouraudShader_add_normalmap shader;
        //GouraudShader_add_normalmap_tangent shader;
        //GouraudShader_add_spec shader;
        //GouraudShader_add_shadow shader;
        //GouraudShader_add_SSAO shader;
        GouraudShader_add_glow shader;
        //GouraudShader shader;
        shader.uniform_m = Projection*ModelView;
        shader.uniform_mti = (Projection*ModelView).invert_transpose();
        shader.uniform_shadow = shadow_m*((Viewport*Projection*ModelView).invert());
        for (int i=0; i<model->nfaces(); i++) {
            Vec4f screen_coords[3];
            for (int j=0; j<3; j++) {
                screen_coords[j] = shader.vertex(i, j);
            }
            //contour
            // for(int j=0; j<3; j++){
            //     line(proj<2>(screen_coords[j]), proj<2>(screen_coords[(j+1)%3]), image, TGAColor(255,255,255));
            // }
            //triangle(screen_coords, shader, image, zbuffer);
            triangle_my(shader.dc_tri, shader, image, zbuffer_f);
            
        }
        // SSAO_frame.flip_vertically();
        // SSAO_frame.write_tga_file("SSAO_diablo3_pose.tga");
        image.  flip_vertically(); // to place the origin in the bottom left corner of the image
        //zbuffer.flip_vertically();
        image.  write_tga_file("diablo3_pose_GouraudShader_add_glow_floor_perspective_correction.tga");
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
        //     zbimage.write_tga_file("debug_zbuffer.tga");
        // }
        delete [] zbuffer_f;
        delete [] shadow_buffer;
        delete [] abiment_buffer;
        delete model;
    }
    
    return 0;
}
