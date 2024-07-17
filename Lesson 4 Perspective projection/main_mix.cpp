#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const int width  = 800;
const int height = 800;
const int depth  = 255;

Model *model = NULL;
float *zbuffer = NULL;
Vec3f light_dir(0,0,-1);
Vec3f camera(0,0,10);

Vec3f m2v(Matrix m) {
    return Vec3f(m[0][0]/m[3][0], m[1][0]/m[3][0], m[2][0]/m[3][0]);
}

Matrix v2m(Vec3f v) {
    Matrix m(4, 1);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity(4);
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}

void triangle(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, TGAImage &image, float intensity, float *zbuffer) {
    if (t0.y==t1.y && t0.y==t2.y) return; // i dont care about degenerate triangles
    if (t0.y>t1.y) { std::swap(t0, t1); std::swap(uv0, uv1); }
    if (t0.y>t2.y) { std::swap(t0, t2); std::swap(uv0, uv2); }
    if (t1.y>t2.y) { std::swap(t1, t2); std::swap(uv1, uv2); }

    int total_height = t2.y-t0.y;
    for (int i=0; i<total_height; i++) {
        bool second_half = i>t1.y-t0.y || t1.y==t0.y;
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; // be careful: with above conditions no division by zero here
        Vec3i A   =               t0  + Vec3f(t2-t0  )*alpha;
        Vec3i B   = second_half ? t1  + Vec3f(t2-t1  )*beta : t0  + Vec3f(t1-t0  )*beta;
        Vec2i uvA =               uv0 +      (uv2-uv0)*alpha;
        Vec2i uvB = second_half ? uv1 +      (uv2-uv1)*beta : uv0 +      (uv1-uv0)*beta;
        if (A.x>B.x) { std::swap(A, B); std::swap(uvA, uvB); }
        for (int j=A.x; j<=B.x; j++) {
            float phi = B.x==A.x ? 1. : (float)(j-A.x)/(float)(B.x-A.x);
            Vec3i   P = Vec3f(A) + Vec3f(B-A)*phi;
            Vec2i uvP =     uvA +   (uvB-uvA)*phi;
            int idx = P.x+P.y*width;
            if (zbuffer[idx]<P.z) {
                zbuffer[idx] = P.z;
                TGAColor color = model->diffuse(uvP);
                image.set(P.x, P.y, TGAColor(color.r*intensity, color.g*intensity, color.b*intensity));
                //image.set(P.x, P.y, TGAColor(255*intensity, 255*intensity, 255*intensity));
            }
        }
    }
}


Vec3f barycentric(Vec3f *pts, Vec3f P) {
    Vec3f v[2];
    for(int i=0; i<2; ++i){
        v[i][0] = pts[2][i] - pts[0][i]; // A->C
        v[i][1] = pts[1][i] - pts[0][i]; // A->B
        v[i][2] = pts[0][i] - P[i];      // P->A
    }
    Vec3f u = v[0] ^ v[1];
    /* `pts` and `P` has integer value as coordinates
       so `abs(u[2])` < 1 means `u[2]` is 0, that means
       triangle is degenerate, in this case return something with negative coordinates */
    if (std::abs(u.z)>1e-2){
        // if u order is (A>C, A>B, P>A) the scalar of A>C is u.x now, according website, it should be z as return value.
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
        // if u order is (A>B, A>C, P>A) 
        //return Vec3f(1.f-(u.x+u.y)/u.z, u.x/u.z, u.y/u.z);   
    } 
    return Vec3f(-1,1,1);

    
}

//my_pratice
void trinagle_texture(Vec3f* pts, Vec2f *v_pts, float* zbuffer, TGAImage &image, const float intensity){
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax( -std::numeric_limits<float>::max(),  -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    // somehow bugged.
    // for (int i=0; i<3; i++) {
    //     for (int j=0; j<2; j++) {
    //         std::cout << pts[i][j] <<" ";
    //         std::cout <<std::min(bboxmin[j], pts[i][j]) << " " << std::max(bboxmax[j], pts[i][j]) << " check ";
    //         bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
    //         bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
    //     }
    //     std::cout << bboxmax.x << " " << bboxmax.y << " ";
    //     std::cout << bboxmin.x << " " << bboxmin.y;
    //     std::cout <<std::endl;
    // }
    
    bboxmax.x = std::max(bboxmax.x, pts[0].x);
    bboxmax.x = std::max(bboxmax.x, pts[1].x);
    bboxmax.x = std::max(bboxmax.x, pts[2].x);
    bboxmax.y = std::max(bboxmax.y, pts[0].y);
    bboxmax.y = std::max(bboxmax.y, pts[1].y);
    bboxmax.y = std::max(bboxmax.y, pts[2].y);

    bboxmin.x = std::min(bboxmin.x, pts[0].x);
    bboxmin.x = std::min(bboxmin.x, pts[1].x);
    bboxmin.x = std::min(bboxmin.x, pts[2].x);
    bboxmin.y = std::min(bboxmin.y, pts[0].y);
    bboxmin.y = std::min(bboxmin.y, pts[1].y);
    bboxmin.y = std::min(bboxmin.y, pts[2].y);
    // std::cout << pts[0] << pts[1] << pts[2];
    // std::cout << bboxmax[0] << " " << bboxmax[1] << " ";
    // std::cout << bboxmin[0] << " " << bboxmin[1] << " <> " << std::endl;
    // std::cout << bboxmax.x << " " << bboxmax.y << " ";
    // std::cout << bboxmin.x << " " << bboxmin.y << std::endl;
    // return ;
    Vec3f P;
    Vec2f T;
    for(P.x=bboxmin.x; P.x<=bboxmax.x; P.x++){
        for(P.y=bboxmin.y; P.y<=bboxmax.y; P.y++){
            // std::cout << P.x << P[0] << " " << P.y << P[1] << std::endl;
            // return ;
            //Vec2i t[3] = {Vec2f(t0.x, t0.y), Vec2i(t1.x, t1.y), Vec2i(t2.x, t2.y)};
            Vec3f bc_screen = barycentric(pts, P);
            if(bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) { continue; }
            P.z = 0.;
            T.x = 0.;
            T.y = 0.;
            for(int i=0; i<3; ++i){
                P.z += pts[i][2] * bc_screen[i];
                T.x += v_pts[i][0] * bc_screen[i];
                T.y += v_pts[i][1] * bc_screen[i];
            }
            if(zbuffer[int(P.x+P.y*width)]<P.z) {
                TGAColor t_color = model->diffuse(Vec2i(T));
                zbuffer[int(P.x+P.y*width)] = P.z;
                image.set(P.x, P.y, TGAColor(intensity*t_color.r, intensity*t_color.g, intensity*t_color.b, 255));
                //image.set(P[0], P[1], TGAColor(255*intensity, 255*intensity, 255*intensity));
            }
        }
    }
}

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    zbuffer = new float[width*height];
    for (int i=0; i<width*height; i++) {
        zbuffer[i] = std::numeric_limits<float>::min();
    }

    //light_dir.normalize();

    { // draw the model
        Matrix Projection = Matrix::identity(4);
        Matrix ViewPort   = viewport(width/8, height/8, width*3/4, height*3/4);
        Projection[3][2] = -1.f/camera.z;

        TGAImage image(width, height, TGAImage::RGB);
        for (int i=0; i<model->nfaces(); i++) {
            std::vector<int> face = model->face(i);
            Vec3f screen_coords[3];
            Vec3f world_coords[3];
            for (int j=0; j<3; j++) {
                Vec3f v = model->vert(face[j]);
                Vec3f tmp = m2v(ViewPort*Projection*v2m(v));
                //screen_coords[j] = Vec3f(int((v.x+1)*width/2.+.5), int((v.y+1)*height/2.+.5), v.z);
                //std::cout << screen_coords[j] << " " << tmp << std::endl;
                screen_coords[j] =  m2v(ViewPort*Projection*v2m(v));
                screen_coords[j].x += .5;
                screen_coords[j].y += .5;
                screen_coords[j].z += .5;
                world_coords[j]  = v;
            }
            Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
            n.normalize();
            float intensity = n*light_dir;
            if (intensity>0) {
                Vec2f uv[3];
                for (int k=0; k<3; k++) {
                    Vec2i tmp = model->uv(i, k);
                    uv[k] = Vec2f(tmp.x, tmp.y);
                    //std::cout << tmp << " uv" << uv[k] << uv[k].x << uv[k].y << std::endl;
                }
                //trinagle_texture(screen_coords, uv, zbuffer, image, intensity);
                triangle(Vec3i(screen_coords[0]), Vec3i(screen_coords[1]), Vec3i(screen_coords[2]), Vec2i(uv[0]), Vec2i(uv[1]), Vec2i(uv[2]), image, intensity, zbuffer);
            }
        }

        image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
        image.write_tga_file("output_2.tga");
    }

    { // dump z-buffer (debugging purposes only)
        TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                zbimage.set(i, j, TGAColor(zbuffer[i+j*width], 1));
            }
        }
        zbimage.flip_vertically(); // i want to have the origin at the left bottom corner of the image
        zbimage.write_tga_file("zbuffer_2.tga");
    }
    delete model;
    delete [] zbuffer;
    return 0;
}

