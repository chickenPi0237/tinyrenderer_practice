#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>

#include "our_gl.h"

Matrix View;
Matrix ViewPort;
Matrix Projection;

IShader::~IShader(){}

Matrix createViewPort(int w, int h, int d){
    //Canonical cubic to screen.
    //so far z-axis value is not much important if z value in range (-1,1) or (0,255) set z to 0~255 is easily for debugging
    Matrix m = Matrix::identity(4);
    m[0][0] = w/2.;
    m[1][1] = h/2.;
    m[2][2] = d/2.;

    m[0][3] = w/2.;
    m[1][3] = h/2.;
    m[2][3] = d/2.;
    return m;
}

Matrix createModelView(Vec3f eye, Vec3f center, Vec3f up){
    Matrix transfer = Matrix::identity(4);
    Matrix rotate = Matrix::identity(4);
    Matrix view;
    //look at -z axis,
    Vec3f z = (eye-center).normalize();
    Vec3f x = (up^z).normalize();
    Vec3f y = (z^x).normalize();
    for(int i=0; i<3; ++i){
        rotate[0][i] = x[i];
        rotate[1][i] = y[i];
        rotate[2][i] = z[i];
        transfer[i][3] = -eye[i];
    }
    std::cout << transfer << rotate;
    view = rotate*transfer;
    //std::cout << view;
    return view;
}
Matrix creatProjection(float coeff){
    Matrix m = Matrix::identity();
    m[3][2] = coeff;
    return m;
}

void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(p0.x-p1.x)<std::abs(p0.y-p1.y)) {
        std::swap(p0.x, p0.y);
        std::swap(p1.x, p1.y);
        steep = true;
    }
    if (p0.x>p1.x) {
        std::swap(p0, p1);
    }

    for (int x=p0.x; x<=p1.x; x++) {
        float t = (x-p0.x)/(float)(p1.x-p0.x);
        int y = p0.y*(1.-t) + p1.y*t;
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
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


// void trinagle_texture_Flat_shading(Vec3f* pts, Vec2f *v_pts, float* zbuffer, TGAImage &image, TGAImage &t_image, float intensity){
//     //std::pair<Vec2i, Vec2i> bbox = foundBoundingBox(t0, t1, t2);
//     Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
//     Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
//     Vec2f clamp(image.get_width()-1, image.get_height()-1);
//     for (int i=0; i<3; i++) {
//         for (int j=0; j<2; j++) {
//             //std::cout << pts[i][j] << " ";
//             bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
//             bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
//         }
//         //std::cout << std::endl;
//     }
//     // std::cout << pts[0] << pts[1] << pts[2];
//     // std::cout << bboxmax[0] << " " << bboxmax[1] << " ";
//     // std::cout << bboxmin[0] << " " << bboxmin[1] << " <> " << std::endl;
//     // std::cout << bboxmax.x << " " << bboxmax.y << " ";
//     // std::cout << bboxmin.x << " " << bboxmin.y << std::endl;
//     //return ;
//     Vec3f P;
//     Vec2f T;
//     for(P.x=bboxmin.x; P.x<=bboxmax.x; P.x++){
//         for(P.y=bboxmin.y; P.y<=bboxmax.y; P.y++){
//             //Vec2i t[3] = {Vec2f(t0.x, t0.y), Vec2i(t1.x, t1.y), Vec2i(t2.x, t2.y)};
//             Vec3f bc_screen = barycentric(pts, P);
//             if(bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) { continue; }
//             P.z = 0.;
//             T.x = 0.;
//             T.y = 0.;
//             for(int i=0; i<3; ++i){
//                 P.z += pts[i][2] * bc_screen[i];
//                 T.x += v_pts[i][0] * bc_screen[i];
//                 T.y += v_pts[i][1] * bc_screen[i];
//             }
//             if(zbuffer[int(P.x+P.y*width_2)]<P.z) {
//                 //TGAColor t_color = t_image.get(int(T.x), int(T.y));
//                 TGAColor t_color = model->diffuse(Vec2i(T)); 
//                 zbuffer[int(P.x+P.y*width_2)] = P.z;
//                 image.set(P.x, P.y, TGAColor(intensity*t_color.r, intensity*t_color.g, intensity*t_color.b, 255));
//             }
//         }
//     }
// }

void triangle_shader(Matrix *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer){
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::min(bboxmin[j], pts[i][j]/pts[i][3]);
            bboxmax[j] = std::max(bboxmax[j], pts[i][j]/pts[i][3]);
        }
    }
    Vec3f P;
    TGAColor color;
    Vec3f vpts[3];
    for(int i=0; i<3; ++i){
        vpts[i] = Vec3f(pts[i]);
    }
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f c = barycentric(vpts, P);
            float z = pts[0][2]*c.x + pts[1][2]*c.y + pts[2][2]*c.z;
            float w = pts[0][3]*c.x + pts[1][3]*c.y + pts[2][3]*c.z;
            int frag_depth = std::max(0, std::min(255, int(z/w+.5)));
            if (c.x<0 || c.y<0 || c.z<0 || zbuffer.get(P.x, P.y)[0]>frag_depth) continue;
            bool discard = shader.fragment(c, color);
            if (!discard) {
                zbuffer.set(P.x, P.y, TGAColor(frag_depth));
                image.set(P.x, P.y, color);
            }
        }
    }

    // Vec3f P;
    // Vec2f T;
    // TGAColor color;
    // Vec3f vpts[3];
    // for(int i=0; i<3; ++i)
    //     vpts[i] = Vec3f(pts[i]);
    // for(P.x=bboxmin.x; P.x<=bboxmax.x; P.x++){
    //     for(P.y=bboxmin.y; P.y<=bboxmax.y; P.y++){
    //         //Vec2i t[3] = {Vec2f(t0.x, t0.y), Vec2i(t1.x, t1.y), Vec2i(t2.x, t2.y)};
    //         Vec3f bc_screen = barycentric(vpts, P);
    //         if(bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) { continue; }
    //         P.z = 0.;
    //         T.x = 0.;
    //         T.y = 0.;
    //         for(int i=0; i<3; ++i){
    //             P.z += pts[i][2] * bc_screen[i];
    //             T.x += v_pts[i][0] * bc_screen[i];
    //             T.y += v_pts[i][1] * bc_screen[i];
    //         }
    //         if(zbuffer[int(P.x+P.y*width_2)]<P.z) {
    //             //TGAColor t_color = t_image.get(int(T.x), int(T.y));
    //             TGAColor t_color = model->diffuse(Vec2i(T)); 
    //             zbuffer[int(P.x+P.y*width_2)] = P.z;
    //             if(v_intensity >= 0)
    //                 image.set(P.x, P.y, TGAColor(v_intensity*t_color.r, v_intensity*t_color.g, v_intensity*t_color.b, 255));
    //         }
    //     }
    // }

}

