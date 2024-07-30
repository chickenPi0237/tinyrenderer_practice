#include <cmath>
#include <limits>
#include <cstdlib>
#include "our_gl.h"

Matrix ModelView;
Matrix Viewport;
Matrix Projection;

IShader::~IShader() {}

void viewport(int x, int y, int w, int h) {
    Viewport = Matrix::identity();
    Viewport[0][3] = x+w/2.f;
    Viewport[1][3] = y+h/2.f;
    Viewport[2][3] = depth/2.f;
    Viewport[0][0] = w/2.f;
    Viewport[1][1] = h/2.f;
    Viewport[2][2] = depth/2.f;
}

void projection(float coeff) {
    Projection = Matrix::identity();
    Projection[3][2] = coeff;
}

void lookat(Vec3f eye, Vec3f center, Vec3f up) {
    // Vec3f z = (eye-center).normalize(); 
    // Vec3f x = cross(up,z).normalize();
    // Vec3f y = cross(z,x).normalize();
    // ModelView = Matrix::identity();
    // for (int i=0; i<3; i++) {
    //     ModelView[0][i] = x[i];
    //     ModelView[1][i] = y[i];
    //     ModelView[2][i] = z[i];
    //     ModelView[i][3] = -center[i];
    // }
    Matrix transfer = Matrix::identity();
    Matrix rotate = Matrix::identity();
    //look at -z axis,
    Vec3f z = (eye-center).normalize();
    Vec3f x = cross(up,z).normalize();
    Vec3f y = cross(z,x).normalize();
    for(int i=0; i<3; ++i){
        rotate[0][i] = x[i];
        rotate[1][i] = y[i];
        rotate[2][i] = z[i];
        transfer[i][3] = -eye[i];
    }
    ModelView = rotate*transfer;
}

Vec3f barycentric(Vec2f A, Vec2f B, Vec2f C, Vec2f P) {
    //std::cout << A << B << C << std::endl;
    Vec3f s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    Vec3f u = cross(s[0], s[1]);
    if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return Vec3f(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void triangle(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer) {
    //std::cout << *pts << std::endl;
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::min(bboxmin[j], pts[i][j]/pts[i][3]);
            bboxmax[j] = std::max(bboxmax[j], pts[i][j]/pts[i][3]);
        }
    }
    Vec2i P;
    TGAColor color;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f c = barycentric(proj<2>(pts[0]/pts[0][3]), proj<2>(pts[1]/pts[1][3]), proj<2>(pts[2]/pts[2][3]), proj<2>(P));
            if (c.x<0 || c.y<0 || c.z<0) continue;
            float z = pts[0][2]*c.x + pts[1][2]*c.y + pts[2][2]*c.z;
            float w = pts[0][3]*c.x + pts[1][3]*c.y + pts[2][3]*c.z;
            //std::cout << "a:" << (pts[0][2]*c.x + pts[1][2]*c.y + pts[2][2]*c.z)/(pts[0][3]*c.x + pts[1][3]*c.y + pts[2][3]*c.z) << std::endl;
            //std::cout << "m:" << (pts[0][2]/pts[0][3])*c.x + (pts[1][2]/pts[0][3])*c.y + (pts[2][2]/pts[0][3])*c.z << std::endl;
            // float z = (pts[0][2]/pts[0][3])*c.x + (pts[1][2]/pts[0][3])*c.y + (pts[2][2]/pts[0][3])*c.z;
            // float w = 1.;
            //int frag_depth = std::max(0, std::min(255, int(z/w+.5)));
            //int frag_depth = int(z/w+.5);
            int frag_depth = int(z/w+.5);
            //frag_depth = std::min(255, frag_depth);
            //frag_depth %= 255;
            if(frag_depth<0) { 
                frag_depth = std::min(0, std::max(-255, frag_depth)); 
                frag_depth += 255;
                //frag_depth = -frag_depth;
            } else{
                frag_depth = std::max(0, std::min(255, frag_depth));
            }
            //std::cout << frag_depth << std::endl;
            if (zbuffer.get(P.x, P.y)[0]>frag_depth) continue;
            bool discard = shader.fragment(c, color);
            if (!discard) {
                
                zbuffer.set(P.x, P.y, TGAColor(frag_depth));
                image.set(P.x, P.y, color);
            }
        }
    }
}

void triangle_my(Vec4f *pts, IShader &shader, TGAImage &image, float* zbuffer) {
    //std::cout << *pts << std::endl;
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::min(bboxmin[j], pts[i][j]/pts[i][3]);
            bboxmax[j] = std::max(bboxmax[j], pts[i][j]/pts[i][3]);
        }
    }
    Vec2i P;
    TGAColor color;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f c = barycentric(proj<2>(pts[0]/pts[0][3]), proj<2>(pts[1]/pts[1][3]), proj<2>(pts[2]/pts[2][3]), proj<2>(P));
            if (c.x<0 || c.y<0 || c.z<0) continue;
            float z = pts[0][2]*c.x + pts[1][2]*c.y + pts[2][2]*c.z;
            float w = pts[0][3]*c.x + pts[1][3]*c.y + pts[2][3]*c.z;
            
            //not this one, after MPV, interploation should be adjust.
            // float z = (pts[0][2]/pts[0][3])*c.x + (pts[1][2]/pts[0][3])*c.y + (pts[2][2]/pts[0][3])*c.z;
            // float w = 1.;
            //int frag_depth = std::max(0, std::min(255, int(z/w+.5)));

            //somehow with offset, casuing rounding problem.
            //int frag_depth = int(z/w+.5);
            float frag_depth = z/w;

            if (zbuffer[P.y*image.get_width()+P.x]>frag_depth) continue;
            bool discard = shader.fragment(c, color);
            if (!discard) {
                //zbuffer.set(P.x, P.y, TGAColor(frag_depth));
                zbuffer[P.y*image.get_width()+P.x] = frag_depth;
                image.set(P.x, P.y, color);
            }
        }
    }
}
