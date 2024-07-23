#ifndef __OUR_GL_H__
#define __OUR_GL_H__

#include "tgaimage.h"
#include "geometry.h"

extern Matrix View;
extern Matrix ViewPort;
extern Matrix Projection;

struct IShader {
    virtual ~IShader();
    virtual Matrix vertex(int iface, int nthvert) = 0;
    virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
};

Matrix createViewPort(int w, int h, int d);
Matrix createModelView(Vec3f eye, Vec3f center, Vec3f up);
Matrix creatProjection(float coeff); //coeff = -1/c;

Vec3f barycentric(Vec3f *pts, Vec3f P);
void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color);

void triangle_shader(Matrix *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer);
//void trinagle_texture_Flat_shading(Vec3f* pts, Vec2f *v_pts, float* zbuffer, TGAImage &image, TGAImage &t_image, float intensity);

#endif