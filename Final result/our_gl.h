#ifndef __OUR_GL_H__
#define __OUR_GL_H__

#include "tgaimage.h"
#include "geometry.h"

extern Matrix ModelView;
extern Matrix Viewport;
extern Matrix Projection;
const float depth = 255.f;

void viewport(int x, int y, int w, int h);
void projection(float coeff=0.f); // coeff = -1/c
void lookat(Vec3f eye, Vec3f center, Vec3f up);

struct IShader {
    virtual ~IShader();
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color) = 0;
};

void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color);
void triangle(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer);
void triangle_my(Vec4f *pts, IShader &shader, TGAImage &image, float* zbuffer);
void triangle_my(mat<4,3,float> &clipc, IShader &shader, TGAImage &image, float* zbuffer);

#endif //__OUR_GL_H__
