#include<stdio.h>
#include<stdlib.h>
#include<algorithm>
#include<math.h>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color){
    //if dx < dy means more points in y direction is better.
    bool steep = false;
    if(std::abs(x1-x0) < std::abs(y1-y0)){
        //after swap, x value actually is y value.
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    //if end point smaller then start point
    if(x1 < x0){
        std::swap(x1, x0);
        std::swap(y1, y0);
    }
    //further optimize loop
    // int dx = x1-x0;
    // int dy = y1-y0;
    // int derror = dy*2;
    // int error = 0;
    // int y = y0;
    // for (int x=x0; x<=x1; x++) {
    //     if(steep){
    //         image.set(y, x, color);
    //     } else{
    //         image.set(x, y, color);
    //     }
    //     error += derror;
    //     if(error>=dx){
    //         y += (y1>y0)?1:-1;
    //         error -= dx*2;
    //     }
    // }
    //optimize loop
    int dx = x1-x0;
    int dy = y1-y0;
    float derror = std::abs(dy/float(dx));
    float error = 0;
    int y = y0;
    for (int x=x0; x<=x1; x++) {
        if(steep){
            image.set(y, x, color);
        } else{
            image.set(x, y, color);
        }
        error += derror;
        if(error>=.5){
            y += (y1>y0)?1:-1;
            error -= 1;
        }
    }
    // for (int x=x0; x<=x1; x++) { 
    //     float t = (x-x0)/(float)(x1-x0);
    //     //int y = y0*(1.-t) + y1*t; use float line won't jiggle.
    //     float y = y0*(1.-t) + y1*t; 
    //     if(steep){
    //         image.set(y, x, color);
    //     } else{
    //         image.set(x, y, color);
    //     }
    // }
    //more step ver
    // for (float x=x0; x<=x1; x+=0.1) { 
    //     float t = (x-x0)/(float)(x1-x0); 
    //     int y = y0*(1.-t) + y1*t; 
    //     if(steep){
    //         image.set(y, x, color);
    //     } else{
    //         image.set(x, y, color);
    //     }
    // } 
}

int main(int argc, char** argv) {
    // TGAImage image(100, 100, TGAImage::RGB);
    // //for(int i=0; i<1000000; ++i){
    //     line(13, 20, 80, 40, image, white); 
    //     line(20, 13, 40, 80, image, red); 
    //     line(80, 40, 13, 20, image, red);
    //     line(60, 40, 60, 20, image, TGAColor(0, 255,   0,   255));
    //     line(1, 90, 100, 90, image, TGAColor(0, 255,   0,   255));
    // //}
    // image.flip_vertically(); // if want left-bottom corner is origion.
    // image.write_tga_file("output_optimized_further.tga");

    //draw wiremesh
    Model* model = new Model("obj/african_head.obj");
    int width = 800;
    int height = 800;
    TGAImage image2(width, height, TGAImage::RGB);
    for (int i=0; i<model->nfaces(); i++) {
        //every face contain a triangle represented by 3 vertics.
        std::vector<int> face = model->face(i);
        for (int j=0; j<3; j++) {
            Vec3f v0 = model->vert(face[j]);
            Vec3f v1 = model->vert(face[(j+1)%3]);
            //scaling
            int x0 = (v0.x+1.)*width/2.;
            int y0 = (v0.y+1.)*height/2.;
            int x1 = (v1.x+1.)*width/2.;
            int y1 = (v1.y+1.)*height/2.;
            line(x0, y0, x1, y1, image2, white);
        }
    }
    image2.flip_vertically();
    image2.write_tga_file("wiremesh.tga");
    delete model;
    return 0;
}