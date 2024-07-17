#include <vector>
#include <cmath>
#include <iostream>
#include <limits.h>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const int width  = 200;
const int height = 200;

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

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
    //sort by y vertic
    if(t0.y>t1.y) {std::swap(t0, t1);}
    if(t1.y>t2.y) {std::swap(t1, t2);}
    if(t0.y>t1.y) {std::swap(t0, t1);}
    // line(t0, t1, image, green);
    // line(t1, t2, image, green);
    // line(t2, t0, image, red);

    //cut B boundary to 2 segments, so A boundary will cut into 2 segments too.

    // int totalHeight = t2.y - t0.y;
    // //botton segment
    // int segmentHeight = t1.y - t0.y;
    // for(int y=t0.y; y<=t1.y; ++y){
    //     float alpha = (float)(y-t0.y) / totalHeight;
    //     float beta = (float)(y-t0.y) / segmentHeight; // be careful with divisions by zero 
    //     Vec2i A = t0 + (t2-t0) * alpha;
    //     Vec2i B = t0 + (t1-t0) * beta;
    //     //found boundary of 2 side, now fill it.
    //     if(A.x > B.x) {std::swap(A,B);}
    //     for(int i=A.x; i<=B.x; ++i){
    //         image.set(i, y, color);         // attention, due to int casts t0.y+i != A.y
    //     }
    // }
    // //top segment
    // segmentHeight = t2.y - t1.y;
    // for(int y=t1.y; y<=t2.y; ++y){
    //     float alpha = (float)(y-t0.y) / totalHeight;
    //     float beta = (float)(y-t1.y) / segmentHeight;
    //     Vec2i A = t0 + (t2-t0) * alpha;
    //     Vec2i B = t1 + (t2-t1) * beta;
    //     //found boundary of 2 side, now fill it.
    //     if(A.x > B.x) {std::swap(A,B);}
    //     for(int i=A.x; i<=B.x; ++i){
    //         image.set(i, y, color);
    //     }
    // }

    //more improvement
    int totalHeight = t2.y - t0.y;
    for(int y=0; y<=totalHeight; ++y){
        bool isSecondHalf = y>t1.y-t0.y || t1.y==t0.y;
        int segmentHeight = isSecondHalf? t2.y-t1.y : t1.y-t0.y;
        float alpha = (float)(y) / totalHeight;
        float beta = isSecondHalf? (float)(y-(t1.y-t0.y)) / segmentHeight : (float)(y) / segmentHeight; // be careful, with above condition, it is no divide by zero now.
        Vec2i A = t0 + (t2-t0) * alpha;
        Vec2i B = isSecondHalf ? t1 + (t2-t1) * beta : t0 + (t1-t0) * beta; ;
        //found boundary of 2 side, now fill it.
        if(A.x > B.x) {std::swap(A,B);}
        for(int i=A.x; i<=B.x; ++i){
            image.set(i, t0.y+y, color);        
        }
    }
}

void triangle_contour(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color){
    if(t0.y>t1.y) {std::swap(t0, t1);}
    if(t1.y>t2.y) {std::swap(t1, t2);}
    if(t0.y>t1.y) {std::swap(t0, t1);}
    line(t0, t1, image, green);
    line(t1, t2, image, green);
    line(t2, t0, image, red);
}
//author
Vec3f barycentric(Vec2i *pts, Vec2i P) { 
    Vec3f u = Vec3f(pts[2][0]-pts[0][0], pts[1][0]-pts[0][0], pts[0][0]-P[0])^Vec3f(pts[2][1]-pts[0][1], pts[1][1]-pts[0][1], pts[0][1]-P[1]);
    /* `pts` and `P` has integer value as coordinates
       so `abs(u[2])` < 1 means `u[2]` is 0, that means
       triangle is degenerate, in this case return something with negative coordinates */
    if (std::abs(u.z)<1) return Vec3f(-1,1,1);
    return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z); 
}

std::pair<Vec2i, Vec2i> foundBoundingBox(Vec2i t0, Vec2i t1, Vec2i t2){
    int xMax=INT_MIN, xMin=INT_MAX;
    int yMax=INT_MIN, yMin=INT_MAX;
    xMax = std::max(t0.x, t1.x);
    xMax = std::max(xMax, t2.x);
    xMin = std::min(t0.x, t1.x);
    xMin = std::min(xMin, t2.x);
    yMax = std::max(t0.y, t1.y);
    yMax = std::max(yMax, t2.y);
    yMin = std::min(t0.y, t1.y);
    yMin = std::min(yMin, t2.y);
    //std::cout << xMin << yMin << " " << xMax << yMax << std::endl;
    return {Vec2i(xMin, yMin), Vec2i(xMax, yMax)};
}

bool inside(int x, int y, Vec2i t0, Vec2i t1, Vec2i t2){
    //https://davidhsu666.com/archives/barycentric-coordinates/
    //三角形ABC內的點座標P可以表示為從A頂點加上A頂點分別到B頂點和C頂點的向量分量的和
    //又可以寫成(u,v,1) dot product (AB, AC, PA) = 0, 又可以拆成x,y分量 為(u,v,1) dot product (ABx, ACx, PAx) = 0, (u,v,1) dot product (ABy, ACy, PAy) = 0
    //可知(u,v,1)垂直(ABx, ACx, PAx), (ABy, ACy, PAy)得(u,v,1)為(ABx, ACx, PAx), (ABy, ACy, PAy)外積結果為(u',v',c')，但外積結果僅代表方向正確的向量，大小需要除以c'以符合假設(u,v,1)z軸為1
    Vec3f uv = Vec3f(t2.x-t0.x, t1.x-t0.x, t0.x-x)^Vec3f(t2.y-t0.y, t1.y-t0.y, t0.y-y);
    uv = Vec3f(uv.x/uv.z, uv.y/uv.z, uv.z/uv.z);
    if(uv.x+uv.y > 1 || uv.x < 0 || uv.y < 0) {return false;}
    return true;

    //same result
    // Vec2i t[3] = {Vec2i(t0.x, t0.y), Vec2i(t1.x, t1.y), Vec2i(t2.x, t2.y)};
    // Vec3f bc_screen = barycentric(t, Vec2i(x,y));
    // if(bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) {return false;}
    // return true;
}

void triangle_2(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color){
    std::pair<Vec2i, Vec2i> bbox = foundBoundingBox(t0, t1, t2);
    for(int x=bbox.first.x; x<=bbox.second.x; ++x){
        for(int y=bbox.first.y; y<=bbox.second.y; ++y){
            if(inside(x, y, t0, t1, t2)){
                image.set(x, y, color);
            }
        }
    }
}


int main(int argc, char** argv) {
    TGAImage image(width, height, TGAImage::RGB);

    Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80)};
    Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)};
    Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)};
    //for(int i=0; i<1000000; ++i){
    triangle_2(t0[0], t0[1], t0[2], image, red);
    triangle_2(t1[0], t1[1], t1[2], image, white);
    triangle_2(t2[0], t2[1], t2[2], image, green);
    //}
    // triangle_contour(t0[0], t0[1], t0[2], image, red);
    // triangle_contour(t1[0], t1[1], t1[2], image, white);
    // triangle_contour(t2[0], t2[1], t2[2], image, green);

    Model* model = new Model("obj/african_head.obj");
    int width_2 = 800;
    int height_2 = 800;
    TGAImage image2(width_2, height_2, TGAImage::RGB);
    Vec3f light_dir(-1,-1,-1);
    for (int i=0; i<model->nfaces(); i++) {
        //every face contain a triangle represented by 3 vertics.
        std::vector<int> face = model->face(i);
        Vec2i screen_coords[3];
        Vec3f world_coords[3];
        for (int j=0; j<3; j++) {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = Vec2i((v.x+1)*width_2/2., (v.y+1)*height_2/2.);
            world_coords[j]  = v;
        }
        //if n = (world_coords[1]-world_coords[0])^(world_coords[2]-world_coords[0]); will show back side of face.
        //because of normal of surface could be two direction determinated by which vector cross product another vector
        //here I guess the cross product is along with -z side because of light is (0, 0, -1) vector and author want intensity > 0 is bright side.
        Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
        //normalize so the intensity will be -1 to 1.
        n.normalize(); 
        float intensity = n*light_dir;
        //Back-face culling, if intensity is negative, it means back side of face. so can be ignored to draw.
        if (intensity>0) {
            triangle_2(screen_coords[0], screen_coords[1], screen_coords[2], image2, TGAColor(intensity*255, intensity*255, intensity*255, 255)); 
        } 
        //random fill color
        //triangle_2(screen_coords[0], screen_coords[1], screen_coords[2], image2, TGAColor(rand()%255, rand()%255, rand()%255, 255));
    }

    // image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    // image.write_tga_file("test.tga");
    image2.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image2.write_tga_file("face_fill_light_4.tga");
    return 0;
}