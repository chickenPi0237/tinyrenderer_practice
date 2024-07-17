#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const int width  = 200;
const int height = 200;
const int width_2 = 800;
const int height_2 = 800;

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
    // https://davidhsu666.com/archives/barycentric-coordinates/
    // 推導沒有錯，程式碼也沒有錯。
    // 問題來自程式碼的第一項是向量AC，由程式碼s[i][0] = C[i]-A[i]可知; 推導的則是向量AB。
    // 也就是說如果照著推導的回傳值應該是
    // Vec3f(1.f-(u.x+u.y)/u.z, u.x/u.z, u.y/u.z)
    // 但按照程式碼的則是
    // Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z)
    // u = (A>C, A>B, P>A)
    Vec3f u = Vec3f(pts[2][0]-pts[0][0], pts[1][0]-pts[0][0], pts[0][0]-P[0])^Vec3f(pts[2][1]-pts[0][1], pts[1][1]-pts[0][1], pts[0][1]-P[1]);
    // u = (A>B, A>C, P>A)
    //Vec3f u = Vec3f(pts[1][0]-pts[0][0], pts[2][0]-pts[0][0], pts[0][0]-P[0])^Vec3f(pts[1][1]-pts[0][1], pts[2][1]-pts[0][1], pts[0][1]-P[1]);

    /* `pts` and `P` has integer value as coordinates
       so `abs(u[2])` < 1 means `u[2]` is 0, that means
       triangle is degenerate, in this case return something with negative coordinates */
    if (std::abs(u.z)<1){
        return Vec3f(-1,1,1);
    } 

    // if u order is (A>C, A>B, P>A) the scalar of A>C is u.x now, according website, it should be z as return value.
    return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    // if u order is (A>B, A>C, P>A) 
    //return Vec3f(1.f-(u.x+u.y)/u.z, u.x/u.z, u.y/u.z); 
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
    // Vec3f uv = Vec3f(t2.x-t0.x, t1.x-t0.x, t0.x-x)^Vec3f(t2.y-t0.y, t1.y-t0.y, t0.y-y);
    // uv = Vec3f(uv.x/uv.z, uv.y/uv.z, uv.z/uv.z);
    // if(uv.x+uv.y > 1 || uv.x < 0 || uv.y < 0) {return false;}
    // return true;

    //same result
    Vec2i t[3] = {Vec2i(t0.x, t0.y), Vec2i(t1.x, t1.y), Vec2i(t2.x, t2.y)};
    Vec3f bc_screen = barycentric(t, Vec2i(x,y));
    if(bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) {return false;}
    return true;
}

void triangle_2(Vec3f* pts, float* zbuffer, TGAImage &image, TGAColor color){
    //std::pair<Vec2i, Vec2i> bbox = foundBoundingBox(t0, t1, t2);
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    Vec3f P;
    for(P.x=bboxmin.x; P.x<=bboxmax.x; P.x++){
        for(P.y=bboxmin.y; P.y<=bboxmax.y; P.y++){
            //Vec2i t[3] = {Vec2f(t0.x, t0.y), Vec2i(t1.x, t1.y), Vec2i(t2.x, t2.y)};
            Vec3f bc_screen = barycentric(pts, P);
            if(bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) { continue; }
            P.z = 0.;
            for(int i=0; i<3; ++i){
                P.z += pts[i][2] * bc_screen[i];
            }
            if(zbuffer[int(P.x+P.y*width_2)]<P.z) {
                zbuffer[int(P.x+P.y*width_2)] = P.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}

void trinagle_texture(Vec3f* pts, Vec2f *v_pts, float* zbuffer, TGAImage &image, TGAImage &t_image, const float intensity){
    //std::pair<Vec2i, Vec2i> bbox = foundBoundingBox(t0, t1, t2);
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            std::cout << pts[i][j] << " ";
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
        std::cout << std::endl;
    }
    std::cout << pts[0] << pts[1] << pts[2];
    std::cout << bboxmax[0] << " " << bboxmax[1] << " ";
    std::cout << bboxmin[0] << " " << bboxmin[1] << " <> " << std::endl;
    std::cout << bboxmax.x << " " << bboxmax.y << " ";
    std::cout << bboxmin.x << " " << bboxmin.y << std::endl;
    return ;
    Vec3f P;
    Vec2f T;
    for(P.x=bboxmin.x; P.x<=bboxmax.x; P.x++){
        for(P.y=bboxmin.y; P.y<=bboxmax.y; P.y++){
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
            if(zbuffer[int(P.x+P.y*width_2)]<P.z) {
                TGAColor t_color = t_image.get(int(T.x), int(T.y));
                //std::cout<<v_pts[0].x<<v_pts[1].y<<std::endl;
                //std::cout << t_image.get(T.x, T.y)[0] << " " << t_color[1] << " " << t_color[2] << std::endl;
                zbuffer[int(P.x+P.y*width_2)] = P.z;
                image.set(P.x, P.y, TGAColor(intensity*t_color.r, intensity*t_color.g, intensity*t_color.b, 255));
                
            }
        }
    }
}


int main(int argc, char** argv) {
    // TGAImage image(width, height, TGAImage::RGB);

    // Vec2i t0[3] = {Vec2i(70, 80),   Vec2i(50, 160),  Vec2i(10, 70)};
    // Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)};
    // Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)};
    // Vec3f a = barycentric(t0, Vec2i(51, 156));
    // std::cout << a.x << " " << a.y << " " << a.z << std::endl;
    // a = barycentric(t0, Vec2i(54, 156));
    // std::cout << a.x << " " << a.y << " " << a.z << std::endl;

    TGAImage texture;
    texture.read_tga_file("african_head_diffuse.tga");
    texture.flip_vertically();

    Model* model = new Model("obj/african_head.obj");
    float* zbuffer = new float[width_2*height_2];
    for (int i=0; i<width_2*height_2; i++) {
        zbuffer[i] = -std::numeric_limits<float>::max();
    }
    TGAImage image2(width_2, height_2, TGAImage::RGB);
    Vec3f light_dir(-1,0,0);
    light_dir.normalize();
    for (int i=0; i<model->nfaces(); i++) {
        //every face contain a triangle represented by 3 vertics.
        std::vector<int> face = model->face(i);
        Vec3f screen_coords[3];
        Vec3f world_coords[3];
        Vec2f texture_coords[3];
        for (int j=0; j<3; j++) {
            Vec3f v = model->vert(face[j]);
            Vec2f t_v = model->t_vert(model->texture(i)[j]);
            //this int() very important.
            screen_coords[j] = Vec3f(int((v.x+1)*width_2/2.+.5), int((v.y+1)*height_2/2.+.5), v.z);
            world_coords[j]  = v;
            texture_coords[j] = Vec2f(int(t_v.x*texture.get_width()), int(t_v.y*texture.get_height()));
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
            //triangle_2(screen_coords, zbuffer, image2, TGAColor(intensity*255, intensity*255, intensity*255, 255));
            trinagle_texture(screen_coords, texture_coords, zbuffer, image2, texture, intensity);
        } 
        //random fill color
        //triangle_2(screen_coords[0], screen_coords[1], screen_coords[2], image2, TGAColor(rand()%255, rand()%255, rand()%255, 255));
    }

    // image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    // image.write_tga_file("test.tga");
    image2.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image2.write_tga_file("face_w_light_normalize.tga");

    { // dump z-buffer (debugging purposes only)
        TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                zbimage.set(i, j, TGAColor(zbuffer[i+j*width], 1));
            }
        }
        zbimage.flip_vertically(); // i want to have the origin at the left bottom corner of the image
        zbimage.write_tga_file("zbuffer.tga");
    }

    delete [] zbuffer;
    delete model;
    return 0;
}