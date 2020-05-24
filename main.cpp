#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "vertex.h"
#include <tuple>
#include <iostream>
#include <algorithm>
#include "utils.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor black = TGAColor(0, 0, 0 ,0);
Model *model = NULL;
const int width  = 800;
const int height = 800;
const int depth  = 255;
Vec3f light_dir(0,0,-1);



using namespace std;

void line(int x0, int y0, int x1, int y1, TGAImage &image, const TGAColor &color) {
    bool steep = false;
    if (std::abs(x0-x1)<std::abs(y0-y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0>x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    for (int x=x0; x<=x1; x++) {
        float t = (x-x0)/(float)(x1-x0);
        int y = floor(y0*(1.-t) + y1*t);
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec2f P) {
    Vec3f s1 (C.x-A.x, B.x-A.x, A.x-P.x);
    Vec3f s2 (C.y-A.y, B.y-A.y, A.y-P.y);


    Vec3f u = cross(s1, s2);
    if (std::abs(u.z)>1e-2)
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return Vec3f(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}


void triangle3D(Vec3f v1, Vec3f v2, Vec3f v3, float *zbuffer, TGAImage &image, TGAColor color) {
    
    Vec2f bboxMin = Vec2f(
        min(v1.x, min(v2.x, v3.x)),
        min(v1.y, min(v2.y, v3.y))
    );
    Vec2f bboxMax = Vec2f(
        max(v1.x, max(v2.x, v3.x)),
        max(v1.y, max(v2.y, v3.y))
    );

    for (int x = bboxMin.x; x <= bboxMax.x; x++) {
        for (int y = bboxMin.y; y <= bboxMax.y; y++) {
            Vec2f point2D = Vec2f(x, y);
            Vec3f barycentricCoords = barycentric(v1, v2, v3, point2D);
            
            if (barycentricCoords.x < 0 || barycentricCoords.y < 0 || barycentricCoords.z < 0)
                continue;

            float zValue = barycentricCoords.x * v1.z
                + barycentricCoords.y * v2.z
                + barycentricCoords.z * v3.z;
            
            int xIndexInZBuffer = point2D.x;
            int yIndexInZBuffer = point2D.y * width;

            if (zbuffer[xIndexInZBuffer + yIndexInZBuffer] > zValue)
                continue;

            zbuffer[xIndexInZBuffer + yIndexInZBuffer] = zValue;
            image.set(point2D.x, point2D.y, color);
        }
    }

}

void triangle2D(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, const TGAColor &color) { 
    if (t0.y > t1.y) swap(t0, t1);
    if (t0.y > t2.y) swap(t0, t2);
    if (t1.y > t2.y) swap(t1, t2);

    float partialHeight = t1.y - t0.y;
    float totalHeight = t2.y - t0.y;
    float diferenceHeight = partialHeight/totalHeight;
    for (int y = t0.y; y <= t1.y; y++) {
        float advance = (y - t0.y) / partialHeight;
        float longAdvance = advance * diferenceHeight;

        Vec2i v1 = t0 + ((t2 - t0) * longAdvance);
        Vec2i v2 = t0 + ((t1 - t0) * advance);
        if (v1.x > v2.x) swap(v1, v2);
        for (int j=v1.x; j<=v2.x; j++) { 
            image.set(j, y, color);
        } 
    }

    partialHeight = t2.y - t1.y;
    diferenceHeight = partialHeight/totalHeight;
    for (int y = t2.y; y >= t1.y; y--) {
        float advance = abs((y - t2.y) / (float)(t2.y - t1.y));
        float longAdvance = advance * diferenceHeight;

        Vec2i v1 = t2 - ((t2 - t0) * longAdvance);
        Vec2i v2 = t2 - ((t2 - t1) * advance);
        if (v1.x > v2.x) swap(v1, v2);
        for (int j = v1.x; j <= v2.x; j++) {
            image.set(j, y, color);
        } 
    }
}



void printFaceInLines(TGAImage &image, Model* model) {
    for (int i=0; i<model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        for (int j=0; j<3; j++) {
            Vec3f v0 = model->vert(face[j]);
            Vec3f v1 = model->vert(face[(j+1)%3]);
            int x0 = (v0.x+1.)*width/2.;
            int y0 = (v0.y+1.)*height/2.;
            int x1 = (v1.x+1.)*width/2.;
            int y1 = (v1.y+1.)*height/2.;
            line(x0, y0, x1, y1, image, white);
        }
    }
}

void printFaceInTringles2D(TGAImage &image, Model* model) {
    Vec3f light_dir(0,0,-1);
        for (int i=0; i<model->nfaces(); i++) { 
        std::vector<int> face = model->face(i); 
        Vec2i screen_coords[3]; 
        Vec3f world_coords[3]; 
        for (int j=0; j<3; j++) { 
            Vec3f v = model->vert(face[j]); 
            screen_coords[j] = Vec2i((v.x+1.)*width/2., (v.y+1.)*height/2.); 
            world_coords[j]  = v; 
        } 
        Vec3f n;// = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]); 
        n.normalize(); 
        float intensity = n*light_dir; 
        if (intensity>0) { 
            triangle2D(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity*255, intensity*255, intensity*255, 255)); 
        } 
    }
}


void test(char* filename) {
    TGAImage image;
    cout<< "file name: " << filename<< endl;
    bool rightRead = image.read_tga_file("./obj/african_head/african_head_diffuse.tga");
    if(!rightRead) {
        cout << "Error al lee TGA image" << endl;
        return;
    }
    image.write_tga_file("prueba.tga");

}

Vec3f world2screen(Vec3f v) {
    return Vec3f(int((v.x+1.)*width/2.+.5), int((v.y+1.)*height/2.+.5), v.z);
}

void printFaceInTringles3D(TGAImage &image, Model* model) {
    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    for (int i=0; i<model->nfaces(); i++) {
        std::vector<int> face = model->face(i);

         Vec3f screen_coords[3];
        Vec3f world_coords[3];
        for (int j=0; j<3; j++) {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = Vec3f((v.x+1.)*width/2., (v.y+1.)*height/2., (v.z+1.)*depth/2.);
            world_coords[j]  = v;
        }
        Vec3f n = cross((world_coords[2]-world_coords[0]), (world_coords[1]-world_coords[0]));
        n.normalize();
        float intensity = n*light_dir;
        if (intensity>0) {
            triangle3D(screen_coords[0], screen_coords[1], screen_coords[2], zbuffer, image, TGAColor(intensity*255, intensity*255, intensity*255, 200));
        }
    }

    TGAImage zImage (width, height, TGAImage::RGB);
    for(int i = 0; i < width; i++) {
        for(int j = 0; j < height; j++) {
            TGAColor c (zbuffer[i + j*width]*255, 1);
            zImage.set(i, j, c);
        }
        zImage.write_tga_file("zbuffer.tga");
    }
     delete [] zbuffer;
}

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head/african_head.obj");
    }
    
    // test("test string");
    

    TGAImage image(width, height, TGAImage::RGB);
    // triangle(Vec2i(0,0), Vec2i(100, 200), Vec2i(200, 100), image, white);

    // printFaceInLines(image, model);
    // printFaceInTringles2D(image, model);

    printFaceInTringles3D(image, model);

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}

