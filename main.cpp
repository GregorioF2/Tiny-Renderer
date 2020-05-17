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

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, const TGAColor &color) { 
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

void printFaceInTringles(TGAImage &image, Model* model) {
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
        Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]); 
        n.normalize(); 
        float intensity = n*light_dir; 
        if (intensity>0) { 
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity*255, intensity*255, intensity*255, 255)); 
        } 
    }
}

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    TGAImage image(width, height, TGAImage::RGB);
    //triangle(Vec2i(0,0), Vec2i(100, 200), Vec2i(200, 100), image, white);

    //printFaceInLines(image, model);
    printFaceInTringles(image, model);
    

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}

