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



void line2(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    int xSign = sign(x1 - x0);
    int ySign = sign(y1 - y0);

    // dx = delta x
    int dx = abs(x1 - x0);
    int dy = abs(y1 - y0);

    // xs = x step
    for (int xs = 0; xs < dx; xs ++) {
        float advanced = xs / (float)dx;
        float missing = 1 - advanced;
        
        int ys = abs(y0 * advanced) + ySign * abs(y1 * missing);
    }
}

void line(int x0, int y0, int x1, int y1, vector< tuple<int, int> > &coordinates) {
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
        int y = y0*(1.-t) + y1*t;
        if (steep) {
            coordinates.push_back(tuple<int, int>(y, x));
        } else {
            coordinates.push_back(tuple<int, int>(x, y));
        }
    }
}

void paintLine(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    vector< tuple<int, int> > coordinates;
    line(x0, y0, x1, y1, coordinates);
    for (int i = 0; i < coordinates.size(); i++) {
        int x = get<0>(coordinates[i]);
        int y = get<1>(coordinates[i]);
        image.set(x, y, color);
    }
}

void pushYValuesForXPosition(int x, vector< tuple<int, int> > &coords, vector<int> &res) {
    const int size = coords.size();
    for (int i = 0; i < size; i++) {
        if (get<0>(coords[i]) == x) {
            res.push_back(get<1>(coords[i]));
        }
    }
    return;
}


void triangle(Vec2i v1, Vec2i v2, Vec2i v3, TGAImage &image, const TGAColor &color) {
    vector< tuple<int, int> > vectorV1V2;
    vector< tuple<int, int> > vectorV1V3;
    vector< tuple<int, int> > vectorV2V3;
    line(v1.x, v1.y, v2.x, v2.y, vectorV1V2);
    line(v1.x, v1.y, v3.x, v3.y, vectorV1V3);
    line(v2.x, v2.y, v3.x, v3.y, vectorV2V3);

    int bboxX0 = min(min(v1.x, v2.x), v3.x);
    int bboxX1 = max(max(v1.x, v2.x), v3.x);

    int bboxY0 = min(min(v1.y, v2.y), v3.y);
    int bboxY1 = max(max(v1.y, v2.y), v3.y);

    for (int x = bboxX0; x <= bboxX1; x++) {
        vector<int> yValues;
        pushYValuesForXPosition(x, vectorV1V2, yValues);
        pushYValuesForXPosition(x, vectorV1V3, yValues);
        pushYValuesForXPosition(x, vectorV2V3, yValues);
        
        int yMax =  maxElem(yValues);
        int yMin =  minElem(yValues);

        for (int y = yMin; y <= yMax; y ++) {
            image.set(x, y, color);
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
            paintLine(x0, y0, x1, y1, image, white);
        }
    }
}

void printFaceInTringles(TGAImage &image, Model* model) {
    for (int i=0; i<model->nfaces(); i++) { 
        std::vector<int> face = model->face(i); 
        Vec2i screen_coords[3]; 
        for (int j=0; j<3; j++) { 
            Vec3f world_coords = model->vert(face[j]); 
            screen_coords[j] = Vec2i((world_coords.x+1.)*width/2., (world_coords.y+1.)*height/2.); 
        } 
        triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(rand()%255, rand()%255, rand()%255, 255)); 
    }
}

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    TGAImage image(width, height, TGAImage::RGB);
    // triangle(Vec2i(0,0), Vec2i(100, 200), Vec2i(200, 100), image, white);

    //printFaceInLines(image, model);
    printFaceInTringles(image, model);
    

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}

