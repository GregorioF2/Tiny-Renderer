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
Vec3f camera(0,0,3);



using namespace std;

Vec3f m2v(Matrix m) {
    return Vec3f(m[0][0]/m[3][0], m[1][0]/m[3][0], m[2][0]/m[3][0]);
}

Matrix v2m(Vec3f v) {
    Matrix m(4, 1);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

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

/*
Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec2f P) {
    Vec3f s1 (C.x-A.x, B.x-A.x, A.x-P.x);
    Vec3f s2 (C.y-A.y, B.y-A.y, A.y-P.y);


    Vec3f u = cross(s1, s2);
    if (std::abs(u.z)>1e-2)
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return Vec3f(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}


void triangle3DBaycentric(Vec3f v1, Vec3f v2, Vec3f v3, float *zbuffer, TGAImage &image, TGAColor color) {
    
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
*/

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

void triangle3D(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, TGAImage &image, float intensity, int *zbuffer) {
    if (t0.y > t1.y) {
        swap(t0, t1);
        swap(uv0, uv1);
    }
    if (t0.y > t2.y) {
        swap(t0, t2);
        swap(uv0, uv2);
    }
    if (t1.y > t2.y) {
        swap(t1, t2);
        swap(uv1, uv2);
    }
    float totalHeight = t2.y - t0.y;
    for (int i = 0; i<= totalHeight; i++) {
        bool firstPart = i <= t1.y - t0.y;
        float heightAlpha = totalHeight;
        float heightBetha = firstPart
            ? t1.y - t0.y
            : t2.y - t1.y;
        
        heightAlpha = heightAlpha == 0 ? 1 : heightAlpha;
        heightBetha = heightBetha == 0 ? 1 : heightBetha;

        float advancedAlpha = i/heightAlpha;
        float advancedBetha = firstPart
            ? (float) i / heightBetha
            : (float) (i - (t1.y - t0.y)) / heightBetha;

        
        Vec3i alpha = t0 + Vec3f(t2 - t0) * advancedAlpha;
        Vec3i betha = firstPart
            ? t0 + Vec3f(t1 - t0) * advancedBetha
            : t1 + Vec3f(t2 - t1) * advancedBetha;
        
        Vec2i uvAlpha = uv0 + (uv2 - uv0)*advancedAlpha;
        Vec2i uvBetha = firstPart
            ? uv0 + (uv1 - uv0)*advancedBetha
            : uv1 + (uv2 - uv1)*advancedBetha;

        if (alpha.x > betha.x) {
            swap(alpha, betha);
            swap(uvAlpha, uvBetha);
        }

        float distance = betha.x - alpha.x;
        for (int x = alpha.x; x <= betha.x; x++) {
            float advanced = distance < 0.05 ? 1. : (x-alpha.x) / distance;


            Vec3i point = Vec3f(alpha) + Vec3f(betha - alpha) * advanced;
            Vec2i texturePoint = uvAlpha + (uvBetha - uvAlpha) * advanced;
            
            int xZindex = point.x;
            int yZindex = point.y * width;

            if (point.z < zbuffer[xZindex + yZindex]) continue;

            zbuffer[xZindex + yZindex] = point.z;
            
            TGAColor color = model -> diffuse(texturePoint);
            image.set(
                point.x,
                point.y,
                TGAColor(color.r*intensity, color.g*intensity, color.b*intensity));
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

Vec3f world2screen(Vec3f v) {
    return Vec3f(int((v.x+1.)*width/2.+.5), int((v.y+1.)*height/2.+.5), v.z);
}

Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity(4);
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}

void printFaceInTringles3D(TGAImage &image, Model* model) {
    int *zbuffer = new int[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<int>::max());

    Matrix Projection = Matrix::identity(4);
    Matrix ViewPort   = viewport(width/8, height/8, width*3/4, height*3/4);
    Projection[3][2] = -1.f/camera.z;

    for (int faceIndex=0; faceIndex<model->nfaces(); faceIndex++) {
        std::vector<int> face = model->face(faceIndex);

        Vec3f literalCoords [3];
        literalCoords[0] = model->vert(face[0]);
        literalCoords[1] = model->vert(face[1]);
        literalCoords[2] = model->vert(face[2]);
        
        Vec3f coordsOnScreen[3];
        coordsOnScreen[0] =  m2v(ViewPort*Projection*v2m(literalCoords[0]));
        coordsOnScreen[1] =  m2v(ViewPort*Projection*v2m(literalCoords[1]));
        coordsOnScreen[2] =  m2v(ViewPort*Projection*v2m(literalCoords[2]));

        Vec3f normal = (literalCoords[2]-literalCoords[0])^(literalCoords[1]-literalCoords[0]);
        normal.normalize();

        float intensity = normal*light_dir;
        if (intensity <= 0) continue;

        Vec2i textures[3];
        for (int k=0; k<3; k++) {
            textures[k] = model->uv(faceIndex, k);
        }
        triangle3D(coordsOnScreen[0], coordsOnScreen[1], coordsOnScreen[2], textures[0], textures[1], textures[2], image, intensity, zbuffer);
    }

    TGAImage zImage (width, height, TGAImage::RGB);
    for(int i = 0; i < width; i++) {
        for(int j = 0; j < height; j++) {
            TGAColor c (zbuffer[i + j*width]*255, 1);
            zImage.set(i, j, c);
        }
    }
    zImage.flip_vertically();
    zImage.write_tga_file("zbuffer.tga");
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

