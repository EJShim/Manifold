#include <stdlib.h>
#include <stdio.h>
#include "Model_OBJ.h"
#include <Python.h>


class Foo{
public:
    void bar(){
        std::cout << "Hell" << std::endl;
    }
};

extern "C" {
    __declspec(dllexport) void CalculateOBJ(char* path){
        std::cout << path << std::endl;
        Model_OBJ obj;
        int resolution = 20000;
        obj.Load(path);

        obj.Process_Manifold(resolution);
        obj.SaveOBJ("temp.obj");
    
    }


    __declspec(dllexport) PyObject* Calculate(double* posPtr, int posSize, int* facePtr, int faceSize, int resolution){ 
        std::vector<int> faces(facePtr, facePtr + faceSize);
        std::vector<double> points(posPtr, posPtr+posSize);

        int nPoints = points.size()/3;
        int nFaces = faces.size()/3;

        // Load OBJ
        Model_OBJ obj;
        obj.vertices.resize(nPoints);
        obj.face_indices.resize(nFaces);

        for(int i=0 ; i<nPoints ; i++){
            obj.vertices[i] = glm::dvec3( points[i*3], points[i*3+1], points[i*3+2] );
        }

        for(int i=0 ; i<nFaces ; i++){
            obj.face_indices[i] = glm::ivec3( faces[i*3], faces[i*3+1], faces[i*3+2] );
        }
        
        //Test Process
        obj.Process_Manifold(resolution);
        

        PyObject *vList = PyList_New(0);
        PyObject *fList = PyList_New(0);

        for(int i=0 ; i<obj.vertices.size() ; i++){
            PyList_Append(vList, Py_BuildValue("d", obj.vertices[i][0]));
            PyList_Append(vList, Py_BuildValue("d", obj.vertices[i][1]));
            PyList_Append(vList, Py_BuildValue("d", obj.vertices[i][2]));
        }

        for(int i=0 ; i<obj.face_indices.size() ; i++){
            PyList_Append(fList, Py_BuildValue("i", obj.face_indices[i][0]));
            PyList_Append(fList, Py_BuildValue("i", obj.face_indices[i][1]));
            PyList_Append(fList, Py_BuildValue("i", obj.face_indices[i][2]));
        }


        PyObject* dict = PyDict_New();
        PyDict_SetItem(dict, Py_BuildValue("s", "vertices"), vList);
        PyDict_SetItem(dict, Py_BuildValue("s", "faces"), fList);

        return dict;
    }
};


extern int g_sharp;
int main(int argc, char** argv)
{
    Model_OBJ obj;
    int resolution = 20000;
    if (argc < 3)
    {
        cout << "./manifold input.obj output.obj [resolution=20000] [-s]\n";
        return 0;
    }
    obj.Load(argv[1]);
    
    
    if (argc > 3)
    {
        if (strcmp(argv[3], "-s") == 0) {
            g_sharp = 1;
        } else {
            sscanf(argv[3], "%d", &resolution);
            if (argc > 4 && strcmp(argv[4], "-s") == 0) {
                g_sharp = 1;
            }
        }
    }
    printf("manifold %s %s %d\n", argv[1], argv[2], resolution);
    
    obj.Process_Manifold(resolution);
    obj.SaveOBJ(argv[2]);
    return 0; 
}
