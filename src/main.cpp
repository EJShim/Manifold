#include <stdlib.h>
#include <stdio.h>
#include "Model_OBJ.h"
#include <Python.h>

//simplify
#include "igl/decimate.h"
#include "igl/qslim.h"
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


    __declspec(dllexport) PyObject* Simplify(double* posPtr, int posSize, int* facePtr, int faceSize, int max_faces=0x7fffffff){

          
        PyObject* dict = PyDict_New();

        using namespace std;
        using namespace Eigen;
        using namespace igl;
      
        double max_cost = 1e30;
        double max_ratio = 2.;
        bool check_manifold = false;

        MatrixXd OV = Map<MatrixXd>(posPtr,  3, posSize/3).transpose();
        MatrixXi OF = Map<MatrixXi>(facePtr, 3, faceSize/3).transpose();

        if (max_faces == 0x7fffffff) {
            max_faces = max_ratio * OF.rows() + 0.5;
        } else {
        if (max_ratio < 1)
            max_faces = std::max(max_faces, (int)(max_ratio * OF.rows() + 0.5));
        }
        
        std::cout << "max faces : " << max_faces << std::endl;

        auto MyDecimate = [&](
            const Eigen::MatrixXd & V,
            const Eigen::MatrixXi & F,
            const size_t max_m,
            const double max_cost,
            Eigen::MatrixXd & U,
            Eigen::MatrixXi & G,
            Eigen::VectorXi & J,
            Eigen::VectorXi & I) {
            // Original number of faces
            const int orig_m = F.rows();
            // Tracking number of faces
            int m = F.rows();
            typedef Eigen::MatrixXd DerivedV;
            typedef Eigen::MatrixXi DerivedF;
            DerivedV VO;
            DerivedF FO;
            igl::connect_boundary_to_infinity(V,F,VO,FO);
            // decimate will not work correctly on non-edge-manifold meshes. By extension
            // this includes meshes with non-manifold vertices on the boundary since these
            // will create a non-manifold edge when connected to infinity.
            if(!is_edge_manifold(FO))
            {
              return false;
            }
        
            auto stopping_condition = 
            [&](
            const Eigen::MatrixXd &V,
            const Eigen::MatrixXi &F,
            const Eigen::MatrixXi &E,
            const Eigen::VectorXi &EMAP,
            const Eigen::MatrixXi &EF,
            const Eigen::MatrixXi &EI,
            const std::set<std::pair<double,int> > &Q,
            const std::vector<std::set<std::pair<double,int> >::iterator > &QIT,
            const Eigen::MatrixXd &C,
            const int e,
            const int e1,
            const int e2,
            const int f1,
            const int f2)->bool
            {
              // Only subtract if we're collapsing a real face
              if(f1 < orig_m) m-=1;
              if(f2 < orig_m) m-=1;
              return (m<=(int)max_m) || (Q.begin()->first >= max_cost);
            };
        
            using namespace igl;
        
            Eigen::VectorXi EMAP;
            Eigen::MatrixXi E,EF,EI;
            edge_flaps(FO,E,EMAP,EF,EI);
            // Quadrics per vertex
            typedef std::tuple<Eigen::MatrixXd,Eigen::RowVectorXd,double> Quadric;
            std::vector<Quadric> quadrics;
            per_vertex_point_to_plane_quadrics(VO,FO,EMAP,EF,EI,quadrics);
            // State variables keeping track of edge we just collapsed
            int v1 = -1;
            int v2 = -1;
            // Callbacks for computing and updating metric
            std::function<void(
              const int e,
              const Eigen::MatrixXd &,
              const Eigen::MatrixXi &,
              const Eigen::MatrixXi &,
              const Eigen::VectorXi &,
              const Eigen::MatrixXi &,
              const Eigen::MatrixXi &,
              double &,
              Eigen::RowVectorXd &)> cost_and_placement;
            std::function<bool(
              const Eigen::MatrixXd &                                         ,/*V*/
              const Eigen::MatrixXi &                                         ,/*F*/
              const Eigen::MatrixXi &                                         ,/*E*/
              const Eigen::VectorXi &                                         ,/*EMAP*/
              const Eigen::MatrixXi &                                         ,/*EF*/
              const Eigen::MatrixXi &                                         ,/*EI*/
              const std::set<std::pair<double,int> > &                        ,/*Q*/
              const std::vector<std::set<std::pair<double,int> >::iterator > &,/*Qit*/
              const Eigen::MatrixXd &                                         ,/*C*/
              const int                                                        /*e*/
              )> pre_collapse;
            std::function<void(
              const Eigen::MatrixXd &                                         ,   /*V*/
              const Eigen::MatrixXi &                                         ,   /*F*/
              const Eigen::MatrixXi &                                         ,   /*E*/
              const Eigen::VectorXi &                                         ,/*EMAP*/
              const Eigen::MatrixXi &                                         ,  /*EF*/
              const Eigen::MatrixXi &                                         ,  /*EI*/
              const std::set<std::pair<double,int> > &                        ,   /*Q*/
              const std::vector<std::set<std::pair<double,int> >::iterator > &, /*Qit*/
              const Eigen::MatrixXd &                                         ,   /*C*/
              const int                                                       ,   /*e*/
              const int                                                       ,  /*e1*/
              const int                                                       ,  /*e2*/
              const int                                                       ,  /*f1*/
              const int                                                       ,  /*f2*/
              const bool                                                  /*collapsed*/
              )> post_collapse;
            qslim_optimal_collapse_edge_callbacks(
              E,quadrics,v1,v2, cost_and_placement, pre_collapse,post_collapse);
            // Call to greedy decimator
            bool ret = decimate(
              VO, FO,
              cost_and_placement,
              stopping_condition,
              pre_collapse,
              post_collapse,
              E, EMAP, EF, EI,
              U, G, J, I);
            // Remove phony boundary faces and clean up
            const Eigen::Array<bool,Eigen::Dynamic,1> keep = (J.array()<orig_m);
            igl::slice_mask(Eigen::MatrixXi(G),keep,1,G);
            igl::slice_mask(Eigen::VectorXi(J),keep,1,J);
            Eigen::VectorXi _1,I2;
            igl::remove_unreferenced(Eigen::MatrixXd(U),Eigen::MatrixXi(G),U,G,_1,I2);
            igl::slice(Eigen::VectorXi(I),I2,1,I);
        
            return ret;
          };

          MatrixXd V;
          MatrixXi F;
          VectorXi J, I;
          MyDecimate(OV, OF, max_faces, max_cost, V, F, J, I);


          std::cout << V.rows() << "," << V.cols() << std::endl;
          std::cout << F.rows() << "," << F.cols() << std::endl;


          write_triangle_mesh("test.obj",V,F);

          PyObject *vList = PyList_New(0);
          PyObject *fList = PyList_New(0);
  
          for(int i=0 ; i< V.rows() ; i++){
              PyList_Append(vList, Py_BuildValue("d", V.row(i)[0]));
              PyList_Append(vList, Py_BuildValue("d", V.row(i)[1]));
              PyList_Append(vList, Py_BuildValue("d", V.row(i)[2]));
          }
  
          for(int i=0 ; i< F.rows() ; i++){
              PyList_Append(fList, Py_BuildValue("i", F.row(i)[0]));
              PyList_Append(fList, Py_BuildValue("i", F.row(i)[1]));
              PyList_Append(fList, Py_BuildValue("i", F.row(i)[2]));
          }
  

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
