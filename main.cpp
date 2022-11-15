#include <cmath>
#include <fstream>
#include <glob.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include "gauss.h"
#include "ShapeFunction.h"
#include "fem_base_mathTool.h"
#include <sys/stat.h>
#include<set>
#include "hdf5_rw.h"
#include "fileIO.h"
#include <H5Cpp.h>

using namespace std;
using namespace H5;

void calc_wall_share_stress(vector<vector<int>> layer_pair, vector<vector<double>> x, vector<double> u, vector<double> v, vector<double> w, vector<double> &wss_element_u, vector<double> &wss_element_v,vector<double> &wss_element_w,vector<vector<vector<double>>> &OSI_t, vector<int> prism_id)
{
    double mu =0.001;
    string str,tmp;
    int numOfLayer = layer_pair.size();

    GaussTriangle gTri(1);

    vector<vector<double>> tmp_t;
    for (int ic = 0; ic < numOfLayer; ic++){
        int numOfNodeInElm = 6;
        vector<vector<double>> x_current(numOfNodeInElm,vector<double>(3));
        vector<double> N(numOfNodeInElm);
        vector<vector<double>> dNdr(numOfNodeInElm, vector<double>(3));
        vector<vector<double>> dNdx(numOfNodeInElm, vector<double>(3));
        double dxdr[3][3],drdx[3][3];
        for(int i=0;i<numOfNodeInElm;i++){
            for(int j=0;j<3;j++){
                x_current[i][j] = x[layer_pair[ic][i]][j];
            }
        }

        ShapeFunction3D::C3D6_N(N, gTri.point[0][0], gTri.point[0][1], gTri.point[0][2],0.0);
        ShapeFunction3D::C3D6_dNdr(dNdr, gTri.point[0][0], gTri.point[0][1], gTri.point[0][2], 0.0);
        FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
        FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
        
        vector<double> t_x(3), t_y(3), t_z(3);
        double A0 = x[layer_pair[ic][0]][0]; double A1 = x[layer_pair[ic][0]][1]; double A2 = x[layer_pair[ic][0]][2];
        double B0 = x[layer_pair[ic][1]][0]; double B1 = x[layer_pair[ic][1]][1]; double B2 = x[layer_pair[ic][1]][2];
        double C0 = x[layer_pair[ic][2]][0]; double C1 = x[layer_pair[ic][2]][1]; double C2 = x[layer_pair[ic][2]][2];
        
        double a1 = B0 - A0; double a2 = B1 - A1; double a3 = B2 - A2;
        double b1 = C0 - A0; double b2 = C1 - A1; double b3 = C2 - A2;
        vector<double> n(3);
        n[0] = a2 * b3 - a3 * b2;
        n[1] = a3 * b1 - a1 * b3;
        n[2] = a1 * b2 - a2 * b1;
        
        double n_a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        n[0] = -n[0] / n_a;  
        n[1] = -n[1] / n_a;  
        n[2] = -n[2] / n_a;  
        for (int i = 0; i < 6; i++){
            t_x[0] += mu*(dNdx[i][0] * u[layer_pair[ic][i]] + dNdx[i][0] * u[layer_pair[ic][i]]);
            t_x[1] += mu*(dNdx[i][1] * u[layer_pair[ic][i]] + dNdx[i][0] * v[layer_pair[ic][i]]);
            t_x[2] += mu*(dNdx[i][2] * u[layer_pair[ic][i]] + dNdx[i][0] * w[layer_pair[ic][i]]);
            t_y[0] += mu*(dNdx[i][0] * v[layer_pair[ic][i]] + dNdx[i][1] * u[layer_pair[ic][i]]);
            t_y[1] += mu*(dNdx[i][1] * v[layer_pair[ic][i]] + dNdx[i][1] * v[layer_pair[ic][i]]);
            t_y[2] += mu*(dNdx[i][2] * v[layer_pair[ic][i]] + dNdx[i][1] * w[layer_pair[ic][i]]);
            t_z[0] += mu*(dNdx[i][0] * w[layer_pair[ic][i]] + dNdx[i][2] * u[layer_pair[ic][i]]);
            t_z[1] += mu*(dNdx[i][1] * w[layer_pair[ic][i]] + dNdx[i][2] * v[layer_pair[ic][i]]);
            t_z[2] += mu*(dNdx[i][2] * w[layer_pair[ic][i]] + dNdx[i][2] * w[layer_pair[ic][i]]);
        }
        vector<double> t(3);
        t[0] = t_x[0] * n[0] + t_x[1] * n[1] + t_x[2] * n[2];
        t[1] = t_y[0] * n[0] + t_y[1] * n[1] + t_y[2] * n[2];
        t[2] = t_z[0] * n[0] + t_z[1] * n[1] + t_z[2] * n[2];
        
        vector<double> tp_1(3);
        vector<double> tproy(3);
        tp_1[0] = t[1] * n[2] - t[2] * n[1];
        tp_1[1] = t[2] * n[0] - t[0] * n[2];
        tp_1[2] = t[0] * n[1] - t[1] * n[0];
        
        tproy[0] = n[1] * tp_1[2] - n[2] * tp_1[1];
        tproy[1] = n[2] * tp_1[0] - n[0] * tp_1[2];
        tproy[2] = n[0] * tp_1[1] - n[1] * tp_1[0];

        wss_element_u[prism_id[ic]] = tproy[0];
        wss_element_v[prism_id[ic]] = tproy[1];
        wss_element_w[prism_id[ic]] = tproy[2];
        vector<double> element_t(3);
        for(int i=0; i<3; i++) element_t[i]=t[i];
        tmp_t.push_back(element_t);
    }
    OSI_t.push_back(tmp_t);
}

void calc_OSI(vector<double> &OSI, vector<vector<vector<double>>> OSI_t, vector<vector<int>> layer_pair, int fourth_CC, vector<int> prism_id)
{
  cout << "check1" << endl;
  int start= fourth_CC/4*2;
  int end =start+fourth_CC/4;

  cout << "check2" << endl;
  vector<vector<double>> int_vec_t(layer_pair.size(), vector<double>(3,0));
  vector<double> abs_int_vec_t(layer_pair.size());
  
  cout << "check3" << endl;
  for(int i=0; i<layer_pair.size(); i++){
    for(int j=start; j<=end; j++){
      cout << i << " " << j << endl;
      int_vec_t[i][0] += OSI_t[j][i][0];
      int_vec_t[i][1] += OSI_t[j][i][1];
      int_vec_t[i][2] += OSI_t[j][i][2];
    }
    abs_int_vec_t[i] = sqrt(int_vec_t[i][0]*int_vec_t[i][0]+int_vec_t[i][1]*int_vec_t[i][1]+int_vec_t[i][2]*int_vec_t[i][2]);
  }
  cout << "check4" << endl;

  vector<double> int_abs_vec_t(layer_pair.size(),0);

  for(int i=0; i<layer_pair.size(); i++){
    for(int j=start; j<=end; j++){
      double abs_vec_t = sqrt(OSI_t[j][i][0]*OSI_t[j][i][0]+OSI_t[j][i][1]*OSI_t[j][i][1]+OSI_t[j][i][2]*OSI_t[j][i][2]);
      int_abs_vec_t[i] += abs_vec_t;
    }
  }

  cout << "check5" << endl;
  for(int i=0; i<layer_pair.size(); i++){
    cout << i << " " << prism_id[i] << endl;
    OSI[prism_id[i]] = 0.5*(1.0-abs_int_vec_t[i]/int_abs_vec_t[i]);
  }
}

int main(int argc, char *argv[]) {

    HDF5RW H5RW;
    fileIO file;
    string input_dir = argv[1];
    string patient_name = argv[2];
    int h5_time_number = stoi(argv[3]);
    if(argc!=4){
      cout << "input error!" << endl;
      exit(1);
    }
    cout << input_dir << endl;
    cout << h5_time_number << endl;
    vector<vector<double>> u,v,w;
    vector<double> element_wss_u,element_wss_v,element_wss_w;
    vector<vector<vector<double>>> x;
    int Directions_x, Directions_y;

    string h5_file = input_dir + "/" + patient_name+".h5";

    for (int i = 1; i <= h5_time_number; i++){
        vector<double> u_tmp = H5RW.readH5_double_1D(h5_file, i, "u");
        vector<double> v_tmp = H5RW.readH5_double_1D(h5_file, i, "v");
        vector<double> w_tmp = H5RW.readH5_double_1D(h5_file, i, "w");
        vector<double> Directions_tmp = H5RW.readH5_double_2D(h5_file, i, "x", Directions_x, Directions_y);
        vector<vector<double>> x_tmp(Directions_x, vector<double>(Directions_y));
        for(int j=0; j<Directions_x; j++){
            for(int k=0; k<3; k++){
                x_tmp[j][k] = Directions_tmp[j*3+k];
            }
        }
        u.push_back(u_tmp);
        v.push_back(v_tmp);
        w.push_back(w_tmp);
        x.push_back(x_tmp);
    }

    string element_file = input_dir + "/element.dat";
    vector<vector<int>> element;
    file.input_element(element_file, element);

    string layer_file = input_dir + "/" + "prism.dat";
    vector<vector<int>> layer_pair;
    vector<int> prism_id;
    file.input_prism(layer_file, layer_pair, prism_id);

    //for(int i=0; i<prism_id.size(); i++){
    //  cout << prism_id[i] << endl;
    //}

    string result_folder = "Result_" + patient_name;
    mkdir(result_folder.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    vector<vector<double>> t;
    vector<double> OSI;

    element_wss_u.resize(element.size());
    element_wss_v.resize(element.size());
    element_wss_w.resize(element.size());
    OSI.resize(element.size());
    //cout << element_wss.size() << endl;
    //exit(1);
    vector<vector<vector<double>>> OSI_t;
    //cout << "check" << endl;
    for(int i=0; i<u.size(); i++){
        calc_wall_share_stress(layer_pair, x[i], u[i], v[i], w[i], element_wss_u,element_wss_v,element_wss_w,OSI_t, prism_id);
        string vtu_file_path = result_folder + "/test_" + to_string(i) + ".vtu";
        cout << vtu_file_path << endl;
        //for(int j=0; j<element_wss.size(); j++){
        //  cout << element_wss[j] << endl;
        //}
        file.export_vtu(vtu_file_path,x[i],element,element_wss_u,element_wss_v,element_wss_w);
    }

    calc_OSI(OSI, OSI_t, layer_pair, h5_time_number, prism_id);
    file.export_vtu_OSI("OSI_test.vtu",x[0],element,OSI);


    //string output_h5_name = patient_name + "_wss.h5";
    //H5std_string FILE_NAME(output_h5_name.c_str());
    //for (int i = 0; i < wss_u.size(); i++){
    //  if (i == 0) {
    //    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    //    std::string dataName;
    //    std::string Gr = "/"+to_string(i);
    //    file.createGroup(Gr.c_str());
    //    Group group = file.openGroup(Gr.c_str());
    //    dataName = Gr + "/wss_u";
    //    H5RW.exportHDF5_double_1D(file, dataName, wss_u[i], wss_u[i].size());
    //    dataName = Gr + "/wss_v";
    //    H5RW.exportHDF5_double_1D(file, dataName, wss_v[i], wss_u[i].size());
    //    dataName = Gr + "/wss_w";
    //    H5RW.exportHDF5_double_1D(file, dataName, wss_w[i], wss_u[i].size());
    //  }
    //  else {
    //      H5File file(FILE_NAME, H5F_ACC_RDWR);
    //      std::string dataName;
    //      std::string Gr = "/" + to_string(i);
    //      file.createGroup(Gr.c_str());
    //      Group group = file.openGroup(Gr.c_str());
    //      dataName = Gr + "/wss_u";
    //      H5RW.exportHDF5_double_1D(file, dataName, wss_u[i], wss_u[i].size());
    //      dataName = Gr + "/wss_v";
    //      H5RW.exportHDF5_double_1D(file, dataName, wss_v[i], wss_u[i].size());
    //      dataName = Gr + "/wss_w";
    //      H5RW.exportHDF5_double_1D(file, dataName, wss_w[i], wss_u[i].size());
    //  }
    //}
}