#include <H5Cpp.h>
#include <cmath>
#include <fstream>
#include <glob.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
#include "gauss.h"
#include "ShapeFunction.h"
#include "fem_base_mathTool.h"
#include <sys/stat.h>

using namespace std;
using namespace H5;

int CountNumbersOfTextLines(const string &filePath )
{
  long i = 0;

  ifstream ifs( filePath );

  if( ifs ){
    string line;

    while( true ){
      getline( ifs, line );
      i++;
      if( ifs.eof() )
        break;
    }
  }
  return i-1;
}

void calc_wall_share_stress(vector<vector<int>> layer_pair, vector<vector<double>> x, vector<double> u, vector<double> v, vector<double> w, vector<double> &wall_share_stress_u, vector<double> &wall_share_stress_v, vector<double> &wall_share_stress_w)
{
    double mu =0.001;
    string str,tmp;
    int numOfLayer = layer_pair.size();

    GaussWedge gWed(1);
    for (int i = 0; i < x.size(); i++){
            wall_share_stress_u[i] = 0.0;
            wall_share_stress_v[i] = 0.0;
            wall_share_stress_w[i] = 0.0;
    }

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

        ShapeFunction3D::C3D6_N(N, gWed.point[0][0], gWed.point[0][1], gWed.point[0][2], gWed.point[0][0]);
        ShapeFunction3D::C3D6_dNdr(dNdr, gWed.point[0][0], gWed.point[0][1], gWed.point[0][2], gWed.point[0][0]);
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
        
        for (int i = 0; i < 3; i++){
            wall_share_stress_u[layer_pair[ic][i]] = tproy[0];
            wall_share_stress_v[layer_pair[ic][i]] = tproy[1];
            wall_share_stress_w[layer_pair[ic][i]] = tproy[2];
        }
    }
}

void export_vtu(const std::string &file, vector<vector<double>> x, vector<vector<int>> element, vector<double> u, vector<double> v, vector<double> w, vector<double> wall_share_stress_u, vector<double> wall_share_stress_v, vector<double> wall_share_stress_w)
{
  FILE *fp;

  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", x.size(), element.size());
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * x.size() * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    for (int j = 0; j < element[i].size(); j++) fprintf(fp, "%d ", element[i][j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < element.size(); i++)
  {
    num += element[i].size();
    fprintf(fp, "%d\n", num);
  }

  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    if(element[i].size()==4) fprintf(fp, "%d\n", 10);
    if(element[i].size()==6) fprintf(fp, "%d\n", 13);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"velocity[m/s]\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * x.size() * 3;

  fprintf(fp, "<DataArray type=\"Float64\" Name=\"wall_share_stress[Pa]\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * x.size() * 3;
  fprintf(fp, "</PointData>\n");


  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);
   

  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  double *data_d = new double[x.size()*3];
  num = 0;
  int size=0;
  for (int ic = 0; ic < x.size(); ic++){
    for(int j=0;j<3;j++){
      data_d[num] = x[ic][j];
      num++;
    }
  }
  size=sizeof(double)*x.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  num=0;
  for (int ic = 0; ic < x.size(); ic++){
      data_d[num]   = u[ic];
      data_d[num+1] = v[ic];
      data_d[num+2] = w[ic];
      num=num+3;
  }
  size=sizeof(double)*x.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  num=0;
  for (int ic = 0; ic < x.size(); ic++){
    data_d[num]   = wall_share_stress_u[ic];
    data_d[num+1] = wall_share_stress_v[ic];
    data_d[num+2] = wall_share_stress_w[ic];
    num=num+3;
  }
  size=sizeof(double)*x.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  delete data_d;
  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

vector<double> importHDF5_double_H5_2D(H5::H5File &file, const std::string &dataName, int &nx_t, int &ny_t) {
  H5::DataSet dataset = file.openDataSet(dataName.c_str());
  H5::DataSpace dataspace = dataset.getSpace();

  hsize_t dims_out[2];
  int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

  nx_t = dims_out[0];
  ny_t = dims_out[1];

  int nx = (unsigned long)(dims_out[0]);
  int ny = (unsigned long)(dims_out[1]);

  double *data;
  data = new double[nx*ny];

  dataset.read(&data[0], H5::PredType::NATIVE_DOUBLE);
  //cout << nx << " " << ny << endl;
  vector<double> d;
  for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++){
          d.push_back(data[ny*i+j]);
          //cout << data[nx*i+j] << endl;
      }
      //exit(1);
  }
  delete[] data;
  return d;
}

vector<double> importHDF5_double_H5_1D(H5::H5File &file, const std::string &dataName) {
    H5::DataSet dataset = file.openDataSet(dataName.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[2];
    int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
    int nx = (unsigned long)(dims_out[0]);
    double *data;
    data = new double[nx];
    
    dataset.read(&data[0], H5::PredType::NATIVE_DOUBLE);

    vector<double> d;
    for (int i = 0; i < nx; i++){
        d.push_back(data[i]);
    }
    delete[] data;

    return d;
}

vector<double> readH5_double_2D(std::string h5file, const int number, string data_n, int &nx, int &ny) {
  string dataName;
  H5::H5File file(h5file, H5F_ACC_RDONLY);

  dataName = to_string(number)+"/" + data_n;
  vector<double> d = importHDF5_double_H5_2D(file, dataName, nx, ny);
  return d;
}

vector<double> readH5_double_1D(std::string h5file, const int number, string data_n) {
  string dataName;
  H5::H5File file(h5file, H5F_ACC_RDONLY);

  dataName = to_string(number)+"/" + data_n;
  vector<double> d=importHDF5_double_H5_1D(file, dataName);
  return d;
}

void exportHDF5_double_1D(H5::H5File &file, const std::string &dataName, vector<double> i_data, int i_dim) {
  H5std_string DATASET_NAME(dataName.c_str());

  hsize_t dim[1] = {i_dim}; // dataset dimensions
  H5::DataSpace dataspace(1, dim);

  double *data;
  data = new double[i_data.size()];

  for (int i = 0; i < i_data.size(); i++) {
    data[i] = i_data[i];
  }

  H5::IntType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);
  H5::DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);
  dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE);
  delete[] data;
}

void input_prism(string layer_file, vector<vector<int>> &layer_pair)
{
  string str,tmp;
  int numOfLayer = CountNumbersOfTextLines(layer_file);
  layer_pair.resize(numOfLayer);
  for(int i=0; i<layer_pair.size(); i++){
    layer_pair[i].resize(6);
  }
  
  for(int ic=0;ic<numOfLayer;ic++){
    for(int j=0;j<6;j++) layer_pair[ic][j] = 0e0;
  }
    
  ifstream file(layer_file);
  if(!file){
      cout << "Error:Input "<< layer_file << " not found" << endl;
      exit(1);
  }
  
  for(int i=0;i<numOfLayer;i++){
      getline(file,str);
      istringstream stream(str);
      for(int j=0;j<6;j++){
          getline(stream,tmp,' ');
          layer_pair[i][j] = stoi(tmp);
      }
  }
  file.close();
}

int main(int argc, char *argv[]) {
    string input_dir = argv[1];
    cout << input_dir << endl;
    vector<vector<double>> u,v,w;
    vector<double> wall_share_stress_u,wall_share_stress_v,wall_share_stress_w;
    vector<vector<vector<double>>> x;
    int h5_time_number = 393;
    int Directions_x, Directions_y;

    string h5_file = input_dir + "/" + input_dir+".h5";

    for (int i = 1; i <= h5_time_number; i++){
        vector<double> u_tmp = readH5_double_1D(h5_file, i, "u");
        vector<double> v_tmp = readH5_double_1D(h5_file, i, "v");
        vector<double> w_tmp = readH5_double_1D(h5_file, i, "w");
        vector<double> Directions_tmp = readH5_double_2D(h5_file, i, "x", Directions_x, Directions_y);
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

    ifstream ifs(element_file);
    int numOfelement = CountNumbersOfTextLines(element_file);
    vector<vector<int>> element;
    for(int i=0; i<numOfelement; i++){
        string str;
        getline(ifs,str);
        istringstream stream(str);
        vector<int> tmp_element;
        while(getline(stream,str,' ')){
            tmp_element.push_back(stoi(str));
        }
        element.push_back(tmp_element);
    }

    wall_share_stress_u.resize(x[0].size());
    wall_share_stress_v.resize(x[0].size());
    wall_share_stress_w.resize(x[0].size());

    string layer_file = input_dir + "/" + "prism.dat";
    vector<vector<int>> layer_pair;
    input_prism(layer_file, layer_pair);
    
    string result_folder = "Result_" + input_dir;
    mkdir(result_folder.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    vector<vector<double>> wss_u, wss_v, wss_w;

    for(int i=0; i<u.size(); i++){
        calc_wall_share_stress(layer_pair, x[i], u[i], v[i], w[i], wall_share_stress_u, wall_share_stress_v, wall_share_stress_w);
        string vtu_file_path = result_folder + "/test_" + to_string(i) + ".vtu";
        cout << vtu_file_path << endl;
        export_vtu(vtu_file_path,x[i],element,u[i],v[i],w[i], wall_share_stress_u, wall_share_stress_v, wall_share_stress_w);
        wss_u.push_back(wall_share_stress_u);
        wss_v.push_back(wall_share_stress_v);
        wss_w.push_back(wall_share_stress_w);
    }

    cout << wss_u.size() << " " << wss_u[0].size() << endl;
    string output_h5_name = input_dir + "_wss.h5";
    H5std_string FILE_NAME(output_h5_name.c_str());
    for (int i = 0; i < wss_u.size(); i++){
      if (i == 0) {
        H5File file(FILE_NAME, H5F_ACC_TRUNC);
        std::string dataName;
        std::string Gr = "/"+to_string(i);
        file.createGroup(Gr.c_str());
        Group group = file.openGroup(Gr.c_str());
        dataName = Gr + "/wss_u";
        exportHDF5_double_1D(file, dataName, wss_u[i], wss_u[i].size());
        dataName = Gr + "/wss_v";
        exportHDF5_double_1D(file, dataName, wss_v[i], wss_u[i].size());
        dataName = Gr + "/wss_w";
        exportHDF5_double_1D(file, dataName, wss_w[i], wss_u[i].size());
      }
      else {
          H5File file(FILE_NAME, H5F_ACC_RDWR);
          std::string dataName;
          std::string Gr = "/" + to_string(i);
          file.createGroup(Gr.c_str());
          Group group = file.openGroup(Gr.c_str());
          dataName = Gr + "/wss_u";
          exportHDF5_double_1D(file, dataName, wss_u[i], wss_u[i].size());
          dataName = Gr + "/wss_v";
          exportHDF5_double_1D(file, dataName, wss_v[i], wss_u[i].size());
          dataName = Gr + "/wss_w";
          exportHDF5_double_1D(file, dataName, wss_w[i], wss_u[i].size());
      }
    }
}