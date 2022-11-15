#include "fileIO.h"

using namespace std;

int fileIO::CountNumbersOfTextLines(const string &filePath )
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

void fileIO::input_prism(string layer_file, vector<vector<int>> &layer_pair, vector<int> &prism_id)
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
      for(int j=0;j<7;j++){
        getline(stream,tmp,' ');
          if(j==6){
            prism_id.push_back(stoi(tmp));
            continue;
          }
          
          layer_pair[i][j] = stoi(tmp);
      }
  }
  file.close();
}

void fileIO::input_element(std::string element_file, std::vector<std::vector<int>> &element)
{
    ifstream ifs(element_file);
    int numOfelement = CountNumbersOfTextLines(element_file);
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
}

void fileIO::export_vtu(const std::string &file, vector<vector<double>> x, vector<vector<int>> element, vector<double> element_wss_u,vector<double> element_wss_v,vector<double> element_wss_w)
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

  fprintf(fp, "</PointData>\n");


  fprintf(fp, "<CellData>\n");
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"wall_share_stress[Pa]\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * element.size()*3;
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
  delete data_d;

  double *data_d2 = new double[element.size()*3];
  num=0;
  for (int ic = 0; ic < element.size(); ic++){
    data_d2[num]   = element_wss_u[ic];
    data_d2[num+1]   = element_wss_v[ic];
    data_d2[num+2]   = element_wss_w[ic];
    num+=3;
  }
  size=sizeof(double)*element.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d2, size);

  delete data_d2;
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

void fileIO::export_vtu_OSI(const std::string &file, std::vector<std::vector<double>> x, std::vector<std::vector<int>> element, std::vector<double> OSI)
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

  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"OSI[-]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * element.size();
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
  delete data_d;

  double *data_d2 = new double[element.size()];
  num=0;
  for (int ic = 0; ic < element.size(); ic++){
    data_d2[num]   = OSI[ic];
    num++;
  }
  size=sizeof(double)*element.size();
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d2, size);

  delete data_d2;
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