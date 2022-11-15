#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
#include<sstream>
#include<set>

using namespace std;

int CountNumbersOfTextLines(const string &filePath )
{
  long i = 0;
  ifstream ifs(filePath);
  if( ifs ){
    string line;
    while( true ){
      getline( ifs, line );
      i++;
      if(ifs.eof())
        break;
    }
  }
  return i-1;
}

int main()
{
    ifstream ifs_velocity("boundary_disp.dat");
    ifstream ifs2_element("element.dat");
    ifstream ifs_pressure("boundary_disp2.dat");

    int velocity_line=CountNumbersOfTextLines("boundary_disp.dat");
    int element_line=CountNumbersOfTextLines("element.dat");
    int pressure_line=CountNumbersOfTextLines("boundary_disp2.dat");

    set<int> boundary_node;
    vector<vector<int>> element;

    string str;
    for(int i=0; i<velocity_line; i++){
        getline(ifs_velocity,str);
        istringstream stream(str);
        vector<int> tmp_node;
        for(int j=0; j<4; j++){
            getline(stream,str,' ');
              if(j==0){
                boundary_node.insert(stoi(str));
              }
            }
            tmp_node.push_back(stoi(str));
    }

    for(int i=0; i<pressure_line; i++){
        getline(ifs_pressure,str);
        istringstream stream(str);
        vector<int> tmp_node;
        for(int j=0; j<4; j++){
            getline(stream,str,' ');
              if(j==0){
                boundary_node.insert(stoi(str));
              }
            }
            tmp_node.push_back(stoi(str));
    }

    for(int i=0; i<element_line; i++){
        getline(ifs2_element,str);
        istringstream stream(str);
        vector<int> tmp_element;
        while(getline(stream, str, ' ')){
            tmp_element.push_back(stoi(str));
        }
        element.push_back(tmp_element);
    }

    //for(int i=0; i<pressure_line; i++){
    //    getline(ifs_pressure,str);
    //    istringstream stream(str);
    //    vector<int> tmp;
    //    for(int j=0; j<3; j++){
    //        getline(stream,str,' ');
    //        if(j==0){
    //            boundary_node.insert(stoi(str));
    //        }
    //    }
    //}

    ofstream ofs("prism.dat");

    for(int i=0; i<element.size(); i++){
        cout << i << endl;
        for(int j=0; j<element[i].size(); j++){
            if(boundary_node.find(element[i][j])!=boundary_node.end()){
                if(element[i].size()==4) continue;
                for(int k=0; k<element[i].size(); k++){
                    ofs << element[i][k] << " ";
                    //cout << element[i][k] << " ";
                }
                //cout << " " << i <<endl;
                ofs << i <<endl;
                break;
            }
        }
    }
}