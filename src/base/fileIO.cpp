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

void fileIO::input_prism(string layer_file, vector<vector<int>> &layer_pair)
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