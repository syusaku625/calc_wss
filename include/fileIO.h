#ifndef _fileIO_H_
#define _fileIO_H_

#include<string>
#include<vector>
#include<sstream>
#include<fstream>
#include<iostream>

class fileIO{
public:
    int CountNumbersOfTextLines(const std::string &filePath );
    void input_element(std::string element_file, std::vector<std::vector<int>> &element);
    void input_prism(std::string layer_file, std::vector<std::vector<int>> &layer_pair, std::vector<int> &prism_id);
};

#endif