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
    void export_vtu(const std::string &file, std::vector<std::vector<double>> x, std::vector<std::vector<int>> element, std::vector<double> element_wss_u,std::vector<double> element_wss_v, std::vector<double> element_wss_w);
    void export_vtu_OSI(const std::string &file, std::vector<std::vector<double>> x, std::vector<std::vector<int>> element, std::vector<double> OSI);


};

#endif