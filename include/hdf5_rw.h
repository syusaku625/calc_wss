#ifndef _HDF5RW_H_
#define _HDF5RW_H_

#include<vector>
#include<string>
#include <H5Cpp.h>

class HDF5RW{
    public:
        std::vector<double> readH5_double_1D(std::string h5file, const int number, std::string data_n);
        std::vector<double> readH5_double_2D(std::string h5file, const int number, std::string data_n, int &nx, int &ny);
        std::vector<double> importHDF5_double_H5_1D(H5::H5File &file, const std::string &dataName);
        std::vector<double> importHDF5_double_H5_2D(H5::H5File &file, const std::string &dataName, int &nx_t, int &ny_t);
        void exportHDF5_double_1D(H5::H5File &file, const std::string &dataName, std::vector<double> i_data, int i_dim);
};

#endif