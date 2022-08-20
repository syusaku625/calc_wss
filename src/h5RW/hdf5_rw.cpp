#include "hdf5_rw.h"

using namespace std;

vector<double> HDF5RW::importHDF5_double_H5_2D(H5::H5File &file, const std::string &dataName, int &nx_t, int &ny_t) {
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

vector<double> HDF5RW::importHDF5_double_H5_1D(H5::H5File &file, const std::string &dataName) {
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

vector<double> HDF5RW::readH5_double_2D(std::string h5file, const int number, string data_n, int &nx, int &ny) {
  string dataName;
  H5::H5File file(h5file, H5F_ACC_RDONLY);

  dataName = to_string(number)+"/" + data_n;
  vector<double> d = importHDF5_double_H5_2D(file, dataName, nx, ny);
  return d;
}

vector<double> HDF5RW::readH5_double_1D(std::string h5file, const int number, string data_n) {
  string dataName;
  H5::H5File file(h5file, H5F_ACC_RDONLY);

  dataName = to_string(number)+"/" + data_n;
  vector<double> d=importHDF5_double_H5_1D(file, dataName);
  return d;
}

void HDF5RW::exportHDF5_double_1D(H5::H5File &file, const std::string &dataName, vector<double> i_data, int i_dim) {
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