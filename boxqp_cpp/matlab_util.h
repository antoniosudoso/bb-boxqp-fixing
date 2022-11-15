#ifndef MATLAB_CALLER_MATLAB_UTIL_H
#define MATLAB_CALLER_MATLAB_UTIL_H


#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#include "config_params.h"
#include <armadillo>

std::unique_ptr<matlab::engine::MATLABEngine> start_matlab(Config *config);
matlab::data::TypedArray<double> vector_pair_to_matlab_matrix(matlab::data::ArrayFactory &factory, std::vector<std::pair<int, int>> &pairs);
matlab::data::TypedArray<double> arma_to_matlab_vector(matlab::data::ArrayFactory &factory, arma::vec &v);
matlab::data::TypedArray<double> arma_to_matlab_matrix(matlab::data::ArrayFactory &factory, arma::mat &X);
matlab::data::SparseArray<double> arma_to_matlab_sparse(matlab::data::ArrayFactory &factory, arma::sp_mat &X);
matlab::data::CellArray arma_to_matlab_cell(matlab::data::ArrayFactory &factory, std::vector<arma::sp_mat> &B_vector);
matlab::data::CellArray arma_to_matlab_vector_cell(matlab::data::ArrayFactory &factory, std::vector<std::vector<arma::sp_mat>> &B_vector);
arma::mat matlab_to_arma_matrix(matlab::data::TypedArray<double> &X_matlab);
arma::vec matlab_to_arma_vector(matlab::data::TypedArray<double> &v_matlab);
std::vector<arma::sp_mat> matlab_to_arma_vector_sp_mat(matlab::data::CellArray &B_matlab);
std::vector<std::vector<arma::sp_mat>> matlab_to_arma_B_cell(matlab::data::CellArray &B_matlab);
arma::sp_mat matlab_to_arma_sparse(matlab::data::SparseArray<double> &X_matlab);
std::vector<std::pair<int, int>> matlab_matrix_to_vector_pair(matlab::data::TypedArray<double> &X_matlab);


#endif //MATLAB_CALLER_MATLAB_UTIL_H
