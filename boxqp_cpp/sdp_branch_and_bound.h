#ifndef CLUSTERING_SDP_BRANCH_AND_BOUND_H
#define CLUSTERING_SDP_BRANCH_AND_BOUND_H

#include <armadillo>
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#include "JobQueue.h"
#include "config_params.h"

typedef struct MatlabStruct {

    std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr;
    matlab::data::ArrayFactory factory;

} MatlabStruct;


typedef struct SharedData {

    // Between workers and main
    std::condition_variable mainConditionVariable;
    std::vector<bool> threadStates;

    // Queue of requests waiting to be processed
    JobAbstractQueue *queue;
    // This condition variable is used for the threads to wait until there is work to do
    std::condition_variable queueConditionVariable;
    // Mutex to protect queue
    std::mutex queueMutex;

    double global_ub;
    arma::vec global_x;
    double gap;
    int n_nodes;

    // for root
    int root_n_sdp_fixing;
    double root_time_sdp_fixing;

} SharedData;

typedef struct InputData {

    arma::mat Q;
    arma::vec c;

} InputData;


arma::vec sdp_branch_and_bound(Config *config, InputData *input_data);
std::vector<JobData *> build_root_problem(MatlabStruct *matlab_struct, InputData *input_data, SharedData *shared_data, Config *config);
std::vector<JobData *> build_child_problem(int type, NodeData *job_data, InputData *input_data, SharedData *shared_data, Config *config);
#endif //CLUSTERING_SDP_BRANCH_AND_BOUND_H
