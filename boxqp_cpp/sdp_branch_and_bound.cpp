#include <thread>
#include "matlab_util.h"
#include "sdp_branch_and_bound.h"
#include "JobQueue.h"
#include "util.h"
#include "config_params.h"
#include "Node.h"
#include "ThreadPool.h"

// root
SDPResult solve_sdp(std::unique_ptr<matlab::engine::MATLABEngine> &matlabPtr,
                    matlab::data::ArrayFactory &factory, InputData *input_data, Config *config) {

    matlab::data::TypedArray<double> Q_matlab = arma_to_matlab_matrix(factory, input_data->Q);
    matlab::data::TypedArray<double> c_matlab = arma_to_matlab_vector(factory, input_data->c);
    // Create StructArray
    std::vector<std::string> f = {"fixing", "fixing_tol", "n_threads", "opt_tol", "sdp_verbose", "sdp_tol",
                                  "cp_maxiter", "cp_tol", "cp_maxineq", "cp_percineq", "cp_epsineq", "cp_activeineq"};

    matlab::data::StructArray struct_matlab = factory.createStructArray({1, 1},f);
    struct_matlab[0]["fixing"] = factory.createScalar<int>(config->branch_and_bound_fixing);
    struct_matlab[0]["fixing_tol"] = factory.createScalar<double>(config->branch_and_bound_fixing_tol);
    struct_matlab[0]["n_threads"] = factory.createScalar<int>(config->matlab_session_threads_root);
    struct_matlab[0]["opt_tol"] = factory.createScalar<double>(config->branch_and_bound_tol);
    struct_matlab[0]["sdp_verbose"] = factory.createScalar<int>(config->sdp_solver_verbose);
    struct_matlab[0]["sdp_tol"] = factory.createScalar<double>(config->sdp_solver_tol);
    struct_matlab[0]["cp_maxiter"] = factory.createScalar<int>(config->cp_max_iter);
    struct_matlab[0]["cp_tol"] = factory.createScalar<double>(config->cp_tol);
    struct_matlab[0]["cp_maxineq"] = factory.createScalar<double>(config->cp_max_ineq);
    struct_matlab[0]["cp_percineq"] = factory.createScalar<double>(config->cp_perc_ineq);
    struct_matlab[0]["cp_epsineq"] = factory.createScalar<double>(config->cp_eps_ineq);
    struct_matlab[0]["cp_activeineq"] = factory.createScalar<double>(config->cp_eps_active);

    std::vector<matlab::data::Array> args({Q_matlab, c_matlab, struct_matlab});

    // Call MATLAB function and return result
    const size_t n_return = 1;
    std::vector<matlab::data::Array> result;
    result = matlabPtr->feval(u"call_solve_boxqp_root",n_return, args);
    matlab::data::StructArray structArray = result[0];
    matlab::data::TypedArray<double> field_best_lb = structArray[0]["best_lb"];
    matlab::data::TypedArray<double> field_best_ub = structArray[0]["best_ub"];
    matlab::data::TypedArray<double> field_best_x_gub = structArray[0]["best_x_ub"];
    matlab::data::TypedArray<double> field_best_x_lb = structArray[0]["best_x_lb"];
    matlab::data::TypedArray<double> field_cp_iter = structArray[0]["cp_iter"];
    matlab::data::TypedArray<double> field_cp_flag = structArray[0]["cp_flag"];
    matlab::data::TypedArray<double> field_ineq_list = structArray[0]["ineq_list"];
    matlab::data::CellArray field_best_B_cell = structArray[0]["best_Bcell"];
    matlab::data::TypedArray<double> field_best_l = structArray[0]["best_l"];
    matlab::data::TypedArray<double> field_idx_i = structArray[0]["idx_i"];
    matlab::data::TypedArray<double> field_n = structArray[0]["n"];
    matlab::data::TypedArray<double> field_n_bin = structArray[0]["N_bin"];
    matlab::data::TypedArray<double> field_init_xfix = structArray[0]["init_xfix"];
    matlab::data::TypedArray<double> field_time_fix = structArray[0]["time_fix"];
    matlab::data::TypedArray<double> field_sdp_fix = structArray[0]["sdp_fix"];
    matlab::data::TypedArray<double> field_n_fix = structArray[0]["n_fixed"];
    matlab::data::TypedArray<double> field_branching_type = structArray[0]["branching_type"];

    double lb = field_best_lb[0];
    double ub = field_best_ub[0];
    arma::vec x_heuristic = matlab_to_arma_vector(field_best_x_gub);
    arma::vec x_sdp = matlab_to_arma_vector(field_best_x_lb);
    arma::vec l_vector = matlab_to_arma_vector(field_best_l);
    int cp_iter = (int) field_cp_iter[0];
    int cp_flag = (int) field_cp_flag[0];
    int n_ineq = (int) field_ineq_list[cp_iter];
    std::vector<std::vector<arma::sp_mat>> B_vector;
    B_vector.reserve(1);
    std::vector<arma::sp_mat> B_inner = matlab_to_arma_vector_sp_mat(field_best_B_cell);
    B_vector.push_back(B_inner);
    int idx_i = (int) field_idx_i[0];
    int branching_index = idx_i;
    int n = (int) field_n[0];
    std::vector<std::pair<int, int>> init_x_fix = matlab_matrix_to_vector_pair(field_init_xfix);
    double time_fix = (double) field_time_fix[0];
    int sdp_fix = (int) field_sdp_fix[0];
    int n_fix = (int) field_n_fix[0];
    int n_bin = (int) field_n_bin[0];
    int branching_type = (int) field_branching_type[0];

    return SDPResult{n, n_bin, init_x_fix, sdp_fix, time_fix, n_fix, lb, ub, x_sdp, x_heuristic, n_ineq, cp_iter, cp_flag, branching_index, branching_type, B_vector, l_vector};
}

// lower bound for children nodes
SDPResult solve_sdp(std::unique_ptr<matlab::engine::MATLABEngine> &matlabPtr, matlab::data::ArrayFactory &factory,
        std::vector<std::vector<arma::sp_mat>> &parent_B_vector, arma::vec &parent_l_vector,
        double global_ub, arma::vec &global_x, std::vector<std::pair<int, int>> &global_x_fix, std::vector<int> &global_x_in,
        InputData *input_data, Config *config) {

    matlab::data::TypedArray<double> Q_matlab = arma_to_matlab_matrix(factory, input_data->Q);
    matlab::data::TypedArray<double> c_matlab = arma_to_matlab_vector(factory, input_data->c);

    // convert data
    matlab::data::TypedArray<double> matlab_x_fix = vector_pair_to_matlab_matrix(factory, global_x_fix);
    arma::vec arma_x_in = arma::conv_to<arma::vec>::from(global_x_in);
    matlab::data::TypedArray<double> matlab_x_in = arma_to_matlab_vector(factory, arma_x_in);
    matlab::data::CellArray matlab_B_cell = factory.createCellArray({0, 0});
    matlab_B_cell = arma_to_matlab_cell(factory, parent_B_vector[0]);
    matlab::data::TypedArray<double> matlab_parent_l = arma_to_matlab_vector(factory, parent_l_vector);
    matlab::data::TypedArray<double> matlab_global_x = arma_to_matlab_vector(factory, global_x);
    matlab::data::TypedArray<double> matlab_gub = factory.createScalar<double>(global_ub);

    // Create StructArray
    std::vector<std::string> f = {"fixing", "fixing_tol", "n_threads", "opt_tol", "sdp_verbose", "sdp_tol",
                                  "cp_maxiter", "cp_tol", "cp_maxineq", "cp_percineq", "cp_epsineq", "cp_activeineq"};

    matlab::data::StructArray struct_matlab = factory.createStructArray({1, 1}, f);
    struct_matlab[0]["fixing"] = factory.createScalar<int>(config->branch_and_bound_fixing);
    struct_matlab[0]["fixing_tol"] = factory.createScalar<double>(config->branch_and_bound_fixing_tol);
    struct_matlab[0]["n_threads"] = factory.createScalar<int>(config->matlab_session_threads_child);
    struct_matlab[0]["opt_tol"] = factory.createScalar<double>(config->branch_and_bound_tol);
    struct_matlab[0]["sdp_verbose"] = factory.createScalar<int>(config->sdp_solver_verbose);
    struct_matlab[0]["sdp_tol"] = factory.createScalar<double>(config->sdp_solver_tol);
    struct_matlab[0]["cp_maxiter"] = factory.createScalar<int>(config->cp_max_iter);
    struct_matlab[0]["cp_tol"] = factory.createScalar<double>(config->cp_tol);
    struct_matlab[0]["cp_maxineq"] = factory.createScalar<double>(config->cp_max_ineq);
    struct_matlab[0]["cp_percineq"] = factory.createScalar<double>(config->cp_perc_ineq);
    struct_matlab[0]["cp_epsineq"] = factory.createScalar<double>(config->cp_eps_ineq);
    struct_matlab[0]["cp_activeineq"] = factory.createScalar<double>(config->cp_eps_active);

    std::vector<matlab::data::Array> args({Q_matlab, c_matlab, matlab_x_fix, matlab_x_in, matlab_gub, matlab_global_x,
                                           matlab_B_cell, matlab_parent_l, struct_matlab});

    // Call MATLAB function and return result
    const size_t n_return = 1;
    std::vector<matlab::data::Array> result;
    result = matlabPtr->feval(u"call_solve_boxqp_child", n_return, args);

    matlab::data::StructArray structArray = result[0];
    matlab::data::TypedArray<double> field_best_lb = structArray[0]["best_lb"];
    matlab::data::TypedArray<double> field_best_ub = structArray[0]["best_ub"];
    matlab::data::TypedArray<double> field_best_x_gub = structArray[0]["best_x_ub"];
    matlab::data::TypedArray<double> field_best_x_lb = structArray[0]["best_x_lb"];
    matlab::data::TypedArray<double> field_cp_iter = structArray[0]["cp_iter"];
    matlab::data::TypedArray<double> field_cp_flag = structArray[0]["cp_flag"];
    matlab::data::TypedArray<double> field_ineq_list = structArray[0]["ineq_list"];
    matlab::data::CellArray field_best_B_cell = structArray[0]["best_Bcell"];
    matlab::data::TypedArray<double> field_best_l = structArray[0]["best_l"];
    matlab::data::TypedArray<double> field_idx_i = structArray[0]["idx_i"];
    matlab::data::TypedArray<double> field_n = structArray[0]["n"];
    matlab::data::TypedArray<double> field_n_bin = structArray[0]["N_bin"];
    matlab::data::TypedArray<double> field_init_xfix = structArray[0]["init_xfix"];
    matlab::data::TypedArray<double> field_time_fix = structArray[0]["time_fix"];
    matlab::data::TypedArray<double> field_sdp_fix = structArray[0]["sdp_fix"];
    matlab::data::TypedArray<double> field_n_fix = structArray[0]["n_fixed"];
    matlab::data::TypedArray<double> field_branching_type = structArray[0]["branching_type"];

    double lb = field_best_lb[0];
    double ub = field_best_ub[0];
    arma::vec x_sdp = matlab_to_arma_vector(field_best_x_lb);
    arma::vec x_heuristic = matlab_to_arma_vector(field_best_x_gub);
    arma::vec l_vector = matlab_to_arma_vector(field_best_l);
    int cp_iter = (int) field_cp_iter[0];
    int cp_flag = (int) field_cp_flag[0];
    int n_ineq = (int) field_ineq_list[cp_iter];
    std::vector<std::vector<arma::sp_mat>> B_vector;
    B_vector.reserve(1);
    std::vector<arma::sp_mat> B_inner = matlab_to_arma_vector_sp_mat(field_best_B_cell);
    B_vector.push_back(B_inner);
    int idx_i = (int) field_idx_i[0];
    int branching_index = idx_i;
    int n = (int) field_n[0];
    std::vector<std::pair<int, int>> init_x_fix = matlab_matrix_to_vector_pair(field_init_xfix);
    double time_fix = (double) field_time_fix[0];
    int sdp_fix = (int) field_sdp_fix[0];
    int n_fix = (int) field_n_fix[0];
    int n_bin = (int) field_n_bin[0];
    int branching_type = (int) field_branching_type[0];

    return SDPResult{n, n_bin, init_x_fix, sdp_fix, time_fix, n_fix, lb, ub, x_sdp, x_heuristic, n_ineq, cp_iter, cp_flag, branching_index, branching_type, B_vector, l_vector};
}



std::vector<JobData *> create_children_jobs(double node_gap, arma::vec &x_sdp, SDPNode *node, int branching_index, int branching_type,
                                            NodeData *parent, SharedData *shared_data, Config *config) {

    if (std::isinf(node->lb) || node_gap <= config->branch_and_bound_tol) {
        delete(node);
        if (parent != nullptr) {
            delete (parent->node);
            delete (parent);
        }
        return std::vector<JobData *>();
    }

    if (branching_index == -1) {

        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        log_file << "PRUNING BY OPTIMALITY " << node->id << "\n";
        if (node->lb - shared_data->global_ub <= -config->branch_and_bound_tol) {
            // update global upper bound
            shared_data->global_ub = node->lb;
            shared_data->global_x = x_sdp;
        }
        delete (node);
        if (parent != nullptr) {
            delete (parent->node);
            delete (parent);
        }
        return std::vector<JobData *>();

        // mutex is automatically released when lock goes out of scope
    }

    std::vector<JobData *> branch;

    auto *child_data_zero = new NodeData();
    child_data_zero->node = new SDPNode(*node);
    child_data_zero->i = branching_index;

    auto *child_data_one = new NodeData();
    child_data_one->node = new SDPNode(*node);
    child_data_one->i = branching_index;

    auto *zero_job_data = new JobData();
    zero_job_data->type = ZERO;
    zero_job_data->node_data = child_data_zero;

    auto *one_job_data = new JobData();
    one_job_data->type = ONE;
    one_job_data->node_data = child_data_one;

    if (branching_type == BINARY) {

        branch.emplace_back(zero_job_data);
        branch.emplace_back(one_job_data);

    } else { // ternary branching

        auto *child_data_in = new NodeData();
        child_data_in->node = new SDPNode(*node);
        child_data_in->i = branching_index;

        auto *in_job_data = new JobData();
        in_job_data->type = IN;
        in_job_data->node_data = child_data_in;

        branch.emplace_back(zero_job_data);
        branch.emplace_back(one_job_data);
        branch.emplace_back(in_job_data);

    }

    if (parent != nullptr) {
        delete (parent->node);
        delete (parent);
    }

    delete (node);

    return branch;

}


std::vector<JobData *> build_child_problem(int type, NodeData *node_data, InputData *input_data, SharedData  *shared_data, Config *config) {

    auto *matlab_struct = new MatlabStruct();
    matlab_struct->matlabPtr = start_matlab(config);

    auto child_node = new SDPNode();
	auto parent = node_data->node;

	double parent_gap = (shared_data->global_ub - parent->lb) / std::abs(shared_data->global_ub);
	if (parent_gap <= config->branch_and_bound_tol)
        return std::vector<JobData *>();

    if (type == ZERO) {
        child_node->x_in = parent->x_in;
        child_node->x_fix = parent->x_fix;
        child_node->x_fix.emplace_back(node_data->i, 0);
    } else if (type == ONE) {
        child_node->x_in = parent->x_in;
        child_node->x_fix = parent->x_fix;
        child_node->x_fix.emplace_back(node_data->i, 1);
    } else if (type == IN) {
        child_node->x_fix = parent->x_fix;
        child_node->x_in = parent->x_in;
        child_node->x_in.emplace_back(node_data->i);
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    SDPResult  sdp_result = solve_sdp(matlab_struct->matlabPtr, matlab_struct->factory,
              parent->B_vector, parent->l_vector,shared_data->global_ub,
              shared_data->global_x, child_node->x_fix, child_node->x_in, input_data, config);

    int n = sdp_result.n;
    int n_bin = sdp_result.n_bin;
    int cp_iter = sdp_result.cp_iter;
    int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    child_node->lb = std::max(sdp_result.lb, parent->lb);
    child_node->ub = sdp_result.ub;
    child_node->B_vector = sdp_result.B_vector;
    child_node->l_vector = sdp_result.l_vector;
    arma::vec x_heuristic = sdp_result.x_heuristic;
    int branching_index = sdp_result.branching_index;
    int branching_type = sdp_result.branching_type;
    arma::vec x_sdp = sdp_result.x_sdp;

    std::vector<std::pair<int, int>> init_x_fix = sdp_result.init_x_fix;
    // add variables after fixing
    for (auto &elem : init_x_fix) {
        child_node->x_fix.emplace_back(elem.first, elem.second);
    }

    int n_sdp_fixing = sdp_result.n_sdp_fixing;
    double time_sdp_fixing = sdp_result.time_sdp_fixing;
    int n_fixed = sdp_result.n_fixed;

    double node_gap;

    {
        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        bool ub_updated = false;
        if (child_node->ub - shared_data->global_ub <= -config->branch_and_bound_tol) {
            // update global upper bound
            shared_data->global_ub = child_node->ub;
            shared_data->global_x = x_heuristic;
            ub_updated = true;
        }

        child_node->id = shared_data->n_nodes;
        shared_data->n_nodes++;
        int open = shared_data->queue->getSize();

        node_gap = (shared_data->global_ub - child_node->lb) / std::abs(shared_data->global_ub);
        double gap = node_gap;
        Node *min_lb_node = shared_data->queue->getMinLb();
        if (min_lb_node != nullptr)
            gap = (shared_data->global_ub - min_lb_node->lb) / std::abs(shared_data->global_ub);

        shared_data->gap = gap;

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        auto time = (double) duration.count();

        print_log_sdp(log_file, n, n_bin, parent->id, child_node->id,
                      parent->lb, child_node->lb, time, cp_iter, cp_flag, n_ineq, n_sdp_fixing, time_sdp_fixing, n_fixed,
                      child_node->ub,shared_data->global_ub, node_data->i,
                      node_gap,shared_data->gap, open, ub_updated);

    }

    delete(matlab_struct);

    return create_children_jobs(node_gap, x_sdp, child_node, branching_index, branching_type, node_data, shared_data, config);
}


std::vector<JobData *> build_root_problem(MatlabStruct *matlab_struct, InputData *input_data, SharedData *shared_data, Config *config) {

    // init root
    SDPNode *root;
    root = new SDPNode();
    root->id = shared_data->n_nodes;

    auto start_time = std::chrono::high_resolution_clock::now();

    SDPResult sdp_result = solve_sdp(matlab_struct->matlabPtr, matlab_struct->factory, input_data, config);

    int n = sdp_result.n;
    int n_bin = sdp_result.n_bin;
    int cp_iter = sdp_result.cp_iter;
    int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    root->lb = sdp_result.lb;
    root->ub = sdp_result.ub;
    root->B_vector = sdp_result.B_vector;
    root->l_vector = sdp_result.l_vector;
    arma::vec x_heuristic = sdp_result.x_heuristic;
    int branching_index = sdp_result.branching_index;
    int branching_type = sdp_result.branching_type;
    arma::vec x_sdp = sdp_result.x_sdp;
    root->x_fix = sdp_result.init_x_fix;

    int n_sdp_fixing = sdp_result.n_sdp_fixing;
    double time_sdp_fixing = sdp_result.time_sdp_fixing;
    int n_fixed = sdp_result.n_fixed;

    // update shared data
    shared_data->global_ub = root->ub;
    shared_data->global_x = x_heuristic;
    shared_data->n_nodes++;

    int open = shared_data->queue->getSize();

    double node_gap = (shared_data->global_ub - root->lb) / std::abs(shared_data->global_ub);
    shared_data->gap = node_gap;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    auto time = (double) duration.count();

    print_log_sdp(log_file, n, n_bin, -1, root->id, -std::numeric_limits<double>::infinity(),
                  root->lb, time, cp_iter, cp_flag, n_ineq, n_sdp_fixing, time_sdp_fixing, n_fixed,
                  root->ub, shared_data->global_ub, -1, node_gap, node_gap, open, true);

    //return std::make_pair(nullptr, nullptr);
    return create_children_jobs(node_gap, x_sdp, root, branching_index, branching_type, nullptr, shared_data, config);
}

bool is_thread_pool_working(std::vector<bool> &thread_state) {
    int count = 0;
    for (auto && i : thread_state) {
        if (i)
            count++;
    }
    if (count == 0)
        return false;
    return true;
}

void save_x_to_file(Config *config, arma::vec &x){

    std::ofstream f;
    f.open(config->result_path);
    for (size_t i = 0; i < x.size(); i++){
        double val = x(i);
        f << val << "\n";
    }
    f.close();
}


arma::vec sdp_branch_and_bound(Config *config, InputData *input_data) {

    int n_thread = config->branch_and_bound_parallel;

    JobAbstractQueue *queue;
    switch (config->branch_and_bound_visiting_strategy) {
        case DEPTH_FIRST:
            queue = new JobStack();
            break;
        case BEST_FIRST:
            queue = new JobPriorityQueue();
            break;
        case BREADTH_FIRST:
            queue = new JobQueue();
            break;
        default:
            queue = new JobPriorityQueue();
    }

    auto *shared_data = new SharedData();
    shared_data->global_ub = std::numeric_limits<double>::infinity();
    shared_data->n_nodes = 0;
    shared_data->queue = queue;

    shared_data->threadStates.reserve(n_thread);
    for (int i = 0; i < n_thread; i++) {
        shared_data->threadStates.push_back(false);
    }

    ThreadPool pool(config, shared_data, input_data, n_thread);
    
    print_header_sdp(log_file);

    auto start_all = std::chrono::high_resolution_clock::now();
    
    auto *matlab_struct = new MatlabStruct();
    matlab_struct->matlabPtr = start_matlab(config);

    std::vector<JobData *> jobs = build_root_problem(matlab_struct, input_data, shared_data, config);

    delete (matlab_struct);
    
    double root_gap = shared_data->gap;

    // submit jobs to the thread pool
    for (auto & job : jobs) {
        pool.addJob(job);
    }

    while (true) {

        {
            std::unique_lock<std::mutex> l(shared_data->queueMutex);
            while (is_thread_pool_working(shared_data->threadStates) && shared_data->n_nodes < config->branch_and_bound_max_nodes) {
                shared_data->mainConditionVariable.wait(l);
            }

            if (shared_data->queue->empty() || shared_data->n_nodes >= config->branch_and_bound_max_nodes)
                break;
        }

    }

    auto end_all = std::chrono::high_resolution_clock::now();
    auto duration_all = std::chrono::duration_cast<std::chrono::seconds>(end_all - start_all);

    pool.quitPool();

    if (queue->empty()) {
        shared_data->gap = 0.0;
    }

    log_file << "\n";
    log_file << "TIME: " << duration_all.count() << " sec\n";
    log_file << "NODES: " << shared_data->n_nodes << "\n";
    log_file << "ROOT_GAP: " << std::max(0.0, root_gap) << "\n";
    log_file << "GAP: " << std::max(0.0, shared_data->gap) << "\n";
    log_file << "OPT: " << shared_data->global_ub << "\n\n";

    arma::vec result = shared_data->global_x;
    save_x_to_file(config, result);

    // free memory

    delete (input_data);
    delete (queue);
    delete (shared_data);

    return result;

}
