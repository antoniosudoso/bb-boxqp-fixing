#include "util.h"
#include <iomanip>

void save_x_to_file(Config *config, arma::vec &x) {

    std::ofstream f;
    f.open(config->result_path);
    for (size_t i = 0; i < x.size(); i++){
        double val = x(i);
        f << val << "\n";
    }
    f.close();
}

void print_header_sdp(std::stringstream &log_file) {

    log_file << "\n" << "|" <<
             std::setw(5) << "N" << "|" <<
             std::setw(5) << "BIN" << "|" <<
             std::setw(8) << "PARENT" << "|" <<
             std::setw(8) << "NODE" << "|" <<
             std::setw(12) << "LB_PAR" << "|" <<
             std::setw(12) << "LB" << "|" <<
             std::setw(10) << "TIME (s)" << "|" <<
             std::setw(8) << "CP_ITER" << "|" <<
             std::setw(8) << "CP_FLAG" << "|" <<
             std::setw(10) << "CP_INEQ" << "|" <<
             std::setw(8) << "SDP_FIX" << "|" <<
             std::setw(9) << "TIME_FIX" << "|" <<
             std::setw(8) << "N_FIX" << "|" <<
             std::setw(12) << "UB" << "|" <<
             std::setw(12) << "GUB" << "|" <<
             std::setw(6) << "I_FIX" << "|" <<
             std::setw(13) << "NODE_GAP" << "|" <<
             std::setw(13) << "GAP" << "|" <<
             std::setw(6) << "OPEN" << "|"
             << std::endl;

}

void print_log_sdp(std::stringstream &log_file, int n, int n_bin, int node_parent, int node, double lb_parent, double lb,
                   double time, int cp_iter, int cp_flag, int n_ineq, int n_sdp_fixing, double time_sdp_fixing, int n_fixed,
                   double ub, double gub, int i, double node_gap, double gap, int open, bool update) {

    if (!update) {

        log_file << "|" <<
                 std::setw(5) << n << "|" <<
                 std::setw(5) << n_bin << "|" <<
                 std::setw(8) << node_parent << "|" <<
                 std::setw(8) << node << "|" <<
                 std::setw(12) << lb_parent << "|" <<
                 std::setw(12) << lb << "|" <<
                 std::setw(10) << time << "|" <<
                 std::setw(8) << cp_iter << "|" <<
                 std::setw(8) << cp_flag << "|" <<
                 std::setw(10) << n_ineq << "|" <<
                 std::setw(8) << n_sdp_fixing << "|" <<
                 std::setw(9) << time_sdp_fixing << "|" <<
                 std::setw(8) << n_fixed << "|" <<
                 std::setw(12) << ub << "|" <<
                 std::setw(12) << gub << "|" <<
                 std::setw(6) << i << "|" <<
                 std::setw(13) << node_gap << "|" <<
                 std::setw(13) << gap << "|" <<
                 std::setw(6) << open << "|"
                 << std::endl;

    } else {

        log_file << "|" <<
                 std::setw(5) << n << "|" <<
                 std::setw(5) << n_bin << "|" <<
                 std::setw(8) << node_parent << "|" <<
                 std::setw(8) << node << "|" <<
                 std::setw(12) << lb_parent << "|" <<
                 std::setw(12) << lb << "|" <<
                 std::setw(10) << time << "|" <<
                 std::setw(8) << cp_iter << "|" <<
                 std::setw(8) << cp_flag << "|" <<
                 std::setw(10) << n_ineq << "|" <<
                 std::setw(8) << n_sdp_fixing << "|" <<
                 std::setw(9) << time_sdp_fixing << "|" <<
                 std::setw(8) << n_fixed << "|" <<
                 std::setw(12) << ub << "|" <<
                 std::setw(11) << gub << "*|" <<
                 std::setw(6) << i << "|" <<
                 std::setw(13) << node_gap << "|" <<
                 std::setw(13) << gap << "|" <<
                 std::setw(6) << open << "|"
                 << std::endl;

    }
}