## Fix and Bound: An efficient approach for solving large-scale quadratic programming problems with box constraints

This repository contains the source code of the B&B algorithm described in the paper ["Fix and Bound: An efficient approach for solving large-scale quadratic programming problems with box constraints"](https://arxiv.org/abs/2211.08911) for solving nonconvex BoxQPs to global optimality:

$$ 
\begin{align*}
\min_{} \quad & \frac{1}{2} {\bf x }^\top {\bf Q} {\bf x} + {\bf c}^\top {\bf x}  \\
\textrm{s.t.} \quad & {\bf x}\in [0, 1]^n,\\
\end{align*}
$$

where ${\bf Q}\in \mathbb{R}^{n\times n}$ is symmetric and not positive semidefinite and ${\bf c}\in \mathbb{R}^n$.

## Installation
The B&B algorithm is implemented in C++/MATLAB and requires the following solvers: 
- [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/) 
- [Gurobi](https://www.gurobi.com/)
- [SNOPT](https://ccom.ucsd.edu/~optimizers/solvers/snopt/)

SDPNAL+ is called by using the [MATLAB Engine API](https://www.mathworks.com/help/matlab/calling-matlab-engine-from-cpp-programs.html) for C++, whereas Gurobi and SNOPT are called from MATLAB functions. 

Ubuntu and Debian instructions:

1) Install CMake and Armadillo:
 ```
sudo apt-get update
sudo apt-get install cmake libarmadillo-dev
```
4) Open the makefile `boxqp_cpp/Makefile` 
	- Set the variable `matlab_path` to the MATLAB installation folder.
5) Compile the code:

```
make
```
This code has been tested on Ubuntu 20.04 LTS with MATLAB R2020b Update 7, SDPNAL+ 1.0, Gurobi 9.5.1 and SNOPT 7.7. To speed up the B&B search, we implemented a configurable pool of POSIX threads.

## Configuration
Various parameters used in the B&B algorithm can be modified via the configuration file `boxqp_cpp/config.txt`:

- `BRANCH_AND_BOUND_TOL` - optimality tolerance of the algorithm
- `BRANCH_AND_BOUND_PARALLEL` -  thread pool size: single thread (1), multi-thread (> 1)
- `BRANCH_AND_BOUND_MAX_NODES` - maximum number of nodes
- `BRANCH_AND_BOUND_VISITING_STRATEGY` - best first (0),  depth first (1), breadth first (2)
- `BRANCH_AND_BOUND_FIXING` - do not fix variables (0),  fix variables to 0 or 1 (1)
- `BRANCH_AND_BOUND_FIXING_TOL` - tolerance for fixing variables (only if fixing is enabled)
- `MATLAB_SESSION_THREADS_ROOT` - number of threads for the MATLAB session at the root
- `MATLAB_SESSION_THREADS_CHILD` - number of threads for the MATLAB session of children nodes
- `MATLAB_SOURCE_FOLDER` - MATLAB folder containing source files of the algorithm
- `SDP_SOLVER_TOL` - accuracy of SDPNAL+ in the relative KKT residual
- `SDP_SOLVER_VERBOSE` - do not display SDP log (0), display SDP log (1)
- `GUROBI_FOLDER` - Gurobi installation folder
- `SNOPT_FOLDER` - SNOPT installation folder
- `SNOPT_LICENSE` - SNOPT license file (.lic)
- `SDPNAL_FOLDER` - SDPNAL+ installation folder
- `CP_MAX_ITER` - maximum number of cutting-plane iterations
- `CP_TOL` - tolerance between two consecutive cutting-plane iterations
- `CP_MAX_INEQ` - maximum number of triangle inequalities to separate
- `CP_PERC_INEQ` - fraction of the most violated triangle inequalities to add
- `CP_EPS_INEQ` - tolerance for checking the violation of triangle inequalities
- `CP_EPS_ACTIVE` - tolerance for detecting active  triangle inequalities

## Usage
```
./bb <config_file> <data_file> <log_file> <sol_file>
```
- `config_file` - configuration file
- `data_file` - problem data ${\bf Q}$, ${\bf c}$
- `log_file` - log file
- `sol_file` - optimal solution ${\bf x^\star}$

File `data_file` contains the problem data in the following format:

```
n
c_1 c_2 ... c_n
Q_11 Q_12 ... Q_1n
Q_21 Q_22 ... Q_2n
...
Q_n1 Q_n2 ... Q_nn
```

## Log

The log file reports the progress of the algorithm at the node level:

- `N` - node size
- `BIN` - number of variables in set $N$
- `PARENT` - parent node ID
- `NODE` - node ID
- `LB_PAR` - lower bound of the parent node
- `LB` - lower bound
- `TIME (s)` - running time in seconds
- `CP_ITER` - number of cutting-plane iterations
- `CP_FLAG` - termination flag of the cutting-plane procedure
    - `-1` - maximum number of iterations
    -  `0` - node must be pruned
    -  `1` - less than $10n$ triangle inequalities
    -  `2` - lower bound not significantly improved
    -  `3` - lower bound not improved
- `CP_INEQ` - number of inequalities added in the last cutting-plane iteration
- `SDP_FIX` - number of SDPs solved for fixing variables
- `TIME_FIX` - time spent in seconds for fixing variables
- `N_FIX` - number of fixed variables
- `UB` - current upper bound
- `GUB` - global upper bound
- `I_FIX` - current branching decision
- `NODE_GAP` - gap at the node
- `GAP` - overall gap 
- `OPEN` - number of open nodes

At the end of the file, it shows the statistics:

- `TIME` - overall computational time in seconds
- `NODES` - overall number of nodes 
- `ROOT_GAP` - gap at the root node
- `GAP` - overall optimality gap
- `OPT` - optimal objective function value

Log file example:

```
DATA_PATH: /boxqp_instances/spar200-075-2.in
LOG_PATH: /log_spar200-075-2.txt

BRANCH_AND_BOUND_TOL: 0.0001
BRANCH_AND_BOUND_PARALLEL: 2
BRANCH_AND_BOUND_MAX_NODES: 1000
BRANCH_AND_BOUND_VISITING_STRATEGY: 0
BRANCH_AND_BOUND_FIXING: 1
BRANCH_AND_BOUND_FIXING_TOL: 0.01
MATLAB_SESSION_THREADS_ROOT: 4
MATLAB_SESSION_THREADS_CHILD: 2
MATLAB_SOURCE_FOLDER: /bb-boxqp-fixing/boxqp_matlab/
SDP_SOLVER_TOL: 0.0001
SDP_SOLVER_VERBOSE: 0
CP_MAX_ITER: 10
CP_TOL: 0.0001
CP_MAX_INEQ: 100000
CP_PERC_INEQ: 0.1
CP_EPS_INEQ: 0.0001
CP_EPS_ACTIVE: 0.00
GUROBI_FOLDER: /gurobi951/
SNOPT_FOLDER: /snopt7_matlab/
SNOPT_LICENSE: /snopt7.lic
SDPNAL_FOLDER: /SDPNAL+/

|    N|  BIN|  PARENT|    NODE|      LB_PAR|          LB|  TIME (s)| CP_ITER| CP_FLAG|   CP_INEQ| SDP_FIX| TIME_FIX|   N_FIX|          UB|         GUB| I_FIX|     NODE_GAP|          GAP|  OPEN|
|  200|  126|      -1|       0|        -inf|    -22266.4|       364|       3|       1|     76669|       1|  58.1208|      19|      -22163|     -22163*|    -1|   0.00466673|   0.00466673|     0|
|  180|  106|       0|       1|    -22266.4|    -22223.6|       110|       1|       1|     59120|       0|        0|       0|      -22062|      -22163|    38|   0.00273636|   0.00273636|     0|
|  180|  106|       0|       2|    -22266.4|    -22210.9|       240|       2|       1|     60048|       1|  40.9975|      10|      -22163|      -22163|    38|   0.00215997|   0.00273636|     1|
|  179|  105|       1|       3|    -22223.6|    -22183.4|       156|       1|       1|     59324|       1|  38.8991|       2|      -22062|      -22163|   199|  0.000921238|   0.00215997|     2|
|  179|  105|       1|       4|    -22223.6|    -22187.3|       161|       1|       1|     59630|       1|  41.3265|       3|    -22099.4|      -22163|   199|    0.0010967|   0.00215997|     3|
|  169|   95|       2|       5|    -22210.9|    -22172.6|       188|       1|       1|     58267|       1|  60.8439|      17|      -22163|      -22163|    81|  0.000432555|    0.0010967|     4|
|  175|  101|       4|       6|    -22187.3|    -22153.1|        62|       0|       0|     56624|       0|        0|       0|    -22099.4|      -22163|    97| -0.000447315|    0.0010967|     5|
|  169|   95|       2|       7|    -22210.9|    -22171.6|       156|       1|       1|     55533|       1|  37.3501|       7|      -22136|      -22163|    81|  0.000385916|  0.000921238|     4|
|  175|  101|       4|       8|    -22187.3|      -22163|        57|       0|       0|     56624|       0|        0|       0|    -22095.4|      -22163|    97|  9.47342e-08|  0.000921238|     5|
|  176|  102|       3|       9|    -22183.4|    -22156.2|        60|       0|       0|     57021|       0|        0|       0|    -22056.5|      -22163|   119|  -0.00030788|  0.000432555|     4|
|  176|  102|       3|      10|    -22183.4|    -22157.1|        67|       0|       0|     57021|       0|        0|       0|    -22029.7|      -22163|   119|   -0.0002681|  0.000432555|     3|
|  151|   77|       5|      11|    -22172.6|    -22163.1|        58|       0|       0|     45312|       0|        0|       0|      -22163|      -22163|    18|  3.10148e-06|  0.000385916|     2|
|  151|   77|       5|      12|    -22172.6|    -22139.5|        55|       0|       0|     45312|       0|        0|       0|      -22138|      -22163|    18|  -0.00106123|  0.000385916|     1|
|  161|   87|       7|      13|    -22171.6|    -22145.9|        63|       0|       0|     49688|       0|        0|       0|    -22124.8|      -22163|    18| -0.000771785| -0.000771785|     0|
|  161|   87|       7|      14|    -22171.6|    -22149.5|        52|       0|       0|     49688|       0|        0|       0|      -22136|      -22163|    18| -0.000609419| -0.000609419|     0|

TIME: 1145 sec
NODES: 15
ROOT_GAP: 0.00466673
GAP: 0
OPT: -22163
```
