#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "MPC.h"

using CppAD::AD;

// N: Number of the timestamps
// dt: How much time elapses between actuations
// N x dt = Durations of over which predictions are made
const size_t N = 20;
const double dt = 0.1;

// Reference verocity
// miles/h => meters/s
const double ref_v = 40 * MPC::TO_METERS_PER_SECOND;

const size_t x_start = 0;
const size_t y_start = x_start + N;
const size_t psi_start = y_start + N;
const size_t v_start = psi_start + N;
const size_t cte_start = v_start + N;
const size_t epsi_start = cte_start + N;
const size_t delta_start = epsi_start + N;
const size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Constructor
  FG_eval(Eigen::VectorXd coeffs) {
    coeffs_ = coeffs;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  // fg: a vector of the cost constraints,
  // vars: a vector of variable values (state & actuators)
  void operator()(ADvector& fg, const ADvector& vars) {

    // The first element of fg is used to store the cost
    fg[0] = 0;

    // Calculate the cost from here
    // To minimize the cross-track error, heading error, velocity error

    const double CTE_WEIGHT = 1.0;
    const double EPSI_WEIGHT = 100;
    const double VELOCITY_WEIGHT = 1.0;

    for (int t = 0; t < N; t++) {
      fg[0] += CTE_WEIGHT * CppAD::pow(vars[cte_start + t], 2);
      fg[0] += EPSI_WEIGHT * CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += VELOCITY_WEIGHT * CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    // To minimize the use of actuators for smooth turns
    // not to make the vehicle velocity change too radically

    const double STEERING_WEIGHT = 500;
    const double THROTTLE_WEIGHT = 1.0;

    for (int t = 0; t < N - 1; t++) {
      fg[0] += STEERING_WEIGHT * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += THROTTLE_WEIGHT * CppAD::pow(vars[a_start + t], 2);
    }

    // To minimize the value gap between sequential actuations
    // to make control decisions more consistent or smoother

    const double STEERING_DIFF_WEIGHT = 1000;
    const double THROTTLE_DIFF_WEIGHT = 1.0;

    for (int t = 0; t < N - 2; t++) {
      fg[0] += STEERING_DIFF_WEIGHT * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += THROTTLE_DIFF_WEIGHT * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    // Store the values of the states and actuators
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    for (int t = 1; t < N; t++) {
      // The state at time t+1
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // The state at time t
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      // Use 3rd degree polynomial
      AD<double> f0 = coeffs_[0] + coeffs_[1] * x0 + coeffs_[2] * CppAD::pow(x0, 2) + coeffs_[3] * CppAD::pow(x0, 3);
      AD<double> psides0 = CppAD::atan(coeffs_[1] + 2 * coeffs_[2] * x0 + 3 * coeffs_[3] * CppAD::pow(x0, 2));

      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / MPC::LF * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / MPC::LF * dt);
    }
  }

private:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs_;
};

//
// MPC class definition implementation.
//

MPC::MPC() {}
MPC::~MPC() {}

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double MPC::LF = 2.67;

const double MPC::TO_METERS_PER_SECOND = 0.44704;

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  typedef CPPAD_TESTVECTOR(double) Dvector;
  bool ok = true;

  // The number '2' comes from the fact that this model is using 2 actuators
  int state_size = state.size();
  size_t n_vars = N * state_size + (N - 1) * 2;
  size_t n_constraints = N * state_size;

  // Initial value of the independent variables.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lowerlimits to the max negative and positive values.
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25 degrees (values in radians).
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // Object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  // Options for IPOPT solver
  std::string options;

  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // Minimize the cost here finally by solving the non-linear programming
  CppAD::ipopt::solve_result<Dvector> solution;
  CppAD::ipopt::solve<Dvector, FG_eval>(
                                        options,
                                        vars,
                                        vars_lowerbound,
                                        vars_upperbound,
                                        constraints_lowerbound,
                                        constraints_upperbound,
                                        fg_eval,
                                        solution
                                        );

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  predicted_x_.clear();
  predicted_y_.clear();

  for (int i = 0; i < N-1; i++) {
    predicted_x_.push_back(solution.x[x_start + i + 1]);
    predicted_y_.push_back(solution.x[y_start + i + 1]);
  }

  return {
    solution.x[delta_start],
    solution.x[a_start]
  };
}

vector<double> MPC::GetPredictedX() {
  return predicted_x_;
}

vector<double> MPC::GetPredictedY() {
  return predicted_y_;
}
