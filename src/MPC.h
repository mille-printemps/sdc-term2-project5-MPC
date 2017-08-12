#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();
  virtual ~MPC();

  // Coefficient to convert Miles Per Hour to Meter Per Second
  static const double TO_METERS_PER_SECOND;
  static const double LF;

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

  // Get the ptedicted x and y of the positions
  vector<double> GetPredictedX();
  vector<double> GetPredictedY();

private:
  vector<double> predicted_x_;
  vector<double> predicted_y_;
};

#endif /* MPC_H */
