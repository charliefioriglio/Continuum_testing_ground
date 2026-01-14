#pragma once

#include <vector>
#include <cmath>
#include <print>

class rk45_dormand_prince {
public:
  rk45_dormand_prince(double tolerance_abs, double tolerance_rel) {
    // constructor, initialize the data members
    atol = tolerance_abs;
    rtol = tolerance_rel;
  }

  template <typename F,  typename StopCondition>
  std::vector<double> integrate(const F &f, double x0, const StopCondition &stop_condition,
                                const std::vector<double> &y0,
                                bool dense_output = false,
                                const std::vector<double> &x_out = {}) {
    // resize the temporary arrays according to the given y0
    n_eq = y0.size();
    y_tmp.resize(n_eq);
    y_err.resize(n_eq);

    k1.resize(n_eq);
    k2.resize(n_eq);
    k3.resize(n_eq);
    k4.resize(n_eq);
    k5.resize(n_eq);
    k6.resize(n_eq);
    k7.resize(n_eq);

    r1.resize(n_eq);
    r2.resize(n_eq);

    // clear the arrays so that we start fresh every time we call this function
    hs.clear();
    xs.clear();
    result.clear();

    // initial conditions
    double h = 0.1; // arbitrary initial step size
    double x = x0;
    std::vector<double> y = y0;
    double err = 1.0;
    // err_prev is for you to implement the more advanced step size control
    double err_prev = 1.0;
    bool step_rejected = false;

    auto y_prev = y; // This is used for computing dense output coefficients
    int xout_idx = 0; // This is the current index for the x_out array

    // Handle initial condition for dense output
    if (dense_output && !x_out.empty()) {
      // Add initial condition if x0 is in the output array
      while (xout_idx < x_out.size() && x_out[xout_idx] <= x0) {
        if (std::abs(x_out[xout_idx] - x0) < 1e-12) {
          xs.push_back(x0);
          result.push_back(y0);
        }
        xout_idx++;
      }
    } else if (!dense_output) {  
      xs.push_back(x);
      result.push_back(y);
    }

    while (!stop_condition(x, y)) {
      y = step(f, h, x, y);
      // Compute the rescaled error
      err = error(y);
      // If err is fortuitously too small, set it to some lower bound
      err = std::max(err, 1.0e-12);

      // If the error is below 1, we accept the step by pushing the new x, y,
      // and h to the arrays, and advancing x by the current value of h
      if (err < 1.0) {
        if (dense_output) {
          // Ensure we have output points to process
          if (!x_out.empty()) {
            compute_dense_rs(y, y_prev, h);
            
            // Find all x_out values in the current step interval [x, x + h]
            while (xout_idx < x_out.size() && x_out[xout_idx] <= x + h) {
              if (x_out[xout_idx] >= x) {
                // Interpolate y value at x_out[xout_idx]
                double theta = (x_out[xout_idx] - x) / h; // normalized position in step
                theta = std::max(0.0, std::min(1.0, theta)); // Clamp to [0,1] for safety
                
                std::vector<double> y_interp(n_eq);
                for (int i = 0; i < n_eq; i++) {
                  // 4th order Dormand-Prince dense output formula
                  // y(x_0 + θh) = y_0 + θ * r1 + θ² * r2
                  y_interp[i] = y_prev[i] + theta * r1[i] + theta * theta * r2[i];
                }
                
                xs.push_back(x_out[xout_idx]);
                result.push_back(y_interp);
              }
              xout_idx++;
            }
          }
        } else {
          // This is what we do when not using dense output
          hs.push_back(h);
          xs.push_back(x + h);
          result.push_back(y);
        }
        x += h;
        y_prev = y;
        step_rejected = false;
      } else {
        step_rejected = true;
      }

      // More robust step size control with better stability
      double S = 0.84; // Slightly more conservative safety factor
      double h_new = h * S * std::pow(err, -0.2) * std::pow(err_prev, 0.1);
      
      // Don't allow step size to increase after rejection
      if (step_rejected && h_new > h) {
        h_new = h;
      }
      
      // More conservative step size limits
      h = std::max(hmin, std::min(h_new, hmax));
      err_prev = err;

      // Uncomment the following line to see the step size and error at each step
      // std::println("x = {}, h = {}, err = {}", x, h, err);
    }
    return y;
  }

  template <typename F>
  std::vector<double> step(const F& f, double h, double x, const std::vector<double> &y) {
    // Compute the next step in y, given x and y of the current step
    std::vector<double> y_next(n_eq);

    // Compute the k coefficients for Dormand-Prince method
    
    // k1 = f(x, y)
    k1 = f(x, y);
    
    // k2 = f(x + c2*h, y + h*(a21*k1))
    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + h * a21 * k1[i];
    }
    k2 = f(x + c2 * h, y_tmp);
    
    // k3 = f(x + c3*h, y + h*(a31*k1 + a32*k2))
    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
    }
    k3 = f(x + c3 * h, y_tmp);
    
    // k4 = f(x + c4*h, y + h*(a41*k1 + a42*k2 + a43*k3))
    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    }
    k4 = f(x + c4 * h, y_tmp);
    
    // k5 = f(x + c5*h, y + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4))
    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    }
    k5 = f(x + c5 * h, y_tmp);
    
    // k6 = f(x + h, y + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5))
    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    }
    k6 = f(x + h, y_tmp);
    
    // Compute 5th order solution: y_next = y + h*(a71*k1 + a72*k2 + a73*k3 + a74*k4 + a75*k5 + a76*k6)
    for (int i = 0; i < n_eq; i++) {
      y_next[i] = y[i] + h * (a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
    }
    
    // k7 = f(x + h, y_next)
    k7 = f(x + h, y_next);
    
    // Compute error estimate as difference between 5th and 4th order solutions
    for (int i = 0; i < n_eq; i++) {
      y_err[i] = h * (e1 * k1[i] + e2 * k2[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
    }

    return y_next;
  }

  double error(const std::vector<double> &y) {
    // Compute the rescaled scalar error from the y_err vector computed in
    // the step function, and return it.
    double err = 0.0;
    for (int i = 0; i < n_eq; i++) {
      double tol = atol + rtol * std::abs(y[i]);
      err += (y_err[i] / tol) * (y_err[i] / tol);
    }
    return std::sqrt(err / n_eq);
  }

  void compute_dense_rs(const std::vector<double> &y,
                        const std::vector<double> &y_prev, double h) {
    // Compute the coefficients r1, r2 for Dormand-Prince dense output
    // The dense output polynomial is: y(x_0 + θh) = y_0 + θ*r1 + θ²*r2
    // where θ = (x - x_0) / h ∈ [0, 1]
    
    for (int i = 0; i < n_eq; i++) {
      // r1 corresponds to the first derivative (slope at start of interval)
      r1[i] = h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
      
      // r2 provides the curvature correction using the dense output coefficients
      r2[i] = h * (d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] + d7 * k7[i]);
    }
  }

  int n_eq; // number of equations
  double atol, rtol; // absolute and relative tolerances

  // We impose a minimum and maximum step size (smaller max for better precision)
  const double hmin = 1.0e-12;
  const double hmax = 0.1;

  // These are temporary variables used to store the coefficients. You are allowed to
  // define additional temporary variables if you need them.
  std::vector<double> k1, k2, k3, k4, k5, k6, k7, y_tmp, y_err;

  // vectors that store the results
  std::vector<double> hs;
  std::vector<double> xs;
  std::vector<std::vector<double>> result;

  // These are temporary variables used to store the coefficients
  // used in dense output (only r1 and r2 are needed for Dormand-Prince)
  std::vector<double> r1, r2;

  // c1 is zero, c6 and c7 are 1.0
  const double c2 = 1.0 / 5.0;
  const double c3 = 3.0 / 10.0;
  const double c4 = 4.0 / 5.0;
  const double c5 = 8.0 / 9.0;

  const double a21 = 1.0 / 5.0;
  const double a31 = 3.0 / 40.0;
  const double a32 = 9.0 / 40.0;
  const double a41 = 44.0 / 45.0;
  const double a42 = -56.0 / 15.0;
  const double a43 = 32.0 / 9.0;
  const double a51 = 19372.0 / 6561.0;
  const double a52 = -25360.0 / 2187.0;
  const double a53 = 64448.0 / 6561.0;
  const double a54 = -212.0 / 729.0;
  const double a61 = 9017.0 / 3168.0;
  const double a62 = -355.0 / 33.0;
  const double a63 = 46732.0 / 5247.0;
  const double a64 = 49.0 / 176.0;
  const double a65 = -5103.0 / 18656.0;

  // Note that a71, a72, a73, a74, a75, a76 are essentially the coefficients of the 5th order solution
  const double a71 = 35.0 / 384.0;
  const double a72 = 0.0;
  const double a73 = 500.0 / 1113.0;
  const double a74 = 125.0 / 192.0;
  const double a75 = -2187.0 / 6784.0;
  const double a76 = 11.0 / 84.0;

  // These coefficients are used to estimate the error in the solution. They are essentially
  // the coefficients of the 5th order solution minus the 4th order solution, i.e. b_i - b_i^*.
  const double e1 = 71.0 / 57600.0;
  const double e2 = 0.0;
  const double e3 = -71.0 / 16695.0;
  const double e4 = 71.0 / 1920.0;
  const double e5 = -17253.0 / 339200.0;
  const double e6 = 22.0 / 525.0;
  const double e7 = -1.0 / 40.0;

  // These are the coefficients for dense output
  const double d1 = -12715105075.0 / 11282082432.0;
  const double d2 = 0.0;
  const double d3 = 87487479700.0 / 32700410799.0;
  const double d4 = -10690763975.0 / 1880347072.0;
  const double d5 = 701980252875.0 / 199316789632.0;
  const double d6 = -1453857185.0 / 822651844.0;
  const double d7 = 69997945.0 / 29380423.0;
};


