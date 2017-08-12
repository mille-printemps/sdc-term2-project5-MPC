# Model

The model used is the same as the one in the lectures.

```
State = [x, y, psi, v]
x = X coordinate of the vehicle position
y = Y coordinate of the vehicle position
psi = Orientation of the vehicle
v = Velocity of the vehicle

Actuators = [delta, a]
delta = Steering angle [-25, 25]
a = Throttle [-1, 1]

Update equations
x_t+1 = x_t + v_t * cos(psi_t) * dt
y_t+1 = y_t + v_t * sin(psi_t) * dt
​​psi_t+1 = psi_t + (v_t/L_f) * delta_t * dt
v_t+1 = v_t + a_t * dt

​where dt is the time that elapses and L_f is the distance between the front of the vehicle and its center of gravity.
```

# Timestep Length and Elapsed Duration

```
const size_t N = 20;
const double dt = 0.1;
```
These values were chosen by trial-and-error. The reference velocity `ref_v` was set to `40` MPH. For this reference velocity, smaller values than `0.1` did not work. However, an appropriate `dt` value could not be found a reference velocity faster than `40`. So eventually, `[N, dt, ref_v] = [20, 0.1, 40]` ended up being the best combination.

# Polynomial Fitting

A 3rd degree polynomial was used as indicated in `main.cpp`.

```
          Eigen::VectorXd coeffs = polyfit(xvals, yvals, 3);

```

# Model Predictive Control with Latency

Much time was taken to figure out this part. Basically, the latency should be able to be incorporated into the model by replacing `dt` with the latency in the update equation. 

However, the problem is where the vehicle is. In the implementation of `main.cpp`, the positions of the vehicle gotten from the server were converted into the vehicle's coordinate once, and then the state values were calculated. 

Namely, since the current state of the vehicle `[x_t, y_t, psi_t, v_t]` became `[0,  0, 0, v_t]` as the result of the conversion, the state values were calculated using `[0, 0, 0, v_t]` as the current state. The coefficients of the polynomial were also calculated based on this state. 

```
          Eigen::VectorXd state(6);
          state <<
                  v * cos(psi) * LATENCY,               // x
                  v * sin(psi) * LATENCY,               // y
                  v * (-steering / MPC::LF) * LATENCY,  // psi
                  v + throttle * LATENCY,               // v
                  cte,
                  epsi;
                  
          vector<double> next_state = mpc.Solve(state, coeffs);
```

