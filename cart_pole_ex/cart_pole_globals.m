global R Q Q_N E1 D N dof g L Mp Mc mu_k Nc_sign;

R = 1;
Q = [10  0 0 0;
      0 10 0 0;
      0  0 1 0;
      0  0 0 1];
Q_N = 100 * eye(4);
E1 = [0  0 0 0;
      0  0 0 0;
      0  0 0 0;
      0  0 0 0];

D = 0.1;
N = 100; % NUMBER OF KNOT POINTS
dof = 2;
n_u = 1; % dimension of u
n_x = 2*dof; % dimension of state vector
n_c = 1;

g = 9.81;
Mp = 0.5; %pendulum mass
L = 0.5;
Mc = 1; %cart mass
mu_k = 0.0; %coefficient of friction
Nc_sign = 1;