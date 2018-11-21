global R Q Q_N E1 D N dof g m L Mp Mc y_des mu_k Nc_sign;

R = 1;
Q = [10  0 0 0;
      0 10 0 0;
      0  0 1 0;
      0  0 0 1];
Q_N = 100 * eye(4);
E1 = zeros(4,4);
D = 4;
N = 50; % NUMBER OF KNOT POINTS
dof = 2;
n_u = 1; % dimension of u
n_x = 2*dof; % dimension of state vector

y_des = 1;

g = 9.81;
m = 1; %pendulum mass
L = 1;
Mc = 1.2; %cart mass
Mp = .4; %pend mass
mu_k = 0.02; %coefficient of friction
Nc_sign = 1;