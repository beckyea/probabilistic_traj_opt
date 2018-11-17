global R Q Q_N E1 D N dof g m L Mp Mc;

R = .1;
Q = [10 0; 0 1];
Q_N = 100 * eye(2);
E1 = [0 0; 0 0];
D = 0.04;
N = 200; % NUMBER OF KNOT POINTS
dof = 2;

g = 9.81;
m = 1; %pendulum mass
L = 1;
Mc = 1; %cart mass
Mp = 1; %pend mass