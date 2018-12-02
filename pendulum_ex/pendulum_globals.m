global R Q Q_N E1 D N g m L;

R = .1;
Q = [10 0; 0 1];
Q_N = 100 * eye(2);
E1 = [0 0; 0 0];
D = 0.0;

N = 50; % NUMBER OF KNOT POINTS

n_x = 2;
n_u = 1;

g = 9.81;
m = .8;
L = 1;