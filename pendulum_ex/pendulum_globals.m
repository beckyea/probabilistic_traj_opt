global R Q Q_N E1 D N g m L n_c n_x n_u;

R = .1;
Q = [10 0; 0 1];
Q_N = 100 * eye(2);
E1 = [0 0; 0 0];
D = 0.0;

N = 150; % NUMBER OF KNOT POINTS

n_x = 2;
n_u = 1;
n_c = 1; 
g = 9.81;
m = 0.4;
L = 1;