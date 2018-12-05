global R Q Q_N E1 D N g m L n_u n_x m_test Q_test Q_N_test R_test 

R = .1;
Q = [10 0; 0 1];
Q_N = 10 * eye(2);
E1 = [0 0; 0 0];
D = 0.25;

N = 70; % NUMBER OF KNOT POINTS

n_x = 2;
n_u = 1;

g = 9.81;
m = 1.25;
L = 1;

m_test = .5;
Q_test = [500 0; 0 10];
Q_N_test = [10000 0; 0 10];
R_test = 0.9;