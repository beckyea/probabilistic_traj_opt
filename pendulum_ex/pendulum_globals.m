global dt R Q Q_N E1 D N g m L;

dt = 0.01;
R = .1;
Q = [10 0; 0 1];
Q_N = 100 * eye(2);
E1 = [0 0; 0 0];
D = 0.04;
N = 50; % NUMBER OF KNOT POINTS


g = 9.81;
m = 1;
L = 1;