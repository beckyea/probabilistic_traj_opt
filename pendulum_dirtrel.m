clear
R = .1;
Q = [10 0; 0 1];
Q_N = 100 * eye(2);
E1 = [0 0; 0 0];
D = 0.04;


% TODO: Solve for x, u from optimal nominal trajectory
x = {[0 0], [0 .1], [.1, .2]};
u = {1, 1};

l_w(x, u, D, E1, Q, R, Q_N, Q, R)


function [l] = l_w(x, u, D, E1, Q_l, R_l, Q_N, Q, R)
% x - a cell array of states {[theta1 thetadot1] [theta2 thetadot2]....}
% u - a cell array of control input {in1 in2 ...}
% D - disturbance (scalar)
% E1 - initial noise (matrix)
% Q_l - matrix
% R_l - matrix
% Q_N - final matrix
% Q - matrix
% R - matrix

% Other parameters defined
h = 0.01;  % is this a decision var?? stay tuned (TODO!)
g = 9.81;
L = 1;
m = 1.0;

N = size(x,2);
if N-1 ~= size(u)
    print("Size of X != Size of U+1")
    return
end
l = 0;
A = cell(N);
B = cell(N);
G = cell(N);
K = cell(N);
H = cell(N);
E = cell(N);
E{1} = E1;
for i = 1:N-1
    theta_i = x{i}(1);
    theta_dot_i = x{i}(2);
    u_i = u{i};
    A{i} = [ 1 h; (-g/L*cos(theta_i)*theta_dot_i*h) 1 ];
    B{i} = [ 0; h/(m*L^2) ];
    G{i} = [ 0 0; -h*u_i/(L^2*m^2) 0];
    K{i} = lqr(A{i}, B{i}, Q, R, 0);
end
H{1} = [0 0; 0 0];
for i = 1:N-1
    l = l + trace((Q_l+K{i}'*R_l*K{i})*E{i});
    E{i+1} = (A{i}-B{i}*K{i})*E{i}*(A{i}-B{i}*K{i})' + ...
        (A{i}-B{i}*K{i})*H{i}*G{i}' + ...
        G{i}*H{1}'*(A{i}-B{i}*K{i})' + ...
        G{i}*D*G{i}';
    H{i+1} = (A{i}-B{i}*K{i})*H{i}+G{i}*D;
end
l = l + trace(Q_N*E{N});

end
