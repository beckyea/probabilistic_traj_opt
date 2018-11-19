function [F] = cart_pole_userfun(X)
%PENDULUM_USERFUN Defines nonlinear part 
global N Q R n_u n_x;
x = X(1:N*n_x);
u = X(N*n_x+1:end-1);
h = X(end);
F = zeros(1 + n_x*(N-1) + 2*n_x*n_x*N + 2*n_u*n_u*(N-1), 1);
[l, E, K] = cart_pole_lw(x, u, h, Q, R);
n_c = 1;           % number of cost constraints
n_d = n_x*(N-1); % number of dynamics constraints
n_e_x = n_x*n_x*2; % number of ellipsoidal constraints on x
n_e_u = n_u*n_u*2; % number of ellipsoidal constraints on u
n_e = n_e_x+n_e_u; % number of ellipsoidal constraints

F(1) = l;
for i = 1:N
    x_i = X(n_x*(i-1)+1:n_x*i);
    if i ~= N
        x_ip1 = X(n_x*i+1:n_x*(i+1));
        u_i = u(i);
        F(1) = F(1) + cost(x_i, u_i, h);
        [theta_fun, thetadot_fun, y_fun, ydot_fun] = f_h(x_i, u_i, h, x_ip1);
        % Dynamics Constraints
        F(n_c + n_x*i-3) = theta_fun;
        F(n_c + n_x*i-2) = thetadot_fun;
        F(n_c + n_x*i-1) = y_fun;
        F(n_c + n_x*i)   = ydot_fun;
        sqrt_KEK = sqrt(K{i}*E{i}*K{i}');
        for j = 1:n_u
            u_plus_col = u_i + sqrt_KEK(:, j);
            u_minus_col = u_i - sqrt_KEK(:, j);
        F(n_c + n_d + n_e*(i-1) + 2*n_x*n_x + 1 : n_c + n_d + n_e*(i-1) + 2*n_x*n_x + n_u) = u_plus_col;
        F(n_c + n_d + n_e*(i-1) + 2*n_x*n_x + n_u + 1 : n_c + n_d + n_e*(i-1) + 2*n_x*n_x + 2*n_u) = u_minus_col;
        %n_c + n_d + n_e*(i-1) + 2*n_x*n_x + 1 : n_c + n_d + n_e*(i-1) + 2*n_x*n_x + 2*n_u
        end
        
    end
    sqrt_E = sqrt(E{i});
    for j = 1:n_x
        x_plus_col = x_i + sqrt_E(:, j);
        x_minus_col = x_i - sqrt_E(:, j);
        F(n_c + n_d + n_e*(i-1) + 2*n_x*(j-1)+1 : n_c + n_d + n_e*(i-1) + 2*n_x*(j-1) + n_x) = x_plus_col;
        F(n_c + n_d + n_e*(i-1) + 2*n_x*(j-1)+ n_x +1: n_c + n_d + n_e*(i-1) + 2*n_x*(j-1) + 2*n_x) = x_minus_col;
        %n_c + n_d + n_e*(i-1) + 2*n_x*(j-1)+1 : n_c + n_d + n_e*(i-1) + 2*n_x*(j-1) + 2*n_x
    end
end
end

function [c] = cost(x, u, h)
global Q R y_des;
%TODO double check this
% dx = [x(1) - pi; x(2); x(3) - y_des; x(4)];
dx = [x(1) - pi; 0; 0; 0];
c = dx'*Q*dx + u*R*u + h;
end

function [theta_fun, thetadot_fun, y_fun, ydot_fun] = f_h(x, u, h, x_ip1)
global g L Mp Mc;
    theta_i = x(1);
    thetadot_i = x(2);
    thetaddot_i = thetadot_i - (h*(L*Mp*cos(theta_i)*sin(theta_i)*thetadot_i^2 +...
        u*cos(theta_i) + Mc*g*sin(theta_i) + Mp*g*sin(theta_i)))/...
        (L*(Mc + Mp - Mp*cos(theta_i)^2));
    y_i = x(3);
    ydot_i = x(4);
    yddot_i = ydot_i + (h*(L*Mp*sin(theta_i)*thetadot_i^2 + u + ...
        Mp*g*cos(theta_i)*sin(theta_i)))/(Mc + Mp - Mp*cos(theta_i)^2);
    
    theta_ip1 = x_ip1(1);
    thetadot_ip1 = x_ip1(2);
    y_ip1 = x_ip1(3);
    ydot_ip1 = x_ip1(4);
    
    theta_fun = theta_ip1 - (theta_i + thetadot_i*h);
    thetadot_fun = thetadot_ip1 - (thetaddot_i);
    y_fun = y_ip1 - (y_i + ydot_i*h);
    ydot_fun = ydot_ip1 - (yddot_i);
%     disp("theta:" + theta_i + " thetad:" + thetadot_i + " theta1:" + theta_ip1...
%         + " thetad1:" + thetadot_ip1 + " theta_fun:" + theta_fun + " thetad_fun" + thetadot_fun);
end
