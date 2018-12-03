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
        end
        
    end
    sqrt_E = sqrt(E{i});
    for j = 1:n_x
        x_plus_col = x_i + sqrt_E(:, j);
        x_minus_col = x_i - sqrt_E(:, j);
        F(n_c + n_d + n_e*(i-1) + 2*n_x*(j-1)+1 : n_c + n_d + n_e*(i-1) + 2*n_x*(j-1) + n_x) = x_plus_col;
        F(n_c + n_d + n_e*(i-1) + 2*n_x*(j-1)+ n_x +1: n_c + n_d + n_e*(i-1) + 2*n_x*(j-1) + 2*n_x) = x_minus_col;
    end
end
end

function [c] = cost(x, u, h)
global Q R
%TODO double check this
% dx = [x(1) - pi; x(2); x(3) - y_des; x(4)];
dx = [x(1) - pi; 0; 0; 0];
c = dx'*dx + 0.1*u*u;
end

function [theta_fun, thetadot_fun, y_fun, ydot_fun] = f_h(x_i, u_i, h, x_ip1)
    global g L Mp Mc;
    theta_i    = x_i(1);
    thetadot_i = x_i(2);
    y_i        = x_i(3);
    ydot_i     = x_i(4);
    
    theta_ip1    = x_ip1(1);
    thetadot_ip1 = x_ip1(2);
    y_ip1        = x_ip1(3);
    ydot_ip1     = x_ip1(4);
    
    
%     thetaddot_i = -(L*Mp*cos(theta_i)*sin(theta_i)*thetadot_i^2 +...
%         u_i*cos(theta_i) + Mc*g*sin(theta_i) + Mp*g*sin(theta_i))/...
%         (L*(Mc + Mp - Mp*cos(theta_i)^2));
%     yddot_i = (L*Mp*sin(theta_i)*thetadot_i^2 + u_i + ...
%         Mp*g*cos(theta_i)*sin(theta_i))/(Mc + Mp - Mp*cos(theta_i)^2);
    

    thetaddot_i = (-(Mp+Mc)*g*sin(theta_i)-cos(theta_i)*(u_i+Mp*L*thetadot_i^2*sin(theta_i)))/...
        (L*(Mp+Mc)-Mp*L*cos(theta_i)*cos(theta_i));
    
    yddot_i = (u_i + Mp*L*(thetadot_i^2*sin(theta_i)-thetaddot_i*cos(theta_i)))/(Mc+Mp);
    
    theta_fun    = theta_ip1    - (theta_i    + thetadot_i*h);
    thetadot_fun = thetadot_ip1 - (thetadot_i + thetaddot_i*h);
    y_fun        = y_ip1        - (y_i        + ydot_i*h);
    ydot_fun     = ydot_ip1     - (ydot_i     + yddot_i*h);
    
%     disp("theta:" + theta_i + " thetad:" + thetadot_i + " theta1:" + theta_ip1...
%         + " thetad1:" + thetadot_ip1 + " theta_fun:" + theta_fun + " thetad_fun" + thetadot_fun);
end
