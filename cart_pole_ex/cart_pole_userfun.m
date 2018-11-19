function [F] = cart_pole_userfun(X)
%PENDULUM_USERFUN Defines nonlinear part 
global N Q R dof n_x;
x = X(1:N*n_x);
u = X(N*n_x+1:end-1);
h = X(end);
F = zeros(N*dof*8-3, 1);
[l, E, K] = cart_pole_lw(x, u, h, Q, R);
F(1) = l;
for i = 1:N
    x_i = X(n_x*(i-1)+1:n_x*i);
    if i ~= N
        x_ip1 = X(n_x*i+1:n_x*(i+1));
        u_i = u(i);
        F(1) = F(1) + cost(x_i, u_i, h);
        [theta_fun, thetadot_fun, y_fun, ydot_fun] = f_h(x_i, u_i, h, x_ip1);
        F(2*dof*i-2) = theta_fun;
        F(2*dof*i-1) = thetadot_fun;
        F(2*dof*i) = y_fun;
        F(2*dof*i+1) = ydot_fun;
%         sqrt_KEK = sqrt(K{i}*E{i}*K{i}');
%         %TODO adjust these indexes, not sure whats going on
%         F(2*N-1+6*i-1) = u_i + sqrt_KEK;
%         F(2*N-1+6*i) = u_i - sqrt_KEK; 
    end
%     sqrt_E = sqrt(E{i});
%     %TODO adjust these indexes, not sure whats going on
%     F(2*N-1+6*i-5) = x_i(1) + sqrt_E(1,1);
%     F(2*N-1+6*i-4) = x_i(1) - sqrt_E(1,1);
%     F(2*N-1+6*i-3) = x_i(1) + sqrt_E(1,2);
%     F(2*N-1+6*i-2) = x_i(1) - sqrt_E(1,2);
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
