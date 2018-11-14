function [F] = pendulum_userfun(X)
%PENDULUM_USERFUN Defines nonlinear part 
global N Q R;
x = X(1:N*2);
u = X(N*2+1:N*3-1);
h = X(end);
F = zeros(N*8-3, 1);
[l, E, K] = pendulum_lw(x, u, h, Q, R);
F(1) = l;
for i = 1:N
    x_i = [x(2*i-1) x(2*i)];
    if i ~= N
        x_ip1 = [x(2*i+1) x(2*i+2)];
        u_i = u(i);
        F(1) = F(1) + cost(x_i, u_i, h);
        [theta_fun, thetadot_fun] = f_h(x_i, u_i, h, x_ip1);
        F(2*i) = theta_fun;
        F(2*i+1) = thetadot_fun;
        sqrt_KEK = sqrt(K{i}*E{i}*K{i}');
        F(2*N-1+6*i-1) = u_i + sqrt_KEK;
        F(2*N-1+6*i) = u_i - sqrt_KEK; 
    end
    sqrt_E = sqrt(E{i});
    F(2*N-1+6*i-5) = x_i(1) + sqrt_E(1,1);
    F(2*N-1+6*i-4) = x_i(1) - sqrt_E(1,1);
    F(2*N-1+6*i-3) = x_i(1) + sqrt_E(1,2);
    F(2*N-1+6*i-2) = x_i(1) - sqrt_E(1,2);
end
end

function [c] = cost(x, u, h)
global Q R;
dx = [x(1) - pi; x(2)];
c = h;%dx'*Q*dx + u*R*u + h;
end

function [theta_fun, thetadot_fun] = f_h(x, u, h, x_ip1)
global g L m;
    theta_i = x(1);
    thetadot_i = x(2);
    theta_ip1 = x_ip1(1);
    thetadot_ip1 = x_ip1(2);
    theta_fun = theta_ip1 - (theta_i + thetadot_i*h);
    thetadot_fun = thetadot_ip1 - (thetadot_i - g*h*sin(theta_i)/L + u*h/(m*L*L));
%     disp("theta:" + theta_i + " thetad:" + thetadot_i + " theta1:" + theta_ip1...
%         + " thetad1:" + thetadot_ip1 + " theta_fun:" + theta_fun + " thetad_fun" + thetadot_fun);
end
