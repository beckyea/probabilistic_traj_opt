function [l, E, K] = cart_pole_lw(x, u, h, Q_l, R_l)
%PENDULUM_LW Robust Cost Function Described in Dirtrel Paper (Algorithm 1)

global R Q Q_N E1 D N g m L Mp Mc;

if N-1 ~= size(u)
    print("Size of X != Size of U+1")
    return
end

l = 0;
A = cell(N,1);
B = cell(N,1);
G = cell(N-1,1);
K = cell(N-1,1);
H = cell(N,1);
E = cell(N,1);
E{1} = E1;
for i = 1:N
    theta_i = x(4*i-3);
    thetadot_i = x(4*i-2);
    y_i = x(4*i-1);
    ydot_i = x(4*i); 
    u_i = u(i);
    
    dthdd_dth = (2*Mp*h*cos(theta_i)*sin(theta_i)*(L*Mp*cos(theta_i)*sin(theta_i)*thetadot_i^2 ...
        + u_i*cos(theta_i) + Mc*g*sin(theta_i) + Mp*g*sin(theta_i)))/(L*(- Mp*cos(theta_i)^2 ...
        + Mc + Mp)^2) - (2*h*(L*Mp*(2*cos(theta_i)^2 - 1)*thetadot_i^2 - u_i*sin(theta_i) ...
    + Mc*g*cos(theta_i) + Mp*g*cos(theta_i)))/(L*(2*Mc + Mp - Mp*(2*cos(theta_i)^2 - 1)));
    
    dthdd_dthd = - (Mp*h*thetadot_i*sin(2*theta_i))/(Mc + Mp/2 - (Mp*cos(2*theta_i))/2);
   
    dydd_dth = (Mp*h*(2*g*cos(theta_i)^2 - g + L*thetadot_i^2*cos(theta_i)))...
        /(Mc + Mp - Mp*cos(theta_i)^2) - (2*Mp*h*cos(theta_i)*sin(theta_i)*...
        (L*Mp*sin(theta_i)*thetadot_i^2 + u_i + Mp*g*cos(theta_i)*sin(theta_i)))...
        /(- Mp*cos(theta_i)^2 + Mc + Mp)^2;
    
    dydd_dthd = (2*L*Mp*h*thetadot_i*sin(theta_i))/(Mp*sin(theta_i)^2 + Mc);
    
    A{i} = [        1,              h, 0, 0;
            dthdd_dth, 1 + dthdd_dthd, 0, 0;
                    0,              0, 1, h;
             dydd_dth,      dydd_dthd, 0, 1];
    
    B{i} = [                                          0;
            -cos(theta)/(L*(Mc + Mp - Mp*cos(theta)^2));
                                                      0;
                          1/(Mc + Mp - Mp*cos(theta)^2)];
    
    if i ~= N
        %TODO
        G{i} = [ 0 0; -h*u(i)/(L*L*m*m) 0]; %TODO - FIX THIS
    end
end
P = cell(N,1);
P{N} = Q_N;
for i = N:-1:2
    P{i-1} = Q + A{i}'*P{i}*A{i}-A{i}'*P{i}*B{i}*...
        inv(R+B{i}'*P{i}*B{i})*(B{i}'*P{i}*A{i});
end
for i = 1:N-1
    K{i} = inv(R+B{i}'*P{i+1}*B{i})*(B{i}'*P{i+1}*A{i});
end
H{1} = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
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

