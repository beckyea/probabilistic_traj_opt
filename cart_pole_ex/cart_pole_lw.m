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
%Note: we are doing a weird thing where we are setting u(N) = 0
% because we need something for A{N}
for i = 1:N
    theta_i = x(4*i-3);
    thetadot_i = x(4*i-2);
    y_i = x(4*i-1);
    ydot_i = x(4*i); 
    if i == N
        u_i = 0;
    else
        u_i = u(i);
    end
    
%     dthdd_dth = (2*Mp*h*cos(theta_i)*sin(theta_i)*(L*Mp*cos(theta_i)*sin(theta_i)*thetadot_i^2 ...
%         + u_i*cos(theta_i) + Mc*g*sin(theta_i) + Mp*g*sin(theta_i)))/(L*(- Mp*cos(theta_i)^2 ...
%         + Mc + Mp)^2) - (2*h*(L*Mp*(2*cos(theta_i)^2 - 1)*thetadot_i^2 - u_i*sin(theta_i) ...
%     + Mc*g*cos(theta_i) + Mp*g*cos(theta_i)))/(L*(2*Mc + Mp - Mp*(2*cos(theta_i)^2 - 1)));
%     
%     dthdd_dthd = - (Mp*h*thetadot_i*sin(2*theta_i))/(Mc + Mp/2 - (Mp*cos(2*theta_i))/2);
%    
%     dydd_dth = (Mp*h*(2*g*cos(theta_i)^2 - g + L*thetadot_i^2*cos(theta_i)))...
%         /(Mc + Mp - Mp*cos(theta_i)^2) - (2*Mp*h*cos(theta_i)*sin(theta_i)*...
%         (L*Mp*sin(theta_i)*thetadot_i^2 + u_i + Mp*g*cos(theta_i)*sin(theta_i)))...
%         /(- Mp*cos(theta_i)^2 + Mc + Mp)^2;
%     
%     dydd_dthd = (2*L*Mp*h*thetadot_i*sin(theta_i))/(Mp*sin(theta_i)^2 + Mc);
%     
%     A{i} = [        1,              h, 0, 0;
%             dthdd_dth, 1 + dthdd_dthd, 0, 0;
%                     0,              0, 1, h;
%              dydd_dth,      dydd_dthd, 0, 1];

    c_th = cos(theta_i);
    s_th = sin(theta_i);
    Mt = Mp+Mc;

    A{i} = [1, h, 0, 0;
           (-(h*(g*c_th+(s_th*(u_i+L*Mp*s_th))/Mt...
           -(L*Mp*c_th^2)/Mt))/(L*((Mp*c_th^2)/Mt-4/3))...
           -(2*Mp*h*c_th*s_th*(g*s_th-(c_th*(u_i+L*Mp*s_th))/Mt))/...
           (L*((Mp*c_th^2)/Mt-4/3)^2*Mt)), 1, 0, 0;
           0, 0, 1, h;
           ((L*Mp*h*(thetadot_i^2*c_th-(s_th*(g*s_th-(c_th*(u_i+L*Mp*s_th))/Mt))...
           /(L*((Mp*c_th^2)/Mt-4/3))+(c_th*(g*c_th+(s_th*(u_i+L*Mp*s_th))/Mt...
           -(L*Mp*c_th^2)/Mt))/(L*((Mp*c_th^2)/Mt-4/3))+(2*Mp*c_th^2*s_th*...
           (g*s_th-(c_th*(u_i+L*Mp*s_th))/Mt))/(L*((Mp*c_th^2)/Mt-4/3)^2*Mt)))/Mt), ...
           ((2*L*Mp*h*thetadot_i*s_th)/Mt), 0, 1];
       
       B{i} = [0;
              h*c_th/(L*((Mp*c_th^2)/Mt-4/3)*Mt);
              0;
              h*(8*Mc + 8*Mp)/(2*(Mc + Mp)*(4*Mc + 4*Mp - 3*Mp*c_th^2))];

% B{i} = [                                          0;
%             -cos(theta_i)/(L*(Mc + Mp - Mp*cos(theta_i)^2));
%                                                       0;
%                           1/(Mc + Mp - Mp*cos(theta_i)^2)];
    
    if i ~= N
        G{i} = [0, 0, 0, 0;
                -(3*Mp*h*c_th*(4*L*Mc-3*u_i*c_th+...
                4*L*Mp-3*L*Mp*c_th^2+3*Mc*g*s_th+...
                3*Mp*g*s_th-3*L*Mp*c_th*s_th))/...
                (L*(-3*Mp*c_th^2+4*Mc+4*Mp)^2), 0, 0, 0;
                0, 0, 0, 0;
                (h*(L*Mp*(c_th*thetadot_i^2+...
                (s_th*(g*s_th-(c_th*(u_i+L*Mp*s_th))/...
                Mt)*(3*Mc+3*Mp))/(L*(-3*Mp*c_th^2+4*Mc+4*Mp)))...
                -g*Mt+(3*Mp^2*c_th^2*(4*L*Mc-3*u_i*c_th+...
                4*L*Mp-3*L*Mp*c_th^2+3*Mc*g*s_th+3*Mp*g*s_th...
                -3*L*Mp*c_th*s_th))/(-3*Mp*c_th^2+...
                4*Mc+4*Mp)^2))/Mt, 0, 0, 0];
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

