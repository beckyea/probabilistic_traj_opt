%% attempt to implement gain scheduled controller for LQR on error dynamics
cart_pole_globals()
global Nc_sign
Nc_sign = 1; %set the sign positive to start
tspan = [0, h];
x0 = [0, 0, 0, 0]; %first initial condition
allX = [];
allT = [];

%calculate all the gains a functions of state
K = cell(1,N);
for i = 1:N-1
    A = A_fun([dirtrel_thetas(i+1); dirtrel_thetadots(i+1); dirtrel_ys(i+1); dirtrel_ydots(i+1)], dirtrel_u(i));
    B = B_fun([dirtrel_thetas(i+1); dirtrel_thetadots(i+1); dirtrel_ys(i+1); dirtrel_ydots(i+1)], dirtrel_u(i));
    K{i} = lqr(A, B, Q, R);
end

for j = 1:N-1
    u_Curr = u(j);
    K_curr = K{i};
    xd = [dirtrel_thetas(i+1); dirtrel_thetadots(i+1); dirtrel_ys(i+1); dirtrel_ydots(i+1)];
    ud = dirtrel_u(i);
    [t, x] = ode45(@(t,x) dynamics2(x, t, K_curr, xd, ud), tspan, x0);
    x0 = x(end, :);
    tspan = [t(end), t(end) + h];
    allX = [allX; x];
    allT = [allT; t];
end

% plotting
thetas = allX(:,1);
thetadots = allX(:,2);
ys = allX(:,3);
ydots = allX(:,4);
% thetas = x(:, 1);
% ys = x(:, 3);
% allT = t;
figure(2)
clf
cart_pole_plot(allT, thetas, ys)

function xdot = dynamics2(x, t, K_curr, x_d, u_d)
global Mp Mc g L mu_k  
theta =  x(1);
thetadot = x(2);
y = x(3);
ydot = x(4);
u = K_curr*(x-x_d) + u_d;
thetaddot = (-(Mp+Mc)*g*sin(theta)-cos(theta)*(u+Mp*L*thetadot^2*sin(theta)))/...
        (L*(Mp+Mc)-Mp*L*cos(theta)*cos(theta));
% 
% yddot = (u + Mp*L*(thetadot^2*sin(theta)-thetaddot*cos(theta)))/(Mc+Mp);
Nc = 0; S = 0;
yddot = (u-mu_k*Nc*S+Mp*L*(thetadot^2*sin(theta)-thetaddot*cos(theta)))/(Mc+Mp);

xdot = [thetadot; thetaddot; ydot; yddot];

end

function Amat = A_fun(x_d, u_d)
cart_pole_globals();
global Mc Mp L
theta = x_d(1);
thetadot = x_d(2);
u = u_d;
Amat = [                                                                                                                                                                                                                                                                 0,                                                           1, 0, 0;
         (2*Mp*cos(theta)*sin(theta)*(cos(theta)*(L*Mp*sin(theta)*thetadot^2 + u) + g*sin(theta)*(Mc + Mp)))/(L*(Mc + Mp - Mp*cos(theta)^2)^2) - (L*Mp*(2*cos(theta)^2 - 1)*thetadot^2 - u*sin(theta) + Mc*g*cos(theta) + Mp*g*cos(theta))/(L*(Mc + Mp - Mp*cos(theta)^2)), -(2*Mp*thetadot*sin(2*theta))/(2*Mc + Mp - Mp*cos(2*theta)), 0, 0;
                                                                                                                                                                                                                                                                         0,                                                           0, 0, 1;
                                                                 (Mp*(4*Mc*g*cos(2*theta) - 4*u*sin(2*theta) - 2*Mp*g + 2*Mp*g*cos(2*theta) + L*Mp*thetadot^2*cos(3*theta) + 4*L*Mc*thetadot^2*cos(theta) - L*Mp*thetadot^2*cos(theta)))/(4*(Mc + Mp - Mp*cos(theta)^2)^2),         (2*L*Mp*thetadot*sin(theta))/(Mc + Mp*sin(theta)^2), 0, 0];
 
end

function Bmat = B_fun(x_d, u_d)
cart_pole_globals();
global Mc Mp L
theta = x_d(1);
 Bmat = [                                            0;
         -cos(theta)/(L*(Mc + Mp) - L*Mp*cos(theta)^2);
                                                     0;
                         1/(Mc + Mp - Mp*cos(theta)^2)];
end