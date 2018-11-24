%simulate cartpole dynamics
%% simulation
cart_pole_globals()
global Nc_sign
Nc_sign = 1; %set the sign positive to start
tspan = [0, 10];
x0 = [pi/3, 0, 0, 0];
u = 0;
% for j = 1:length(traj)
[t, x] = ode45(@(t,x) compute_cart_pole_dyn(x, t, u), tspan, x0);

%% plotting
thetas = x(:,1);
theatadots = x(:,2);
ys = x(:,3);
ydots = x(:,4);

cart_pole_plot(t, thetas, ys)