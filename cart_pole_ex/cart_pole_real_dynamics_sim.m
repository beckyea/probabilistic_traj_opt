%simulate cartpole dynamics
%% compute control with dirtrel
%call dirtrel
cart_pole_dirtrel
%% real dynamics simulation
cart_pole_globals()
global Nc_sign
Nc_sign = 1; %set the sign positive to start
tspan = [0, h];
% tspan = [0, 10];
x0 = [0, 0, 0, 0]; %first initial condition
allX = [];
allT = [];
u
for j = 1:N-1
    u_Curr = u(j);
    [t, x] = ode45(@(t,x) compute_cart_pole_dyn(x, t, u_Curr), tspan, x0);
    x0 = x(end, :);
    tspan = [t(end), t(end) + h];
    allX = [allX; x];
    allT = [allT; t];
end
[t, x] = ode45(@(t,x) compute_cart_pole_dyn(x, t, u), tspan, x0);
%% plotting
thetas = allX(:,1);
theatadots = allX(:,2);
ys = allX(:,3);
ydots = allX(:,4);
% thetas = x(:, 1);
% ys = x(:, 3);
% allT = t;
cart_pole_plot(allT, thetas, ys)