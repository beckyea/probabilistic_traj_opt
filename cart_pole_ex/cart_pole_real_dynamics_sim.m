%simulate cartpole dynamics
%% compute control with dirtrel
%call dirtrel
figure(1)
cart_pole_dirtrel
%%
dirtrel_x = x;
dirtrel_u = u;
dirtrel_thetas = x(1:4:end);
dirtrel_thetadots = x(2:4:end);
dirtrel_ys = x(3:4:end);
dirtrel_ydots = x(4:4:end);

%% real dynamics simulation
cart_pole_globals()
global Nc_sign
Nc_sign = 1; %set the sign positive to start
tspan = [0, h];
%tspan = [0, 10];
x0 = [0, 0, 0, 0]; %first initial condition
allX = [];
allT = [];
for j = 1:N-1
    u_Curr = u(j);
    [t, x] = ode45(@(t,x) compute_cart_pole_dyn(x, t, u_Curr), tspan, x0);
    x0 = x(end, :);
    tspan = [t(end), t(end) + h];
    allX = [allX; x];
    allT = [allT; t];
end
% [t, x] = ode45(@(t,x) compute_cart_pole_dyn(x, t, 0), tspan, x0);
%% plotting
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

figure(3)
subplot(5,1,1)
hold on
plot(0:h:h*N-h, dirtrel_thetas)
plot(allT, thetas)
legend('$\theta_{dirtrel}$','$\theta_{dyn}$', 'Interpreter','latex');
hold off
subplot(5,1,2)
hold on
plot(0:h:h*N-h, dirtrel_thetadots)
plot(allT, thetadots)
legend('$\dot{\theta}_{dirtrel}$','$\dot{\theta}_{dyn}$', 'Interpreter','latex');
hold off
subplot(5,1,3)
hold on
plot(0:h:h*N-h, dirtrel_ys)
plot(allT, ys)
legend('$y_{dirtrel}$','$y_{dyn}$', 'Interpreter','latex');
hold off
subplot(5,1,4)
hold on
plot(0:h:h*N-h, dirtrel_ydots)
plot(allT, ys)
legend('$\dot{y}_{dirtrel}$','$\dot{y}_{dyn}$', 'Interpreter','latex');
hold off
subplot(5,1,5)
hold on
plot(1:N-1, dirtrel_u,'-')
plot(1:N-1, u,'--')
xlim([0 N]);
legend('$u_{dirtrel}$','$u_{dyn}$', 'Interpreter','latex');
hold off