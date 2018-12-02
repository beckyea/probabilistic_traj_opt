%% attempt to implement gain scheduled controller for LQR on error dynamics
pendulum_globals()


tspan = [0, h];
x0 = [0, 0]; %first initial condition
allX = [];
allT = [];

%calculate all the gains a functions of state
dirtrelX = [dirtrel_thetas; dirtrel_thetadots]';
% K = cell(1,N);
% for i = 1:N-1
%     A = A_fun([dirtrel_thetas(i); dirtrel_thetadots(i); dirtrel_ys(i); dirtrel_ydots(i)], dirtrel_u(i));
%     B = B_fun([dirtrel_thetas(i); dirtrel_thetadots(i); dirtrel_ys(i); dirtrel_ydots(i)], dirtrel_u(i));
%     K{i} = lqr(A, B, Q, R);
% end

for j = 1:N-1
    ud = dirtrel_u(j);
    [t, x] = ode45(@(t,x) dynamics2(x, t, dirtrelX, ud, j, Q, R), tspan, x0);
    x0 = x(end, :);
    tspan = [t(end), t(end) + h];
    allX = [allX; x];
    allT = [allT; t];
end

% plotting
thetas = allX(:,1);
thetadots = allX(:,2);
% thetas = x(:, 1);
% ys = x(:, 3);
% allT = t;
clf
pendulum_plot(allT, thetas, [-1; zeros(size(thetas,1)-3,1); 1])

function xdot = dynamics2(x, t, dirtrelX, u_d, j, Q, R)
global g L m
theta =  x(1);
thetadot = x(2);

epsilon = (x' - dirtrelX(j,:)) ./ (dirtrelX(j+1, :) - dirtrelX(j, :));
epsilon(isnan(epsilon)) = 1;
epsilon(isinf(epsilon)) = 0;
x_interp = epsilon.*dirtrelX(j,:) + (1-epsilon).*dirtrelX(j+1, :);
%calculate K for this timestep
A = A_fun(x_interp, u_d);
B = B_fun(x_interp, u_d);
K = lqr(A,B,Q,R);
u = K*(x-x_interp') + u_d;
thetaddot =-g/L*sin(theta)+u/(m*L*L);

xdot = [thetadot; thetaddot];

end

function Amat = A_fun(x_d, u_d)
global L g
theta = x_d(1);
Amat = [0 1; -g/L*cos(theta) 0];
 
end

function Bmat = B_fun(x_d, u_d)
global m L
Bmat = [0; 1/(m*L*L)];
end