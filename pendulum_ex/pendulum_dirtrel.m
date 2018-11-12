clear

% Constants/Matrices inherent to the problem
global N;
pendulum_globals();

% Define State Variable (x)
x = zeros(1, 2*N); % dimension of x * N
x_lb = repmat([-pi, -inf], 1, N);
x_ub = repmat([pi, inf], 1, N);
u = zeros(1, 1*N-1); % dimension of u * N
u_lb = ones(1, 1*N-1) * -3;
u_ub = ones(1, 1*N-1) * 3;
h = 0; % time decision variable
h_lb = 0;
h_ub = 10;

X = [x u h]; % state variable for snopt
X_lb = [x_lb u_lb h_lb]; % xlow in SNOPT documentation
X_ub = [x_ub u_ub h_ub]; % xupp in SNOPT documentation
xmul = []; xstate = [];
Fmul = []; Fstate = [];

F_lb = [-inf, zeros(1, 2*N)];
F_ub = [inf, zeros(1, 2*N)];

% Set SNOPT Parameters
snprint('pendulum_snopt.out');

pendulum_snopt.spc = which('pendulum_snopt.spc');
snspec(pendulum_snopt.spc);

[X,F,INFO] = snopt (X, x_lb, x_ub, xmul, xstate,...
                   F_lb, F_ub, Fmul, Fstate, @pendulum_userfun);

% Plotting Purposes Only
thetas = zeros(1,N);
t = 0:h:h*(N-1);
for i = 1:N
    thetas(i) = X(2*i-1);
end
pendulum_plot(t,thetas)