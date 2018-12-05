clear

% Constants/Matrices inherent to the problem
global N  n_x n_u;
pendulum_globals();

% Define State Variable (X)
n = (n_u+n_x)*N; % size of state vector
nF = zeros(1 + n_x*(N-1) + 2*n_x*n_x*N + 2*n_u*n_u*(N-1), 1); % size of constraint vector
   % cost, dynamics, x ellipsoids, y ellipsoids

x = zeros(n_x*N, 1); % dimension of x * N (ie 2*dof*N)
x_const_lower = [-pi; -inf]; % state constraint lower bound
x_const_upper = [ pi;  inf]; % state constraint upper bound
x_lb = repmat(x_const_lower, N, 1);
x_ub = repmat(x_const_upper, N, 1);

x(1) = 0;      x_lb(1) = 0;      x_ub(1) = 0;
x(2) = 0;      x_lb(2) = 0;      x_ub(2) = 0;
x(end-1) = pi; x_lb(end-1) = pi; x_ub(end-1) = pi;
x(end) = 0;    x_lb(end) = 0;    x_ub(end) = 0;

u = zeros(n_u*(N-1), 1); % dimension of u * N
u_const_lower = -10; % u constraint lower bound
u_const_upper = 10;  % u constraint upper bound
u_lb = repmat(u_const_lower, n_u*(N-1), 1);
u_ub = repmat(u_const_upper, n_u*(N-1), 1);

h = 0.01; % time decision variable
h_lb = 0.001;
h_ub = 0.1;

ObjAdd = 0;
ObjRow = 1;

X = [x; u; h]; % state variable for snopt
X_lb = [x_lb; u_lb; h_lb]; % xlow in SNOPT documentation
X_ub = [x_ub; u_ub; h_ub]; % xupp in SNOPT documentation
xmul = []; 
xstate = [];

% Note: +/- 1e-15 used in leiu of 0 for floating point errors
F_lb = [-inf; ones(n_x*(N-1),1)*1e-20; ...
    repmat([repmat(x_const_lower, n_x*2, 1); repmat(u_const_lower, n_u*2, 1)], N-1, 1);...
    repmat(x_const_lower, n_x*2, 1)];
F_ub = [inf; ones(n_x*(N-1),1)*1e-20; ...
    repmat([repmat(x_const_upper, n_x*2, 1); repmat(u_const_upper, n_u*2, 1)], N-1, 1);...
    repmat(x_const_upper, n_x*2, 1)];
Fmul = [];
Fstate = [];


% Set SNOPT Parameters
snscreen on;
snprint('pendulum_snopt.out');

pendulum_snopt.spc = which('pendulum_snopt.spc');
snspec(pendulum_snopt.spc);

options.name = 'pendulum';
snset ('Minimize');

[X,F,INFO,xmul,Fmul,xstate,Fstate,output] = snopt (X, X_lb, X_ub, xmul, xstate,...
                   F_lb, F_ub, Fmul, Fstate, @pendulum_userfun,...
                   ObjAdd, ObjRow, options);

snprint off; % Closes the file and empties the print buffer
snend

x = X(1:N*n_x);
u = X(N*n_x+1:end-1);

% Plotting Purposes Only
dirtrel_thetas = zeros(N,1);
t = 0:h:h*(N-1);
for i = 1:N
    dirtrel_thetas(i) = X(2*i-1);
end
%%
dirtrel_thetadots = zeros(N,1);
for i = 1:N
    dirtrel_thetadots(i) = X(2*i);
end
dirtrel_u = zeros(N-1,1);
for i = 1:N-1
    dirtrel_u(i) = X(2*N+i);
end
pendulum_plot(t,dirtrel_thetas, u, 1)