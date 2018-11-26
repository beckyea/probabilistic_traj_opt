clear

% Constants/Matrices inherent to the problem
global N n_x n_u;
cart_pole_globals();

% Define State Variable (X)
n = (n_u+n_x)*N; % size of state vector
nF = zeros(1 + n_x*(N-1) + 2*n_x*n_x*N + 2*n_u*n_u*(N-1), 1); % size of constraint vector
   % cost, dynamics, x ellipsoids, y ellipsoids

x = zeros(n_x*N, 1); % dimension of x * N (ie 2*dof*N)
x_const_lower = [-pi; -inf; -2; -inf]; % state constraint lower bound
x_const_upper = [ pi;  inf;  2;  inf]; % state constraint upper bound
x_lb = repmat(x_const_lower, N, 1);
x_ub = repmat(x_const_upper, N, 1);

x(1) = 0;      x_lb(1) = 0;      x_ub(1) = 0;
x(2) = 0;      x_lb(2) = 0;      x_ub(2) = 0;
x(3) = 0;      x_lb(3) = 0;      x_ub(3) = 0;
x(4) = 0;      x_lb(4) = 0;      x_ub(4) = 0;
%x(end-1) = 1;  x_lb(end-1) = 1;  x_ub(end-1) = 1;
x(end) = 0;    x_lb(end) = 0;    x_ub(end) = 0;
x(end-3) = pi; x_lb(end-3) = pi; x_ub(end-3) = pi;
x(end-2) = 0;  x_lb(end-2) = 0;  x_ub(end-2) = 0;

u = zeros(n_u*(N-1), 1); % dimension of u * N
u_const_lower = -1; % u constraint lower bound
u_const_upper = 1;  % u constraint upper bound
u_lb = repmat(u_const_lower, n_u*(N-1), 1);
u_ub = repmat(u_const_upper, n_u*(N-1), 1);

h = 0.001; % time decision variable
h_lb = 0.0001;
h_ub = 0.05;

ObjAdd = 0;
ObjRow = 1;

X = [x; u; h]; % state variable for snopt
X_lb = [x_lb; u_lb; h_lb]; % xlow in SNOPT documentation
X_ub = [x_ub; u_ub; h_ub]; % xupp in SNOPT documentation
xmul = []; 
xstate = [];

% Note: +/- 1e-15 used in leiu of 0 for floating point errors
F_lb = [-inf; ones(n_x*(N-1),1)*-1e-15; ...
    repmat([repmat(x_const_lower, n_x*2, 1); repmat(u_const_lower, n_u*2, 1)], N-1, 1);...
    repmat(x_const_lower, n_x*2, 1)];
F_ub = [inf; ones(n_x*(N-1),1)*1e-15; ...
    repmat([repmat(x_const_upper, n_x*2, 1); repmat(u_const_upper, n_u*2, 1)], N-1, 1);...
    repmat(x_const_upper, n_x*2, 1)];
Fmul = []; 
Fstate = [];


% Set SNOPT Parameters
snscreen on;
snprint('cart_pole_snopt.out');

cart_pole_snopt.spc = which('cart_pole_snopt.spc');
snspec(cart_pole_snopt.spc);
options.name = 'cart_pole';
snset ('Minimize');

[X,F,INFO,xmul,Fmul,xstate,Fstate,output] = snopt (X, X_lb, X_ub, xmul, xstate,...
                   F_lb, F_ub, Fmul, Fstate, @cart_pole_userfun,...
                   ObjAdd, ObjRow, options);

snprint off; % Closes the file and empties the print buffer
snend;

u = X(N*n_x:end-1);

% Plotting Purposes Only
thetas = zeros(1,N);
y_pos = zeros(1,N);
h = X(end);
t = 0:h:h*(N-1);
% thetas = X(1:4:4*dof*N);
% y_pos = X(3:4:4*dof*N); %double check this indexing

for i = 1:N
    thetas(i) = X(4*i-3);
    y_pos(i) = X(4*i-1);
end
cart_pole_plot(t,thetas,y_pos)