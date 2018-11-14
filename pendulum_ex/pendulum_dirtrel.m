clear

% Constants/Matrices inherent to the problem
global N;
pendulum_globals();

% Define State Variable (X)
n = 3*N; % size of state vector
nF = 8*N-3; % size of constraint vector

x = zeros(2*N, 1); % dimension of x * N
x_lb = repmat([-pi; -inf], N, 1);
x_ub = repmat([pi; inf], N, 1);

x(1) = 0; x_lb(1) = 0; x_ub(1) = 0;
x(2) = 0; x_lb(2) = 0; x_ub(2) = 0;
x(end-1) = pi; x_lb(end-1) = pi; x_ub(end-1) = pi;
x(end) = 0; x_lb(end) = 0; x_ub(end) = 0;

u = ones(1*N-1, 1) * 0; % dimension of u * N
u_lb = ones(1*N-1, 1) * -1;
u_ub = ones(1*N-1, 1) * 1;

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

F_lb = [-inf; ones(2*(N-1),1)*-1e-15; repmat([-pi; -pi; -pi; -pi; -3; -3], N-1, 1); -pi; -pi; -pi; -pi];
F_ub = [inf; ones(2*(N-1),1)*1e-15; repmat([pi; pi; pi; pi; 3; 3], N-1, 1); pi; pi; pi; pi];
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
snend;

% Plotting Purposes Only
thetas = zeros(1,N);
t = 0:h:h*(N-1);
for i = 1:N
    thetas(i) = X(2*i-1);
end
pendulum_plot(t,thetas)