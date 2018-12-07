function [Amat, iAfun, jAvar, iGfun, jGvar] = pendulum_const_derivs()

global n_x n_u n_c N
pendulum_globals();

A = [];
G = [];

for i = 1:N-1
    % Returns A in the format [Ai, Aj, A]
    thi = n_x*i-1;
    thdi = n_x*i;
    thip1 = n_x*i+1;
    thdip1 = n_x*i+2;
    ui = n_x*N+n_u*(i-1)+1;
    hi = n_x*N+n_u*(N-1)+1;

    th_f = n_c + n_x*i-1;
    thd_f = n_c + n_x*i;


    A((i-1)*4+1:i*4,:) = [thi,     th_f,  -1;
                          thdi,    thd_f, -1;
                          thip1,   th_f,   1;
                          thdip1,  thd_f,  1];
                      
   
    G((i-1)*5+1:i*5,:)= [thdi, th_f;
                         thi, thd_f;
                         ui, thd_f;
                         hi, th_f;
                         hi, thd_f];

end

for i = 1:N*n_x
    G(end+1,:) = [i, 1];
end


Amat = A(:,3);
iAfun = A(:,2);
jAvar = A(:,1);
iGfun = G(:,2);
jGvar = G(:,1);

end