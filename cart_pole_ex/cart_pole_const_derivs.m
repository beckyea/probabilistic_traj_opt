function [Amat, iAfun, jAvar, iGfun, jGvar] = cart_pole_const_derivs()

global n_x n_u n_c N
cart_pole_globals();

A = [];
G = [];

for i = 1:N-1
    % Returns A in the format [Ai, Aj, A]
    thi = n_x*(i-1)+1;
    thdi = n_x*(i-1)+2;
    yi = n_x*(i-1)+3;
    ydi = n_x*(i-1)+4;
    thip1 = n_x*(i-1)+5;
    thdip1 = n_x*(i-1)+6;
    yip1 = n_x*(i-1)+7;
    ydip1 = n_x*(i-1)+8;
    ui = n_x*N+n_u*(i-1)+1;
    hi = n_x*N+n_u*(N-1)+1;

    th_f = n_c + n_x*i-3;
    thd_f = n_c + n_x*i-2;
    y_f = n_c + n_x*i-1;
    yd_f= n_c + n_x*i;

    A((i-1)*8+1:i*8,:) = [thi,     th_f,  -1;
                          thdi,    thd_f, 0;
                          yi,      y_f,   -1;
                          ydi,     yd_f,  -1;
                          thip1,   th_f,   1;
                          thdip1,  thd_f,  1;
                          yip1,    y_f,    1;
                          ydip1,   yd_f,   1];
                      
                      
%     [ 0, v, 0, 0, 0, 0, 0, 0;
%       v, v, 0, 0, 0, 0, 0, 0;
%       0, 0, 0, v, 0, 0, 0, 0;
%       v, v, 0, 0, 0, 0, 0, 0];

    G((i-1)*12+1:i*12,:)= [thdi, th_f;
                           thi, thd_f;
                           thdi, thd_f;
                           ydi, y_f;
                           thi, yd_f;
                           thdi, yd_f;
                           ui, thd_f; 
                           ui, yd_f;
                           hi, th_f;
                           hi, thd_f;
                           hi, y_f;
                           hi, yd_f];
end

for i = 1:N*n_x
    G(end+1,:) = [i, 1];
end


Amat = []; %A(:,3);
iAfun = []; %A(:,2);
jAvar = []; %A(:,1);
iGfun = []; %G(:,2);
jGvar = []; %G(:,1);

end