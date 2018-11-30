function [A,B] = cart_pole_error_dyn_lin()
syms theta thetadot theta y ydot h w u g Mp Mc L theta_D thetadot_D y_D ydot_D u_D mu_k

x = [theta; thetadot; y; ydot];
x_D = [theta_D; thetadot_D; y_D; ydot_D];
e = x - x_D;
v = u - u_D;

e_dot = f(x, u) - f(x_D, u_D);

    
    function xdot = f(x, u)
        syms g Mp Mc L mu_k
        theta1 = x(1);
        thetadot1 = x(2);
        y1 = x(3);
        ydot1 = x(4);
        
        thetaddot1 = (-(Mp+Mc)*g*sin(theta1)-cos(theta1)*(u+Mp*L*thetadot1^2*sin(theta1)))/...
            (L*(Mp+Mc)-Mp*L*cos(theta1)*cos(theta1));
        Nc = 0; S = 0;
        yddot1 = (u-mu_k*Nc*S+Mp*L*(thetadot1^2*sin(theta1)-thetaddot1*cos(theta1)))/(Mc+Mp);
        xdot = [thetadot1; thetaddot1; ydot1; yddot1];
    end

end
