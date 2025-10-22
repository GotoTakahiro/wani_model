clear all

syms theta l m dtheta dl g
syms Pgain Igain Dgain
syms target_length
syms error int_error derror
syms k c x_sp y_sp l_sp

q = [theta; l];
dq = [dtheta; dl];

x = l*sin(theta);
y = -l*cos(theta);

dx = jacobian(x, q)*dq;
dy = jacobian(y, q)*dq;

T = 1/2*m*(dx^2 + dy^2);
U_g = m*g*y;
l_dis = ((x - x_sp)^2 + (y - y_sp)^2)^(1/2) - l_sp;
U_sp = 1/2*k*l_dis^2;
U = U_g + U_sp;

D = 1/2*c*(dtheta^2 + dl^2);

L = T - U;

M = jacobian(jacobian(L, dq), dq);
C = jacobian(jacobian(L, dq), q)*dq;
G = jacobian(L, q).';
D = jacobian(D, dq).';

F_l = Pgain*error + Igain*int_error;
F_theta = 0;
F = [F_theta; F_l];

matlabFunction(M, 'File', 'calc_M_pend');
matlabFunction(C, 'File', 'calc_C_pend');
matlabFunction(G, 'File', 'calc_G_pend');
matlabFunction(F, 'File', 'calc_F_pend');
matlabFunction(D, 'File', 'calc_D_pend');