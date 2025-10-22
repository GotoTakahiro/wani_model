clear;
close all;
clearvars

global force_list now_force dt data_dt data_prev_t
force_list = [];
force_list(1) = 0;
now_force = 0;
data_dt = [];
data_prev_t = [];

m = 1;
l = 1;
g = 9.81;
theta = pi/4;
dtheta = 0;
dl = 0;

l_sp = 0.15;
x_sp = l*sin(theta);
y_sp = -l*cos(theta)+l_sp;

default_l = 1;
l_target = 1;
error = l_target - l;
prev_error = 0;

tmax = 10;
tspace = 0.01;
tspan = 0:tspace:tmax;

Pgain = 5000;
Igain = 10;
% Dgain = 0.1;
c = 0;
k = 5000;
dt = 0;
prev_t = 0;
error = 0;
prev_error = 0;

q = [theta l dtheta dl error prev_error prev_t];
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'OutputFcn', @get_force);

f = @(t, q) calc_EOM(t,q,Pgain,Igain,m,g,default_l,l_sp,x_sp,y_sp,k,c);

[t, q] = ode45(f, tspan, q, options);

power = zeros(size(force_list));
for i = 1:length(force_list)
    power(i) = force_list(i) * q(i,4);
end

% figure(1)
% plot(t, power);
% xlabel('Time [s]');
% ylabel('Power [W]');

% %台形積分で仕事を求める
% work = trapz(t, power);
% disp(work);

% new_filename = 'output.mp4';
% v = VideoWriter(new_filename,'MPEG-4');
% v.Quality = 100;
% v.FrameRate = 100;
% open(v)
% ax = gca();

% for i = 1:size(q,1)
%     theta = q(i,1);
%     l = q(i,2);
%     x = l*sin(theta);
%     y = -l*cos(theta);
%     plot(ax, [0 x], [0 y], 'r', 'LineWidth', 2);
%     hold on;
%     plot(ax, [x x_sp], [y y_sp], 'b', 'LineWidth', 2);
%     axis(ax, 'equal');
%     axis(ax, [-1.5 1.5 -1.5 1.5]);
%     drawnow;
%     hold off
%     frame = getframe(ax);
%     writeVideo(v,frame);
% end
% close(v)


function dq = calc_EOM(t,q,Pgain,Igain,m,g,default_l,l_sp,x_sp,y_sp,k,c)
    global now_force dt prev_t

    % 状態変数の取得
    l_target = default_l - 0.05*t;
    theta = q(1);
    l = q(2);
    dtheta = q(3);
    dl = q(4);
    int_error = q(5);
    previous_error = q(6);

    prev_t = q(7);
    % dt = t - prev_t;
    dt = max(t - prev_t, 1e-5);
    % disp(dt)

    % 誤差の計算
    % disp(l_target);
    error_ = l_target - l;
    drror = (error_ - previous_error)/dt;
    
    
    M = calc_M_pend(l, m, theta);
    C = calc_C_pend(dl, dtheta, l, m, theta);
    G = calc_G_pend(dl,dtheta,g,k,l,l_sp,m,theta,x_sp,y_sp);
    D = calc_D_pend(c,dl,dtheta);
    F = calc_F_pend(Igain,Pgain,error_,int_error);
    % if error_ < 0
    %     disp(F);
    % end
    now_force = F(2);

    prev_error = error_;

    dq = [q(3); q(4); M\(-C+G+D+F); prev_error; drror; 1];
end

function status = get_force(t,q,flag)
    global force_list now_force dt data_dt prev_t data_prev_t
    status = 0;
    if isempty(flag) || strcmp(flag,'')
        force_list(end+1) = now_force;
        data_dt(end+1) = dt;
        data_prev_t(end+1) = prev_t;
        disp(t);
    end
end
