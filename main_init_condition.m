%ワニモデルシミュレーション
clear;
close all;
clearvars

%重力加速度
g = 9.81;

%初期値設定

%リンクの長さを定義
L_fem = 160.563/1000;
L_tib = 166.190/1000;
L_met = 90.739/1000;
L_frame = 0.30;
L_4th_troch = 49.597/1000;
L_GE_origin = 152.265/1000;
r = 0.015; %プーリ半径
l_link_list = zeros(7,1);
l_link_list(1) = L_fem;
l_link_list(2) = L_tib;
l_link_list(3) = L_met;
l_link_list(4) = L_4th_troch;
l_link_list(5) = L_GE_origin;
l_link_list(6) = L_frame;
l_link_list(7) = r;

%各リンクと質点の質量
M1 = 0.01;
M2 = 0.01;
M3 = 0.01;
M_hip = 2.0;
M_frame = 0.5;
M_fem = 0.2;
M_tib = 0.2;
M_met_pulley = 0.3;
M_met_beam = 0.1;
M_met = M_met_pulley + M_met_beam;
m_list = zeros(10,1);
m_list(1) = M1;
m_list(2) = M2;
m_list(3) = M3;
m_list(4) = M_hip;
m_list(5) = M_frame;
m_list(6) = M_fem;
m_list(7) = M_tib;
m_list(8) = M_met_pulley;
m_list(9) = M_met_beam;
m_list(10) = M_met;

%ばね定数
default_wire_k = 100000; %筋腱
k_frame = 100000; %caudal vertebral columnを水平に保つばね
k_ground = 50000; %地面の硬さ

%筋腱減衰係数
default_wire_c = 200; %筋腱
c_frame = 50; %caudal vertebral column
c_ground = 100; %地面の

%摩擦係数
mu = 4;

%Ankle joint limit
limit_ankle_ex = 120/180*pi;
% limit_ankle_ex = 140/180*pi;
limit_ankle_flex = 10/180*pi;
ankle_limit = zeros(2,1);
ankle_limit(1) = limit_ankle_ex;
ankle_limit(2) = limit_ankle_flex;

default_frame_angle = 90/180*pi;
limit_list = zeros(3,1);
limit_list(1) = limit_ankle_ex;
limit_list(2) = limit_ankle_flex;
limit_list(3) = default_frame_angle;

%シミュレーション
tmax = 13;
tspace = 0.01;
tspan = 0:tspace:tmax; % シミュレーションの時間範囲
t_CFL = 1; %CFL収縮開始
end_CFL = 9; %CFL収縮量減衰開始
t_end_exp = 10; %CFLの収縮終了
CFL_alpha = 0.01; %単位時間あたりのCFLの収縮量


%ワイヤ長さ
L_CFL = 0.35;
% L_Ci = [0.036 0.038 0.04 0.042 0.044];
L_Ci = 0.044;
% L_CFLT = [0.095 0.098 0.101 0.104 0.107];
L_CFLT = 0.098;
% L_GEo = [0.035 0.037 0.039 0.041 0.043];
L_GEo = 0.035;
% L_GE = [0.176 0.179 0.182 0.185 0.188 0.191 0.194 0.197 0.200];
% L_GE = 0.180;
L_GE = 0.197;

%ゲイン設定
Pgain_list = 50000;
Igain_list = 50;
Dgain_list = 550;

%組み合わせを保存
[CFL, Ci, CFLT, GEo, GE, Pgain, Igain, Dgain] = ndgrid(L_CFL, L_Ci, L_CFLT, L_GEo, L_GE, Pgain_list, Igain_list, Dgain_list);
length_and_gain_combination = [CFL(:), Ci(:), CFLT(:), GEo(:), GE(:), Pgain(:), Igain(:), Dgain(:)];
save('exp20240726_init_condition_test_MinWork_per2mm_length_and_gain_combination.mat','length_and_gain_combination');

% load('exp20240712_MuscleLengthTest_PID_length_and_gain_combination.mat');
% one_therd_num = 375;
% max_num = size(length_and_gain_combination,1);

%初期条件設定．
init_theta2_ = 10:2:70; %股関節
init_theta2_ = -init_theta2_/180*pi;
init_theta3_ = 30:2:100; %膝関節
init_theta3_ = -init_theta3_/180*pi;
[init_theta2, init_theta3] = ndgrid(init_theta2_, init_theta3_);
init_condition_list = [init_theta2(:), init_theta3(:)];
disp(size(init_condition_list))

for i = size(init_condition_list,1):-1:1
    theta2_ = init_condition_list(i,1);
    theta3_ = init_condition_list(i,2);
    theta4_ = 90/180*pi-(90/180*pi+theta2_+theta3_)-(asin(r/L_met)-0.03); %足関節
    if theta4_ > limit_ankle_ex %初期の足関節角度がlimit_ankle_exよりも大きくなる時は除外
        init_condition_list(i,:) = [];
    end
end
disp(size(init_condition_list))
% 初期姿勢の組み合わせを保存
save('exp20240723_HighWalk_init_condition.mat','init_condition_list');

disp(['Length and gain combination: ', size(length_and_gain_combination,1)]);
disp(['Initial condition: ', size(init_condition_list,1)]);

sim_num = size(length_and_gain_combination,1)*size(init_condition_list,1);
disp(['Number of simulations: ', num2str(sim_num)]);


tic
for i = 1:size(length_and_gain_combination,1)
    L_CFL = length_and_gain_combination(i,1);
    L_Ci = length_and_gain_combination(i,2);
    L_CFLT = length_and_gain_combination(i,3);
    L_GEo = length_and_gain_combination(i,4);
    L_GE = length_and_gain_combination(i,5);
    l_muscle_list = zeros(4,1);
    l_muscle_list(1) = L_Ci;
    l_muscle_list(2) = L_CFLT;
    l_muscle_list(3) = L_GEo;
    l_muscle_list(4) = L_GE;
    default_CFL = L_CFL;
    error_CFL = 0;
    Pgain = length_and_gain_combination(i,6);
    Igain = length_and_gain_combination(i,7);
    Dgain = length_and_gain_combination(i,8);
    gain_list = [Pgain; Igain; Dgain];
    
    for j = 1:sim_num
        disp(['Simulation ', num2str(i), '-', num2str(j), ' is start.']);
        %初期値設定
        theta1 = 90/180*pi;
        theta2 = init_condition_list(j,1);
        theta3 = init_condition_list(j,2);
        theta4 = 90/180*pi-(theta1+theta2+theta3)-(asin(r/L_met)-0.03);
        theta_CFL = 85/180*pi+theta2/5;
        L_CFL = 0.35;

        x1 = -0.30;
        y1 = 0;

        x2 = x1 + L_CFL*sin(theta_CFL);
        y2 = y1 - L_CFL*cos(theta_CFL);
        x3 = x2 + L_CFLT*sin(theta1+theta2);
        y3 = y2 - L_CFLT*cos(theta1+theta2);
        
        dx1 = 0;
        dy1 = 0;
        dx3 = 0;
        dy3 = 0;
        dtheta1 = 0;
        dtheta2 = 0;
        dtheta3 = 0;
        dtheta4 = 0;
        dtheta_CFL = 0;
        dl_CFL = 0;
        
        y_heel_init = calc_y_heel(l_link_list, [x1; y1; x3; y3; theta1; theta2; theta3; theta4; theta_CFL; L_CFL]);

        %固定点の座標
        x_fixed = zeros(1,1); %frame_stopper
        y_fixed = zeros(2,1); %y_stopper, ground

        x_fixed(1) = 0.0;
        y_fixed(1) = y1-0.005;
        y_fixed(2) = y_heel_init - 0.002;

        initial_condition = [x1 y1 x3 y3 theta1 theta2 theta3 theta4 theta_CFL L_CFL dx1 dy1 dx3 dy3 dtheta1 dtheta2 dtheta3 dtheta4 dtheta_CFL dl_CFL error_CFL];

        disp(['Hip: ', num2str(rad2deg(theta2)), ', Knee: ', num2str(rad2deg(theta3))]);

        filename = sprintf('exp20240723_init_condition_test_%d_Hip%d_Knee%d_CFL%d_Ci%d_CFLT%d_GEo%d_GE%d.mat',j,-rad2deg(theta2),-rad2deg(theta3),L_CFL*1000,L_Ci*1000,L_CFLT*1000,L_GEo*1000,L_GE*1000);

        simulation = solve_EOM(tmax,tspace,tspan,initial_condition,x_fixed,y_fixed,k_ground,c_ground,mu,l_link_list,l_muscle_list,limit_list,m_list,default_wire_k,default_wire_c,g,t_CFL,k_frame,c_frame,default_frame_angle,filename,end_CFL,CFL_alpha,t_end_exp,default_CFL,gain_list);
        disp(['Simulation ', num2str(i), '-', num2str(j), ' is done.']);
    end
end
toc