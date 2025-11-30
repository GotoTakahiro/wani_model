%data_T_allの形式(:,113)
%絶対角を座標としたときの角度，角速度，角加速度，並進加速度を導出
%筋腱張力によってx,y方向にかかる力を導出
%その結果から関節間の接触力を計算(とりあえず並進の釣り合いを確認)
% clear;
 close all;
% clearvars
load('/Users/goto/Documents/Matlab_goto/crocodile_sim_PID-main/results/noGE/exp20251028noGE_CFL350_Ci44_CFLT100_GEo35_GE185_2.mat');
% load('results/20240712_MuscleLengthTest_PID/exp20240712_MuscleLengthTest_PID_1125_P50000_I50_D550_CFL350_Ci44_CFLT107_GEo43_GE200.mat');
% load('results/20240726_init_condition_test_per2mm/exp20240726_init_condition_test_per2mm_223_Hip20_Knee44_CFL350_Ci44_CFLT107_GEo37_GE188.mat');
% load('results/20240822_for_nolta_paper_rev/knee83/exp20240822_for_nolta_paper_noGE_knee83_CFL350_Ci44_CFLT100_GEo35_GE185.mat')
% load('results/20240718_HighWalk_init/exp20240718_HighWalk_init_PID_885_P50000_I50_D550_CFL350_Ci44_CFLT98_GEo35_GE197.mat')
% load('results/20250114/exp20250114_NoltaInit_test_noPull_noCFLT_2_CFL350_Ci44_CFLT100_GEo35_GE185.mat')
[filepath,name,ext] = fileparts(filename);
% save_path = 'results/20240712_MuscleLengthTest_PID/20250114_for_check/';
% save_path = 'results/20250114/';
save_path = 'results/';
% save_path = '';
new_filename = fullfile([save_path name '.mp4']);

% プロットの設定．
graph_save = true;
graph_view = true;

time_lim = max(t(:,1));

noCFLT = false;
noGE = false;
%読み込んだファイル名に含まれていた場合，対応したリンクの描画をしない
if contains(name,'noCFLT')
    noCFLT = true;
end
if contains(name,'noGE')
    noGE = true;
end

phi = linspace(0,2*pi,100);
r = l_link_list(7);

% プロット用に座標を計算
[coordinates_x, coordinates_y, angle_wire_pulley] = calc_coordinate_for_plot(q,r,l_link_list); %1:hip, 2:4th trochanter, 3:GE origin, 4:knee, 5:ankle, 6:toem 7:CFTLT branch, 8:Y shaped branch
%データ格納用の変数の定義
muscle_tension = zeros(size(q,1),4);
torque_muscle = zeros(size(q,1),10);
torque_all = zeros(size(q,1),10);
COM = zeros(size(q,1),2);
T_all = zeros(size(q,1),10);

momentum_frame = zeros(size(q,1),2);
momentum_femur = zeros(size(q,1),2);
momentum_tibia = zeros(size(q,1),2);
momentum_metatarsal = zeros(size(q,1),2);

momentum_list = zeros(size(q,1),4);
angular_moment_list = zeros(size(q,1),4);
momentum_COM_list = zeros(size(q,1),8);

tension_flag = false;

for j = 1:size(q,1)
    %qの1~10までを抜き出し
    general_q = q(j,1:10).';
    %重心座標を求める
    COM(j,:) = (calc_COM(m_list,l_link_list,general_q))';
end

phi=zeros(length(q(:,1)),4);
dphi=zeros(length(q(:,1)),4);
ddphi=zeros(length(q(:,1)),4);

% 仕事の計算
power = q(50:end,20).*data_Q(50:end,11);
work = trapz(t(50:end,1),power);
disp(['Work: ', num2str(work)]);
%絶対角の導出
phi(:,1)=q(:,5);
phi(:,2)=q(:,5)+q(:,6);
phi(:,3)=q(:,5)+q(:,6)+q(:,7);
phi(:,4)=q(:,5)+q(:,6)+q(:,7)+q(:,8);

dphi(:,1)=q(:,15);
dphi(:,2)=q(:,15)+q(:,16);
dphi(:,3)=q(:,15)+q(:,16)+q(:,17);
dphi(:,4)=q(:,15)+q(:,16)+q(:,17)+q(:,18);

ddphi(:,1)=data_accel(:,6);
ddphi(:,2)=data_accel(:,6)+data_accel(:,7);
ddphi(:,3)=data_accel(:,6)+data_accel(:,7)+data_accel(:,8);
ddphi(:,4)=data_accel(:,6)+data_accel(:,7)+data_accel(:,8)+data_accel(:,9);

%重力項
%重力加速度
g = 9.81;
%F_ankle_M2= func_F_ankle(q,l_link_list,data_T_gravity_all(:,1:10));
%F_ankle_M3= func_F_ankle(q,l_link_list,data_T_gravity_all(:,11:20));
F_frame= (m_list(4)+m_list(5))*g;
F_fem= m_list(6)*g;
F_tib= m_list(7)*g;
F_met= m_list(10)*g;
F_gravity=F_frame+F_fem+F_tib+F_met;

%各重心の並進加速度を計算
accel_CoM = func_accel(phi,dphi,ddphi,m_list,l_link_list);

%各筋腱張力のxy方向成分を計算
F_dist = func_F_dist(q,r,l_link_list,data_F_all);

%各外力を計算
%[hip_x,hip_y,heel_x,heel_y,toe_x,toe_y]
GRF=[data_force(:,2)+data_force(:,4),data_force(:,3)+data_force(:,5),data_force(:,6),data_force(:,7),data_force(:,8),data_force(:,9)];

%各関節の関節力を計算
F_connect = func_F_connect(q,r,l_link_list,m_list,accel_CoM,F_dist,GRF);
M_connect = func_M_connect(q,r,l_link_list,m_list,accel_CoM,F_dist,GRF,data_F_all,ddphi,data_T_all);


%外力による力の和を計算
GRF_x=data_force(:,2)+data_force(:,4)+data_force(:,6)+data_force(:,8);
GRF_y=data_force(:,3)+data_force(:,5)+data_force(:,7)+data_force(:,9);

if graph_view == true

    %Figure1：筋腱張力によってx方向にかかる力を表示
    figure(1)
    plot(t(:,1),F_connect(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),F_connect(:,3),'LineWidth',2);
    plot(t(:,1),F_connect(:,5),'LineWidth',2);
    plot(t(:,1),F_connect(:,7),'LineWidth',2);
    %plot(t(:,1),F_dist(:,1)+F_dist(:,3)+F_dist(:,5)+F_dist(:,7),'LineWidth',2); 
    %plot(t(:,1),GRF_x,'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-15 15]);
    legend('$ankle$','$knee$','$hip$','$frame$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_connectfx[1].pdf'], 'ContentType', 'vector');
    end

    %Figure2：筋腱張力によってy方向にかかる力を表示
    figure(2)
    plot(t(:,1),F_connect(:,2),'LineWidth',2);
    hold on
    plot(t(:,1),F_connect(:,4),'LineWidth',2);
    plot(t(:,1),F_connect(:,6),'LineWidth',2);
    plot(t(:,1),F_connect(:,8),'LineWidth',2);
    %plot(t(:,1),F_dist(:,2)+F_dist(:,4)+F_dist(:,6)+F_dist(:,8),'LineWidth',2);    
    %plot(t(:,1),F_gravity(:,1),'LineWidth',2);
    %plot(t(:,1),GRF_y,'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    ylim([-120 40]);
    legend('$ankle$','$knee$','$hip$','$frame$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_connectfy[1].pdf'], 'ContentType', 'vector');
    end

end
