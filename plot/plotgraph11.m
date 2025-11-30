%data_T_allの形式(:,113)
%絶対角を座標としたときの角度，角速度，角加速度，並進加速度を導出
% clear;
 close all;
% clearvars
addpath('/Users/goto/Documents/Matlab_goto/crocodile_sim_PID-main');

load('/Users/goto/Documents/Matlab_goto/crocodile_sim_PID-main/plot/exp20251028_CFL350_Ci44_CFLT100_GEo35_GE190.mat');
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
COM_link = zeros(size(q,1),2);
moment_gravity = zeros(size(q,1),1);
T_all = zeros(size(q,1),10);

momentum_frame = zeros(size(q,1),2);
momentum_femur = zeros(size(q,1),2);
momentum_tibia = zeros(size(q,1),2);
momentum_metatarsal = zeros(size(q,1),2);

momentum_list = zeros(size(q,1),4);
angular_moment_list = zeros(size(q,1),4);
momentum_COM_list = zeros(size(q,1),8);

tension_flag = false;
M_all=m_list(1)+m_list(2)+m_list(3)+m_list(4)+m_list(5)+m_list(6)+m_list(7)+m_list(8)+m_list(9);
g=9.81;
for j = 1:size(q,1)
    %qの1~10までを抜き出し
    general_q = q(j,1:10).';
    %リンクの重心座標を求める
    COM(j,:) = (calc_COM(m_list,l_link_list,general_q))';
    COM_link(j,:) = (calc_COM_link(m_list,l_link_list,general_q))';
    moment_gravity(j,1) = M_all*g*(COM_link(j,1)-COM(j,1));

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

%各重心の並進加速度を計算
accel_CoM = func_accel(phi,dphi,ddphi,m_list,l_link_list);

%総質量
M_link=m_list(4)+m_list(5)+m_list(6)+m_list(7)+m_list(8)+m_list(9);
moment_tale = func_moment2(q,m_list,l_link_list,data_accel_GRF,r,data_F_all);
%外力によるモーメントの和
moment_externalforce=moment_tale(:,1)+moment_tale(:,2)+moment_tale(:,3)+moment_tale(:,4)+moment_tale(:,5)+moment_tale(:,6);
%筋肉によるモーメントの和
moment_muscle=moment_tale(:,9)+moment_tale(:,10)+moment_tale(:,11)+moment_tale(:,12)+moment_tale(:,15)+moment_tale(:,16);


%一般化力から股関節にかかる力を計算

if graph_view == true

    %Figure9：CFL起始点周りのモーメントの表示
    figure(9)
    plot(t(:,1),moment_tale(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),moment_tale(:,2),'LineWidth',2);
    plot(t(:,1),moment_tale(:,3),'LineWidth',2);
    plot(t(:,1),moment_tale(:,4),'LineWidth',2);
    plot(t(:,1),moment_tale(:,5),'LineWidth',2);
    plot(t(:,1),moment_tale(:,6),'LineWidth',2);
    plot(t(:,1),moment_externalforce,'LineWidth',2);
    hold off
    %ylim([-4 6]);
    %xlim([0 time_lim]);
    %ylim([-2 2]);
    legend('$heel_x$','$heel_y$','$toe_x$','$toe_y$','$hip_x$','$hip_y$','$all$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Moment[Nm]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_moment[1].pdf'], 'ContentType', 'vector');
    end

    %Figure10：CFL起始点周りのモーメントの表示
    figure(10)
    hold on
    plot(t(:,1),moment_tale(:,9),'LineWidth',2);
    plot(t(:,1),moment_tale(:,10),'LineWidth',2);
    plot(t(:,1),moment_tale(:,11),'LineWidth',2);
    plot(t(:,1),moment_tale(:,12),'LineWidth',2);
    plot(t(:,1),moment_tale(:,13),'LineWidth',2);
    plot(t(:,1),moment_tale(:,14),'LineWidth',2);
    hold off
    %ylim([-4 6]);
    %xlim([0 time_lim]);
    %ylim([-2 2]);
    legend('$Ci_x$','$Ci_y$','$GEo_x$','$GEo_y$','$GE_x$','$GE_y$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Moment[Nm]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_moment[2].pdf'], 'ContentType', 'vector');
    end

    %Figure11：CFL起始点周りのモーメントの表示
    figure(11)
    plot(t(:,1),moment_externalforce,'LineWidth',2);
    hold on
    plot(t(:,1),moment_muscle,'LineWidth',2);
    %plot(t(:,1),moment_gravity,'LineWidth',2);
    plot(t(:,1),moment_externalforce+moment_muscle+moment_tale(:,8),'LineWidth',2); %
    hold off
    ylim([-40 40]);
    %xlim([0 time_lim]);
    %ylim([-2 2]);
    legend('$GRF$','$muscle$','$all$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Moment[Nm]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_moment[3].pdf'], 'ContentType', 'vector');
    end

end
