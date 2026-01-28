%data_T_allの形式(:,113)
%絶対角を座標としたときの角度，角速度，角加速度，並進加速度を導出
%筋腱張力によってx,y方向にかかる力を導出
%その結果から関節間の関節間力，モーメントを計算
% clear;
 close all;
% clearvars

addpath('/Users/goto/Documents/Matlab/crocodile_sim');

load('/Users/goto/Documents/MATLAB/crocodile_sim/plot/results/noGE/exp20251028noGE_CFL350_Ci44_CFLT100_GEo35_GE185.mat');
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
M_connect = func_M_connect(q,r,l_link_list,m_list,accel_CoM,F_dist,GRF,data_F_all,ddphi,data_T_all,F_connect,data_accel);
%各筋腱の長さを計算
musclelength = func_musclelength(q,r,l_link_list);


if graph_view == true

    %Figure1：関節力のx方向成分の力を表示
    figure(1)
    plot(t(:,1),F_connect(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),F_connect(:,3),'LineWidth',2);
    plot(t(:,1),F_connect(:,5),'LineWidth',2);
    %plot(t(:,1),F_connect(:,7),'LineWidth',2);
    %plot(t(:,1),F_dist(:,1)+F_dist(:,3)+F_dist(:,5)+F_dist(:,7),'LineWidth',2); 
    %plot(t(:,1),GRF_x,'LineWidth',2);
    hold off
    ax = gca; % 現在の軸を取得
    ax.FontSize = 24;         % 軸の数値（目盛）のサイズ
    ax.FontName = 'Times New Roman'; % フォントの種類を統一（任意）
    xlim([0 time_lim]);
    ylim([-40 140]);
    legend('ankle','knee','hip','M0','Interpreter', 'latex', 'FontSize',20,'Location','northeast');
    xlabel('Time [s]', 'FontSize', 28);
    ylabel('Force [N]', 'FontSize', 28);
    if graph_save == true
        exportgraphics(gca, [save_path name '_connectfx[1].pdf'], 'ContentType', 'vector');
    end

    %Figure2：関節力のy方向成分の力を表示
    figure(2)
    plot(t(:,1),F_connect(:,2),'LineWidth',2);
    hold on
    plot(t(:,1),F_connect(:,4),'LineWidth',2);
    plot(t(:,1),F_connect(:,6),'LineWidth',2);
    %plot(t(:,1),F_connect(:,8),'LineWidth',2);
    %plot(t(:,1),F_dist(:,2)+F_dist(:,4)+F_dist(:,6)+F_dist(:,8),'LineWidth',2);    
    %plot(t(:,1),F_gravity(:,1),'LineWidth',2);
    %plot(t(:,1),GRF_y,'LineWidth',2);
    hold off
    ax = gca; % 現在の軸を取得
    ax.FontSize = 24;         % 軸の数値（目盛）のサイズ
    ax.FontName = 'Times New Roman'; % フォントの種類を統一（任意）

    xlim([0 time_lim]);
    ylim([-120 40]);
    legend('ankle','knee','hip','M0','Interpreter', 'latex', 'FontSize',20,'Location','northeast');
    xlabel('Time [s]', 'FontSize', 28);
    ylabel('Force [N]', 'FontSize', 28);
    if graph_save == true
        exportgraphics(gca, [save_path name '_connectfy[1].pdf'], 'ContentType', 'vector');
    end

    %Figure3：関節力のモーメントを表示
    figure(3)
    plot(t(:,1),M_connect(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),M_connect(:,2),'LineWidth',2);
    plot(t(:,1),M_connect(:,3),'LineWidth',2);
    %plot(t(:,1),M_connect(:,4),'LineWidth',2);
    %plot(t(:,1),F_dist(:,2)+F_dist(:,4)+F_dist(:,6)+F_dist(:,8),'LineWidth',2);    
    %plot(t(:,1),F_gravity(:,1),'LineWidth',2);
    %plot(t(:,1),GRF_y,'LineWidth',2);
    hold off
    ax = gca; % 現在の軸を取得
    ax.FontSize = 24;         % 軸の数値（目盛）のサイズ
    ax.FontName = 'Times New Roman'; % フォントの種類を統一（任意）

    xlim([0 time_lim]);
    ylim([-1 1]);
    legend('ankle','knee','hip','M0','Interpreter', 'latex', 'FontSize',20,'Location','northeast');
    xlabel('Time [s]', 'FontSize', 28);
    ylabel('Torque [Nm]', 'FontSize', 28);
    if graph_save == true
        exportgraphics(gca, [save_path name '_connectmom[1].pdf'], 'ContentType', 'vector');
    end

    %Figure4：ankle周りのモーメントを表示
    figure(4)
    plot(t(:,1),M_connect(:,8),'LineWidth',2);
    hold on
    plot(t(:,1),M_connect(:,9),'LineWidth',2);
    plot(t(:,1),M_connect(:,10),'LineWidth',2);
    plot(t(:,1),M_connect(:,5),'LineWidth',2);
    plot(t(:,1),M_connect(:,6),'LineWidth',2);
    plot(t(:,1),data_T_all(:,93),'LineWidth',2);
    plot(t(:,1),M_connect(:,7),'LineWidth',2);    
    %plot(t(:,1),M_connect(:,8)+M_connect(:,9)+M_connect(:,10)+M_connect(:,5)+M_connect(:,6)+data_T_all(:,93)+M_connect(:,7),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-120 40]);
    legend('$grav$','$int_x$','$int_y$','$N_x$','$N_y$','$T_{ankleex}$','$FGE$','$All$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Torque [Nm]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_momankle[1].pdf'], 'ContentType', 'vector');
    end
    % %Figure5：knee周りのモーメントを表示
    % figure(5)
    % plot(t(:,1),M_connect(:,11),'LineWidth',2);
    % hold on
    % plot(t(:,1),M_connect(:,12),'LineWidth',2);
    % plot(t(:,1),M_connect(:,13),'LineWidth',2);
    % plot(t(:,1),M_connect(:,14),'LineWidth',2);
    % plot(t(:,1),M_connect(:,15),'LineWidth',2);
    % %plot(t(:,1),M_connect(:,11)+M_connect(:,12)+M_connect(:,13)+M_connect(:,14)+M_connect(:,15),'LineWidth',2);
    % %plot(t(:,1),GRF_y,'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % %ylim([-120 40]);
    % legend('$joint_x$','$joint_y$','$grav$','$int_x$','$int_y$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Time [s]');
    % ylabel('Torque [Nm]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_momknee[1].pdf'], 'ContentType', 'vector');
    % end
    % %Figure6：hip周りのモーメントを表示
    % figure(6)
    % plot(t(:,1),M_connect(:,16),'LineWidth',2);
    % hold on
    % plot(t(:,1),M_connect(:,17),'LineWidth',2);
    % plot(t(:,1),M_connect(:,18),'LineWidth',2);
    % plot(t(:,1),M_connect(:,19),'LineWidth',2);
    % plot(t(:,1),M_connect(:,20),'LineWidth',2);
    % %plot(t(:,1),GRF_y,'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % ylim([-2 12]);
    % legend('$joint_x$','$joint_y$','$grav$','$int_x$','$int_y$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Time [s]');
    % ylabel('Torque [Nm]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_momhip[1].pdf'], 'ContentType', 'vector');
    % end
    % 
    % %Figure7：hip周りのモーメントを表示2
    % figure(7)
    % plot(t(:,1),M_connect(:,21),'LineWidth',2);
    % hold on
    % plot(t(:,1),M_connect(:,22),'LineWidth',2);
    % plot(t(:,1),M_connect(:,23),'LineWidth',2);
    % plot(t(:,1),M_connect(:,24),'LineWidth',2);
    % %plot(t(:,1),GRF_y,'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % ylim([-12 2]);
    % legend('$FCi_x$','$FCi_y$','$FGEo_x$','$FGEo_y$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Time [s]');
    % ylabel('Torque [Nm]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_momhip[2].pdf'], 'ContentType', 'vector');
    % end
    % %Figure8：body周りのモーメントを表示
    % figure(8)
    % plot(t(:,1),M_connect(:,25),'LineWidth',2);
    % 
    % hold on
    % plot(t(:,1),M_connect(:,26),'LineWidth',2);
    % plot(t(:,1),M_connect(:,27),'LineWidth',2);
    % plot(t(:,1),M_connect(:,28),'LineWidth',2);
    % plot(t(:,1),M_connect(:,29),'LineWidth',2);
    % plot(t(:,1),data_T_all(:,91),'LineWidth',2);
    % 
    % %plot(t(:,1),GRF_y,'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % %ylim([-120 40]);
    % legend('$joint_x$','$joint_y$','$grav$','$int_x$','$int_y$','$T_{frame}$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Time [s]');
    % ylabel('Torque [Nm]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_momframe[1].pdf'], 'ContentType', 'vector');
    % end
    % 
    % %Figure9：body周りのモーメントを表示
    % figure(9)
    % axes1 = axes('Parent',gcf);
    % %plot(phi(:,2),musclelength(:,1),'LineWidth',2);
    % hold on
    % plot(phi(:,2),(musclelength(:,2)-l_muscle_list(1))*1000,'LineWidth',2);
    % plot(phi(:,2),(musclelength(:,3)-l_muscle_list(3))*1000,'LineWidth',2);
    % plot(phi(:,2),(musclelength(:,4)-l_muscle_list(2))*1000,'LineWidth',2);
    % %plot(phi(:,2),musclelength(:,5),'LineWidth',2);
    % 
    % %plot(t(:,1),GRF_y,'LineWidth',2);
    % hold off
    % 
    % xlim([-3.2 3.2]);
    % xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    % set(axes1,'FontName','Times New Roman','FontSize',26,'XGrid','off','XTick',[-pi -pi/2 0 pi/2 pi]);
    % ylim([-1 1]);
    % legend('$Ci$','$GEo$','$CFLT$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Angle [rad]');
    % ylabel('dLength [mm]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_musclelength[1].pdf'], 'ContentType', 'vector');
    % end
    % 
    % %Figure10：力による各関節周りのモーメントを表示
    % figure(10)
    % plot(t(:,1),M_connect(:,32),'LineWidth',2);
    % 
    % hold on
    % plot(t(:,1),M_connect(:,33),'LineWidth',2);
    % plot(t(:,1),M_connect(:,34),'LineWidth',2);
    % plot(t(:,1),M_connect(:,35),'LineWidth',2);
    % 
    % %plot(t(:,1),GRF_y,'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % %ylim([-120 40]);
    % legend('$ankle$','$knee$','$hip$','$body$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Time [s]');
    % ylabel('Torque [Nm]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_momforce[1].pdf'], 'ContentType', 'vector');
    % end
    % 
    % %Figure11：つま先周りのモーメントを表示
    % figure(11)
    % plot(t(:,1),M_connect(:,36),'LineWidth',2);
    % 
    % hold on
    % plot(t(:,1),M_connect(:,37),'LineWidth',2);
    % 
    % %plot(t(:,1),GRF_y,'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % ylim([-1 8]);
    % legend('$clock$','$unclock$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Time [s]');
    % ylabel('Torque [Nm]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_momtoe[1].pdf'], 'ContentType', 'vector');
    % end
    % %Figure12：つま先周りのモーメントを表示
    % figure(12)
    % plot(t(:,1),M_connect(:,38),'LineWidth',2);
    % 
    % hold on
    % plot(t(:,1),M_connect(:,39),'LineWidth',2);
    % 
    % %plot(t(:,1),GRF_y,'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % %ylim([-120 40]);
    % legend('$clock$','$unclock$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Time [s]');
    % ylabel('Torque [Nm]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_momtoe[2].pdf'], 'ContentType', 'vector');
    % end

    



end
