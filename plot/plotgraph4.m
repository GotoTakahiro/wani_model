%data_T_allの形式(:,113)
%一般化力からhip，COMにかかるx,y方向の力を計算しプロット（Ci,CFLT,GEoによる力を統合）
 clear;
 close all;
% clearvars
addpath('/Users/goto/Documents/Matlab/crocodile_sim_PID-main');
load('/Users/goto/Documents/MATLAB/crocodile_sim_PID-main/plot/results/m5/exp20251028_CFL350_Ci44_CFLT100_GEo35_GE185.mat');% load('results/20240712_MuscleLengthTest_PID/exp20240712_MuscleLengthTest_PID_1125_P50000_I50_D550_CFL350_Ci44_CFLT107_GEo43_GE200.mat');
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

for j = 1:size(q,1)
    %qの1~10までを抜き出し
    general_q = q(j,1:10).';
    %重心座標を求める
    COM(j,:) = (calc_COM(m_list,l_link_list,general_q))';
end



% 仕事の計算
power = q(50:end,20).*data_Q(50:end,11);
work = trapz(t(50:end,1),power);
disp(['Work: ', num2str(work)]);

% プロットの色を固定
CFL_Color = '#0072BD';
Ci_Color = '#D95319';
CFLT_Color = '#EDB120';
% Ci_Color = '#EDB120';    %lab meeting
% CFLT_Color = '#D95319';
GEo_Color = '#7E2F8E';
GE_Color = '#77AC30';


%全ての力を一般化力として扱った時の総トルクを計算
torque_all=(data_T_all(:,1:10)+data_T_all(:,11:20)+data_T_all(:,21:30)+data_T_all(:,31:40)+data_T_all(:,41:50)+data_T_all(:,51:60)+data_T_all(:,71:80));
torque_all=torque_all+data_T_gravity_all(:,1:10)+data_T_gravity_all(:,11:20)+data_T_gravity_all(:,21:30)+data_T_gravity_all(:,31:40)+data_T_gravity_all(:,41:50)+data_T_gravity_all(:,51:60);
torque_all(:,5)=torque_all(:,5)+data_T_all(:,91);
torque_all(:,8)=torque_all(:,8)+data_T_all(:,93)+data_T_all(:,92);
torque_all=torque_all+data_Q_hip_pull(:,1:10)+data_Q_heel(:,1:10)+data_Q_toe(:,1:10)+data_Q_hip_up(:,1:10);
%慣性行列の逆行列を掛け合わせたものを計算
for i = 1:size(q,1)
    M=Inertial_matrix(m_list,l_link_list,transpose(q(i,1:10)));
    accel=M\transpose(torque_all(i,1:10));
    T_all(i,1:10) =transpose(accel);
end

%一般化力からhipの質点に働く力を逆算するための関数
%ワイヤー項
F_hip_Ci_4th= func_F_hip(q,l_link_list,data_T_all(:,1:10));
F_hip_Ci_M2= func_F_hip(q,l_link_list,data_T_all(:,11:20));
F_hip_Ci=F_hip_Ci_4th+F_hip_Ci_M2;
F_hip_CFLT_M2= func_F_hip(q,l_link_list,data_T_all(:,21:30));
F_hip_CFLT_M3= func_F_hip(q,l_link_list,data_T_all(:,31:40)); %0
F_hip_CFLT=F_hip_CFLT_M2+F_hip_CFLT_M3;
F_hip_GEo_M3= func_F_hip(q,l_link_list,data_T_all(:,41:50));
F_hip_GEo_GEorigin= func_F_hip(q,l_link_list,data_T_all(:,51:60));
F_hip_GEo=F_hip_GEo_M3+F_hip_GEo_GEorigin;
%F_hip_GE_M3= func_F_hip(q,l_link_list,data_T_all(:,61:70));
F_hip_GE= func_F_hip(q,l_link_list,data_T_all(:,71:80));
F_hip_GE_pulley= func_F_hip(q,l_link_list,data_T_all(:,81:90));
F_hip_CFL_M1= func_F_hip(q,l_link_list,data_T_all(:,94:103));
F_hip_CFL_M2= func_F_hip(q,l_link_list,data_T_all(:,104:113));
F_hip_CFL=F_hip_CFL_M1+F_hip_CFL_M2;
%角度拘束による力
Q_frame=[zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),data_T_all(:,91),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1)];
Q_ex=[zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),data_T_all(:,93),zeros(size(q,1),1),zeros(size(q,1),1)];
F_hip_frame_fix= func_F_hip(q,l_link_list,Q_frame);
F_hip_ex= func_F_hip(q,l_link_list,Q_ex);
%重力項
F_hip_M2= func_F_hip(q,l_link_list,data_T_gravity_all(:,1:10));
F_hip_M3= func_F_hip(q,l_link_list,data_T_gravity_all(:,11:20));
F_hip_frame= func_F_hip(q,l_link_list,data_T_gravity_all(:,21:30));
F_hip_fem= func_F_hip(q,l_link_list,data_T_gravity_all(:,31:40));
F_hip_tib= func_F_hip(q,l_link_list,data_T_gravity_all(:,41:50));
F_hip_met= func_F_hip(q,l_link_list,data_T_gravity_all(:,51:60));

%外力
F_hip_hipup= func_F_hip(q,l_link_list,data_Q_hip_up(:,1:10));
F_hip_hippull= func_F_hip(q,l_link_list,data_Q_hip_pull(:,1:10));
F_hip_hipheel= func_F_hip(q,l_link_list,data_Q_heel(:,1:10));
F_hip_hiptoe= func_F_hip(q,l_link_list,data_Q_toe(:,1:10));
%股関節の支持力を除いた総和を計算
F_hip_muscletendon=F_hip_Ci_4th+F_hip_Ci_M2+F_hip_CFLT_M2+F_hip_CFLT_M3+F_hip_GEo_M3+F_hip_GEo_GEorigin+F_hip_GE+F_hip_frame_fix+F_hip_ex+F_hip_CFL_M1+F_hip_CFL_M2;
F_hip_gravity=F_hip_M2+F_hip_M3+F_hip_frame+F_hip_fem+F_hip_tib+F_hip_met;
F_hip_externalforce=F_hip_hipheel+F_hip_hiptoe;%+F_hip_hippull;
F_hip_all=F_hip_muscletendon+F_hip_gravity+F_hip_externalforce;%+F_hip_hipup;

%一般化力からCOMの質点に働く力を逆算するための関数
%ワイヤー項
F_COM_Ci_4th= func_F_COM(q,m_list,l_link_list,data_T_all(:,1:10));
F_COM_Ci_M2= func_F_COM(q,m_list,l_link_list,data_T_all(:,11:20));
F_COM_Ci=F_COM_Ci_4th+F_COM_Ci_M2;
F_COM_CFLT_M2= func_F_COM(q,m_list,l_link_list,data_T_all(:,21:30));
F_COM_CFLT_M3= func_F_COM(q,m_list,l_link_list,data_T_all(:,31:40)); %0
F_COM_CFLT=F_COM_CFLT_M2+F_COM_CFLT_M3;
F_COM_GEo_M3= func_F_COM(q,m_list,l_link_list,data_T_all(:,41:50));
F_COM_GEo_GEorigin= func_F_COM(q,m_list,l_link_list,data_T_all(:,51:60));
F_COM_GEo=F_COM_GEo_M3+F_COM_GEo_GEorigin;
%F_COM_GE_M3= func_F_COM(q,m_list,l_link_list,data_T_all(:,61:70));
F_COM_GE= func_F_COM(q,m_list,l_link_list,data_T_all(:,71:80));
F_COM_GE_pulley= func_F_COM(q,m_list,l_link_list,data_T_all(:,81:90));
F_COM_frame_fix= func_F_COM(q,m_list,l_link_list,Q_frame(:,1:10));
F_COM_ex= func_F_COM(q,m_list,l_link_list,Q_ex(:,1:10));
F_COM_CFL_M1= func_F_COM(q,m_list,l_link_list,data_T_all(:,94:103));
F_COM_CFL_M2= func_F_COM(q,m_list,l_link_list,data_T_all(:,104:113));
F_COM_CFL=F_COM_CFL_M1+F_COM_CFL_M2;
%重力項
F_COM_M2= func_F_COM(q,m_list,l_link_list,data_T_gravity_all(:,1:10));
F_COM_M3= func_F_COM(q,m_list,l_link_list,data_T_gravity_all(:,11:20));
F_COM_frame= func_F_COM(q,m_list,l_link_list,data_T_gravity_all(:,21:30));
F_COM_fem= func_F_COM(q,m_list,l_link_list,data_T_gravity_all(:,31:40));
F_COM_tib= func_F_COM(q,m_list,l_link_list,data_T_gravity_all(:,41:50));
F_COM_met= func_F_COM(q,m_list,l_link_list,data_T_gravity_all(:,51:60));

%外力
F_COM_hipup= func_F_COM(q,m_list,l_link_list,data_Q_hip_up(:,1:10));
F_COM_hippull= func_F_COM(q,m_list,l_link_list,data_Q_hip_pull(:,1:10));
F_COM_hipheel= func_F_COM(q,m_list,l_link_list,data_Q_heel(:,1:10));
F_COM_hiptoe= func_F_COM(q,m_list,l_link_list,data_Q_toe(:,1:10));
%股関節の支持力を除いた総和を計算
F_COM_muscletendon=F_COM_Ci_4th+F_COM_Ci_M2+F_COM_CFLT_M2+F_COM_CFLT_M3+F_COM_GEo_M3+F_COM_GEo_GEorigin+F_COM_GE++F_COM_frame_fix+F_COM_ex+F_COM_CFL_M1+F_COM_CFL_M2;
F_COM_gravity=F_COM_M2+F_COM_M3+F_COM_frame+F_hip_fem+F_COM_tib+F_COM_met;
F_COM_externalforce=F_COM_hippull+F_COM_hipheel+F_COM_hiptoe;
F_COM_all=F_COM_muscletendon+F_COM_gravity+F_COM_externalforce+F_COM_hipup;

%一般化力から股関節にかかる力を計算

if graph_view == true

    %Figure1：股関節のy方向にかかる力をそれぞれ表示1
    figure(1)
    plot(t(:,1),F_hip_Ci(:,2),'LineWidth',2);
    hold on
    plot(t(:,1),F_hip_CFLT(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_GEo(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_GE(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_CFL(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_muscletendon(:,2),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    ylim([-60 60]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$Ci$','$CFLT$','$GEo$','$GE$','$CFL$','$muscle$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forcehipy[3].pdf'], 'ContentType', 'vector');
    end
    %Figure2：股関節のy方向にかかる力をそれぞれ表示2
    figure(2)
    plot(t(:,1),F_hip_hipup(:,2),'LineWidth',2);
    hold on
    plot(t(:,1),F_hip_hippull(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_hipheel(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_hiptoe(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_frame_fix(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_ex(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_all(:,2),'LineWidth',2);
    
    hold off
    xlim([0 time_lim]);
    ylim([-50 50]);
    legend('$up$','$pull$','$heel$','$toe$','$frame_{fix}$','$ex$','$all$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forcehipy[2].pdf'], 'ContentType', 'vector');
    end
    % 
    % %Figure3：COMのy方向にかかる力をそれぞれ表示1
    % figure(3)
    % plot(t(:,1),F_COM_Ci(:,2),'LineWidth',2);
    % hold on
    % plot(t(:,1),F_COM_CFLT(:,2),'LineWidth',2);
    % plot(t(:,1),F_COM_GEo(:,2),'LineWidth',2);
    % plot(t(:,1),F_COM_GE(:,2),'LineWidth',2);
    % plot(t(:,1),F_COM_CFL(:,2),'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % %ylim([-6 6]);
    % legend('$Ci$','$CFLT$','$GEo$','$GE$','$CFL$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Time [s]');
    % ylabel('Force [N]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_forceCOMy[3].pdf'], 'ContentType', 'vector');
    % end
    % %Figure4：COMのy方向にかかる力をそれぞれ表示2
    % figure(4)
    % plot(t(:,1),F_COM_hipup(:,2),'LineWidth',2);
    % hold on
    % plot(t(:,1),F_COM_hippull(:,2),'LineWidth',2);
    % plot(t(:,1),F_COM_hipheel(:,2),'LineWidth',2);
    % plot(t(:,1),F_COM_hiptoe(:,2),'LineWidth',2);
    % plot(t(:,1),F_COM_frame_fix(:,2),'LineWidth',2);
    % plot(t(:,1),F_COM_ex(:,2),'LineWidth',2);
    % plot(t(:,1),F_COM_all(:,2),'LineWidth',2);
    % 
    % hold off
    % xlim([0 time_lim]);
    % %ylim([-6 6]);
    % legend('$up$','$pull$','$heel$','$toe$','$frame_{fix}$','$ex$','$all$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Time [s]');
    % ylabel('Force [N]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_forceCOMy[2].pdf'], 'ContentType', 'vector');
    % end

    %Figure5：股関節のx方向にかかる力をそれぞれ表示1
    figure(5)
    plot(t(:,1),F_hip_Ci(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),F_hip_CFLT(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_GEo(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_GE(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_CFL(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_muscletendon(:,1),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    ylim([-100 100]);
    legend('$Ci$','$CFLT$','$GEo$','$GE$','$CFL$','$muscle$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forcehipx[3].pdf'], 'ContentType', 'vector');
    end
    %Figure6：股関節のx方向にかかる力をそれぞれ表示2
    figure(6)
    plot(t(:,1),F_hip_hipup(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),F_hip_hippull(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_hipheel(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_hiptoe(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_frame_fix(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_ex(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_all(:,1),'LineWidth',2);
    
    hold off
    xlim([0 time_lim]);
    ylim([-15 15]);
    legend('$up$','$pull$','$heel$','$toe$','$frame_{fix}$','$ex$','$all$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forcehipx[2].pdf'], 'ContentType', 'vector');
    end

    % %Figure7：COMのx方向にかかる力をそれぞれ表示1
    % figure(7)
    % plot(t(:,1),F_COM_Ci(:,1),'LineWidth',2);
    % hold on
    % plot(t(:,1),F_COM_CFLT(:,1),'LineWidth',2);
    % plot(t(:,1),F_COM_GEo(:,1),'LineWidth',2);
    % plot(t(:,1),F_COM_GE(:,1),'LineWidth',2);
    % plot(t(:,1),F_COM_CFL(:,1),'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % %ylim([-6 6]);
    % legend('$Ci$','$CFLT$','$GEo$','$GE$','$CFL$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Time [s]');
    % ylabel('Force [N]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_forceCOMx[3].pdf'], 'ContentType', 'vector');
    % end
    % %Figure8：COMのx方向にかかる力をそれぞれ表示2
    % figure(8)
    % plot(t(:,1),F_COM_hipup(:,1),'LineWidth',2);
    % hold on
    % plot(t(:,1),F_COM_hippull(:,1),'LineWidth',2);
    % plot(t(:,1),F_COM_hipheel(:,1),'LineWidth',2);
    % plot(t(:,1),F_COM_hiptoe(:,1),'LineWidth',2);
    % plot(t(:,1),F_COM_frame_fix(:,1),'LineWidth',2);
    % plot(t(:,1),F_COM_ex(:,1),'LineWidth',2);
    % plot(t(:,1),F_COM_all(:,1),'LineWidth',2);
    % 
    % hold off
    % xlim([0 time_lim]);
    % %ylim([-6 6]);
    % legend('$up$','$pull$','$heel$','$toe$','$frame_{fix}$','$ex$','$all$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    % xlabel('Time [s]');
    % ylabel('Force [N]');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_forceCOMx[2].pdf'], 'ContentType', 'vector');
    % end

    %Figure9：トルクの総和をそれぞれ表示
    figure(9)
    plot(t(:,1),torque_all(:,5),'LineWidth',2);
    hold on
    plot(t(:,1),torque_all(:,6),'LineWidth',2);
    plot(t(:,1),torque_all(:,7),'LineWidth',2);
    plot(t(:,1),torque_all(:,8),'LineWidth',2);
    plot(t(:,1),torque_all(:,9),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    ylim([-1 1]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$\theta_{1}$','$\theta_{2}$','$\theta_{3}$','$\theta_{4}$','$\theta_{CFL}$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Torque [Nm]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_torqueall.pdf'], 'ContentType', 'vector');
    end

    %Figure10：リンク角の加速度をそれぞれ表示2
    figure(10)
    plot(t(:,1),T_all(:,5),'LineWidth',2);
    hold on
    plot(t(:,1),T_all(:,6),'LineWidth',2);
    plot(t(:,1),T_all(:,7),'LineWidth',2);
    plot(t(:,1),T_all(:,8),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-1 1]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$\theta_{1}$','$\theta_{2}$','$\theta_{3}$','$\theta_{4}$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Accel [rad/s^2]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_accelall[1].pdf'], 'ContentType', 'vector');
    end

    %Figure11：リンク角の加速度をそれぞれ表示2
    figure(11)
    plot(t(:,1),data_accel(:,5),'LineWidth',2);
    hold on
    plot(t(:,1),data_accel(:,6),'LineWidth',2);
    plot(t(:,1),data_accel(:,7),'LineWidth',2);
    plot(t(:,1),data_accel(:,8),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    ylim([-5 5]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$\theta_{1}$','$\theta_{2}$','$\theta_{3}$','$\theta_{4}$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Accel [rad/s^2]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_accelall[2].pdf'], 'ContentType', 'vector');
    end
    
    %Figure12：リンク角の加速度をそれぞれ表示2
    figure(12)
    plot(t(:,1),data_T_all(:,93),'LineWidth',2);

    xlim([0 time_lim]);
    %ylim([-5 5]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$\theta_{ankle,ex}$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Torque [Nm]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_torque.pdf'], 'ContentType', 'vector');
    end

end
