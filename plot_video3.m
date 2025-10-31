% clear;
 close all;
% clearvars
load('/Users/goto/Documents/Matlab_goto/crocodile_sim_PID-main/results/plotall3/exp20251028_CFL350_Ci44_CFLT100_GEo35_GE185.mat');
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

% 筋腱の張力，トルク，運動量，角運動量の時系列データを計算
for i = 1:size(q,1)
    k_wire = data_k_c_wire(i,2:5); %k_Ci, k_CFLT, k_GEo, k_GE
    c_wire = data_k_c_wire(i,6:10); %c_Ci, c_CFLT, c_GEo, c_GE
    
    general_q = q(i,1:10).';
    general_dq = q(i,11:20).';

    % 張力，トルク，運動量，角運動量の計算．運動量と角運動量はcal_EOM_work.mで計算した関数を用いている
    muscle_tension(i,:) = (calc_muscle_tension(l_link_list,l_muscle_list,k_wire.',c_wire.',general_q, general_dq))';
    torque_muscle(i,:) = (calc_torque_muscle(l_link_list,l_muscle_list,k_wire.',c_wire.',general_q, general_dq))';
    momentum_list(i,:) = calc_momentum_list(m_list,l_link_list,q(i,1:10).',q(i,11:20).')';
    angular_moment_list(i,:) = calc_angular_moment_list(m_list,l_link_list,q(i,1:10).',q(i,11:20).')';
    momentum_COM_list(i,:) = calc_momentum_COM_list(m_list,l_link_list,q(i,1:10).',q(i,11:20).')';

    % -data_Q(:,11)が5以上になったところから最後までを取得
    if tension_flag == false && -data_Q(i,11) > 5 && i > 200
        tension_start_index = i;
        tension_flag = true;
    end
end

% 骨格と筋腱の間の角度を計算するための座標を計算．（結局あまり使わなかったので無視してもOK）
coordinates_x_for_angle = coordinates_x(tension_start_index:end,:);
coordinates_y_for_angle = coordinates_y(tension_start_index:end,:);

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
%theta1
torque_all=(data_T_all(:,1:10)+data_T_all(:,11:20)+data_T_all(:,21:30)+data_T_all(:,31:40)+data_T_all(:,41:50)+data_T_all(:,51:60)+data_T_all(:,71:80));%+data_T_all(:,81:90));
torque_all=torque_all+data_T_gravity_all(:,1:10)+data_T_gravity_all(:,11:20)+data_T_gravity_all(:,21:30)+data_T_gravity_all(:,31:40)+data_T_gravity_all(:,41:50)+data_T_gravity_all(:,51:60);
torque_all(:,5)=torque_all(:,5)+data_T_all(:,91);
torque_all(:,8)=torque_all(:,8)+data_T_all(:,93);
torque_all=torque_all+data_Q_hip_pull(:,1:10)+data_Q_heel(:,1:10)+data_Q_toe(:,1:10)+data_Q_hip_up(:,1:10);

%一般化力からhipの質点に働く力を逆算するための関数
%ワイヤー項
F_hip_Ci_4th= func_F_hip(q,l_link_list,data_T_all(:,1:10));
F_hip_Ci_M2= func_F_hip(q,l_link_list,data_T_all(:,11:20));
F_hip_CFLT_M2= func_F_hip(q,l_link_list,data_T_all(:,21:30));
F_hip_CFLT_M3= func_F_hip(q,l_link_list,data_T_all(:,31:40)); %0
F_hip_GEo_M3= func_F_hip(q,l_link_list,data_T_all(:,41:50));
F_hip_GEo_GEorigin= func_F_hip(q,l_link_list,data_T_all(:,51:60));
%F_hip_GE_M3= func_F_hip(q,l_link_list,data_T_all(:,61:70));
F_hip_GE= func_F_hip(q,l_link_list,data_T_all(:,71:80));
F_hip_GE_pulley= func_F_hip(q,l_link_list,data_T_all(:,81:90));
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
F_hip_all=F_hip_Ci_4th+F_hip_Ci_M2+F_hip_CFLT_M2+F_hip_CFLT_M3+F_hip_GEo_M3+F_hip_GEo_GEorigin+F_hip_GE;
F_hip_all=F_hip_all+F_hip_M2+F_hip_M3+F_hip_frame+F_hip_fem+F_hip_tib+F_hip_met;
F_hip_all=F_hip_all+F_hip_hippull+F_hip_hipheel+F_hip_hiptoe;

%一般化力からCOMの質点に働く力を逆算するための関数
%ワイヤー項
F_COM_Ci_4th= func_F_COM(q,m_list,l_link_list,data_T_all(:,1:10));
F_COM_Ci_M2= func_F_COM(q,m_list,l_link_list,data_T_all(:,11:20));
F_COM_CFLT_M2= func_F_COM(q,m_list,l_link_list,data_T_all(:,21:30));
F_COM_CFLT_M3= func_F_COM(q,m_list,l_link_list,data_T_all(:,31:40)); %0
F_COM_GEo_M3= func_F_COM(q,m_list,l_link_list,data_T_all(:,41:50));
F_COM_GEo_GEorigin= func_F_COM(q,m_list,l_link_list,data_T_all(:,51:60));
%F_COM_GE_M3= func_F_COM(q,m_list,l_link_list,data_T_all(:,61:70));
F_COM_GE= func_F_COM(q,m_list,l_link_list,data_T_all(:,71:80));
F_COM_GE_pulley= func_F_COM(q,m_list,l_link_list,data_T_all(:,81:90));
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
F_COM_all=F_COM_Ci_4th+F_COM_Ci_M2+F_COM_CFLT_M2+F_COM_CFLT_M3+F_COM_GEo_M3+F_COM_GEo_GEorigin+F_COM_GE;
F_COM_all=F_COM_all+F_COM_M2+F_COM_M3+F_COM_frame+F_hip_fem+F_COM_tib+F_COM_met;
F_COM_all=F_COM_all+F_COM_hippull+F_COM_hipheel+F_COM_hiptoe;

%一般化力から股関節にかかる力を計算

if graph_view == true

    %Figure1：質点M1にかかる力をそれぞれ表示
    figure(1)
    hold on
    plot(t(:,1),data_T_all(:,2),'LineWidth',2);
        plot(t(:,1),data_T_all(:,12),'LineWidth',2); %0
        plot(t(:,1),data_T_all(:,22),'LineWidth',2); %0
    %    plot(t(:,1),data_T_all(:,32),'LineWidth',2); %0
    %    plot(t(:,1),data_T_all(:,42),'LineWidth',2); %0
    plot(t(:,1),data_T_all(:,52),'LineWidth',2);
    %    plot(t(:,1),data_T_all(:,62),'LineWidth',2); %0
    plot(t(:,1),data_T_all(:,72),'LineWidth',2);
    %plot(t(:,1),data_T_all(:,82),'LineWidth',2);
    %     plot(t(:,1),data_T_gravity_all(:,2),'LineWidth',2); %0
    %     plot(t(:,1),data_T_gravity_all(:,12),'LineWidth',2); %0
    % plot(t(:,1),data_T_gravity_all(:,22),'LineWidth',2);
    % plot(t(:,1),data_T_gravity_all(:,32),'LineWidth',2);
    % plot(t(:,1),data_T_gravity_all(:,42),'LineWidth',2);
    % plot(t(:,1),data_T_gravity_all(:,52),'LineWidth',2);
    %     plot(t(:,1),data_Q_hip_pull(:,2),'LineWidth',2); %ほぼ0作用線的に
    % plot(t(:,1),data_Q_hip_up(:,2),'LineWidth',2);
    % plot(t(:,1),data_Q_heel(:,2),'LineWidth',2);
    % plot(t(:,1),data_Q_toe(:,2),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-15 55]);
    legend('$Ci_{4th}$','$Ci_{M2}$','$CFLT_{M2}$','$GEo_{GEori}$','$GE_{wire}$','$GE_{ankle}$', ...
            '$g_{M2}$','$g_{M3}$','$g_{frame}$','$g_{fem}$','$g_{tib}$','$g_{met}$','$F_{pull}$','$F_{up}$','$F_{heel}$','$F_{toe}$',...
            'Interpreter', 'latex', 'FontSize',12,'Location','northeast', 'NumColumns',5);
    % legend('$Ci_{4th}$','$Ci_{M2}$','$CFLT_{M2}$','$GEo_{GEori}$','$GE_{M3,p}$','$GE_{t,p}$', ...
    %         '$g_{M2}$','$g_{M3}$','$g_{frame}$','$g_{fem}$','$g_{tib}$','$g_{met}$','$F_{pull}$','$F_{up}$','$F_{heel}$','$F_{toe}$',...
    %         'Interpreter', 'latex', 'FontSize',12,'Location','northeast', 'NumColumns',5);
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_force_M2y[1].pdf'], 'ContentType', 'vector');
    end
    %Figure2：質点M1にかかる力をそれぞれ表示
    figure(2)
    hold on
    % plot(t(:,1),-data_T_all(:,2),'LineWidth',2);
    %     plot(t(:,1),data_T_all(:,12),'LineWidth',2); %0
    %     plot(t(:,1),data_T_all(:,22),'LineWidth',2); %0
    %     plot(t(:,1),data_T_all(:,32),'LineWidth',2); %0
    %     plot(t(:,1),data_T_all(:,42),'LineWidth',2); %0
    % plot(t(:,1),-data_T_all(:,52),'LineWidth',2);
    %     plot(t(:,1),data_T_all(:,62),'LineWidth',2); %0
    % plot(t(:,1),-data_T_all(:,82),'LineWidth',2);
    % plot(t(:,1),-data_T_all(:,72),'LineWidth',2);
    %   plot(t(:,1),data_T_gravity_all(:,2),'LineWidth',2); %0
    %    plot(t(:,1),data_T_gravity_all(:,12),'LineWidth',2); %0
    %plot(t(:,1),data_T_gravity_all(:,22),'LineWidth',2);
    %plot(t(:,1),data_T_gravity_all(:,32),'LineWidth',2);
    %plot(t(:,1),data_T_gravity_all(:,42),'LineWidth',2);
    %plot(t(:,1),data_T_gravity_all(:,52),'LineWidth',2);
    %    plot(t(:,1),data_Q_hip_pull(:,2),'LineWidth',2); %ほぼ0作用線的に
    plot(t(:,1),data_Q_hip_up(:,2),'LineWidth',2);
    plot(t(:,1),data_Q_heel(:,2),'LineWidth',2);
    plot(t(:,1),data_Q_toe(:,2),'LineWidth',2);
    plot(t(:,1),data_Q_toe(:,2),'LineWidth',2);
    plot(t(:,1),torque_all(:,2),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-15 55]);

    legend( ...
            '$F_{up}$','$F_{heel}$','$F_{toe}$','$F_{all}$',...
            'Interpreter', 'latex', 'FontSize',12,'Location','northeast', 'NumColumns',5);
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_force_M2y[2].pdf'], 'ContentType', 'vector');
    end
    %Figure3：股関節のy方向にかかる力の総和をそれぞれ表示
    figure(3)
    plot(t(:,1),F_hip_Ci_4th(:,2),'LineWidth',2);
    hold on
    plot(t(:,1),F_hip_Ci_M2(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_CFLT_M2(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_CFLT_M3(:,2),'LineWidth',2);
    %plot(t(:,1),F_hip_GEo_M3(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_GEo_GEorigin(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_GE(:,2),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-6 6]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$Ci_{4th}$','$Ci_{M2}$','$CFLT_{M2}$','$CFLT_{M3}$','$GEo_{origin}$','$GE$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forcehipy[1].pdf'], 'ContentType', 'vector');
    end
    %Figure4：股関節のy方向にかかる力の総和をそれぞれ表示
    figure(4)
    plot(t(:,1),F_hip_hipup(:,2),'LineWidth',2);
    hold on
    plot(t(:,1),F_hip_hippull(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_hipheel(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_hiptoe(:,2),'LineWidth',2);
    plot(t(:,1),F_hip_all(:,2),'LineWidth',2);
    
    hold off
    xlim([0 time_lim]);
    %ylim([-6 6]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$up$','$pull$','$heel$','$toe$','$all$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forcehipy[2].pdf'], 'ContentType', 'vector');
    end

    %Figure5：COMのy方向にかかる力の総和をそれぞれ表示
    figure(5)
    plot(t(:,1),F_COM_Ci_4th(:,2),'LineWidth',2);
    hold on
    plot(t(:,1),F_COM_Ci_M2(:,2),'LineWidth',2);
    plot(t(:,1),F_COM_CFLT_M2(:,2),'LineWidth',2);
    plot(t(:,1),F_COM_CFLT_M3(:,2),'LineWidth',2);
    plot(t(:,1),F_COM_GEo_M3(:,2),'LineWidth',2);
    plot(t(:,1),F_COM_GEo_GEorigin(:,2),'LineWidth',2);
    plot(t(:,1),F_COM_GE(:,2),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-6 6]);
    legend('$Ci_{4th}$','$Ci_{M2}$','$CFLT_{M2}$','$CFLT_{M3}$','$GEo_{M3}$','$GEo_{origin}$','$GE$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forceCOMy[1].pdf'], 'ContentType', 'vector');
    end
    %Figure6：COMのy方向にかかる力の総和をそれぞれ表示
    figure(6)
    plot(t(:,1),F_COM_hipup(:,2),'LineWidth',2);
    hold on
    plot(t(:,1),F_COM_hippull(:,2),'LineWidth',2);
    plot(t(:,1),F_COM_hipheel(:,2),'LineWidth',2);
    plot(t(:,1),F_COM_hiptoe(:,2),'LineWidth',2);
    plot(t(:,1),F_COM_all(:,2),'LineWidth',2);
    
    hold off
    xlim([0 time_lim]);
    %ylim([-6 6]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$up$','$pull$','$heel$','$toe$','$all$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forceCOMy[2].pdf'], 'ContentType', 'vector');
    end

    %Figure7：股関節のx方向にかかる力の総和をそれぞれ表示
    figure(7)
    plot(t(:,1),F_hip_Ci_4th(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),F_hip_Ci_M2(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_CFLT_M2(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_CFLT_M3(:,1),'LineWidth',2);
    %plot(t(:,1),F_hip_GEo_M3(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_GEo_GEorigin(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_GE(:,1),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-6 6]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$Ci_{4th}$','$Ci_{M2}$','$CFLT_{M2}$','$CFLT_{M3}$','$GEo_{origin}$','$GE$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forcehipx[1].pdf'], 'ContentType', 'vector');
    end
    %Figure8：股関節のx方向にかかる力の総和をそれぞれ表示
    figure(8)
    plot(t(:,1),F_hip_hipup(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),F_hip_hippull(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_hipheel(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_hiptoe(:,1),'LineWidth',2);
    plot(t(:,1),F_hip_all(:,1),'LineWidth',2);
    
    hold off
    xlim([0 time_lim]);
    %ylim([-6 6]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$up$','$pull$','$heel$','$toe$','$all$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forcehipx[2].pdf'], 'ContentType', 'vector');
    end

    %Figure9：COMのx方向にかかる力の総和をそれぞれ表示
    figure(9)
    plot(t(:,1),F_COM_Ci_4th(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),F_COM_Ci_M2(:,1),'LineWidth',2);
    plot(t(:,1),F_COM_CFLT_M2(:,1),'LineWidth',2);
    plot(t(:,1),F_COM_CFLT_M3(:,1),'LineWidth',2);
    plot(t(:,1),F_COM_GEo_M3(:,1),'LineWidth',2);
    plot(t(:,1),F_COM_GEo_GEorigin(:,1),'LineWidth',2);
    plot(t(:,1),F_COM_GE(:,1),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-6 6]);
    legend('$Ci_{4th}$','$Ci_{M2}$','$CFLT_{M2}$','$CFLT_{M3}$','$GEo_{M3}$','$GEo_{origin}$','$GE$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forceCOMy[1].pdf'], 'ContentType', 'vector');
    end
    %Figure10：COMのx方向にかかる力の総和をそれぞれ表示
    figure(10)
    plot(t(:,1),F_COM_hipup(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),F_COM_hippull(:,1),'LineWidth',2);
    plot(t(:,1),F_COM_hipheel(:,1),'LineWidth',2);
    plot(t(:,1),F_COM_hiptoe(:,1),'LineWidth',2);
    plot(t(:,1),F_COM_all(:,1),'LineWidth',2);
    
    hold off
    xlim([0 time_lim]);
    %ylim([-6 6]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$up$','$pull$','$heel$','$toe$','$all$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_forceCOMx[2].pdf'], 'ContentType', 'vector');
    end

end
