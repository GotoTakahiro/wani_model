%data_T_allの形式(:,113)
%一般化力各力，トルクがリンク角theta1~4の加速度に対して与える寄与度をプロット
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
%torque_all = zeros(size(q,1),10);
COM = zeros(size(q,1),2);
Accel_all = zeros(size(q,1),10);
Accel_Ci = zeros(size(q,1),10);
Accel_CFLT = zeros(size(q,1),10);
Accel_GEo = zeros(size(q,1),10);
Accel_GE = zeros(size(q,1),10);
Accel_frame_fix = zeros(size(q,1),10);
Accel_ex = zeros(size(q,1),10);
Accel_M2 = zeros(size(q,1),10);
Accel_M3 = zeros(size(q,1),10);
Accel_frame = zeros(size(q,1),10);
Accel_fem = zeros(size(q,1),10);
Accel_tib = zeros(size(q,1),10);
Accel_met = zeros(size(q,1),10);
Accel_up = zeros(size(q,1),10);
Accel_pull = zeros(size(q,1),10);
Accel_heel = zeros(size(q,1),10);
Accel_toe = zeros(size(q,1),10);


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


%全ての力を一般化力として扱った時の総トルクを計算
torque_all=(data_T_all(:,1:10)+data_T_all(:,11:20)+data_T_all(:,21:30)+data_T_all(:,31:40)+data_T_all(:,41:50)+data_T_all(:,51:60)+data_T_all(:,71:80));
torque_all=torque_all+data_T_gravity_all(:,1:10)+data_T_gravity_all(:,11:20)+data_T_gravity_all(:,21:30)+data_T_gravity_all(:,31:40)+data_T_gravity_all(:,41:50)+data_T_gravity_all(:,51:60);
torque_all(:,5)=torque_all(:,5)+data_T_all(:,91);
torque_all(:,8)=torque_all(:,8)+data_T_all(:,93)+data_T_all(:,92);
torque_all=torque_all+data_Q_hip_pull(:,1:10)+data_Q_heel(:,1:10)+data_Q_toe(:,1:10)+data_Q_hip_up(:,1:10);

%角度拘束に一般化力
Q_frame=[zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),data_T_all(:,91),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1)];
Q_ex=[zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),zeros(size(q,1),1),data_T_all(:,93),zeros(size(q,1),1),zeros(size(q,1),1)];

%慣性行列の逆行列を掛け合わせたものを計算
for i = 1:size(q,1)
    M=Inertial_matrix(m_list,l_link_list,transpose(q(i,1:10)));
    Accel_all(i,1:10) =torque_all(i,1:10)/M;
    Accel_Ci(i,1:10)=(data_T_all(i,1:10)+data_T_all(i,11:20))/M;
    Accel_CFLT(i,1:10)=(data_T_all(i,21:30)+data_T_all(i,31:40))/M;
    Accel_GEo(i,1:10)=(data_T_all(i,41:50)+data_T_all(i,51:60))/M;
    Accel_GE(i,1:10)=data_T_all(i,71:80)/M;
    Accel_frame_fix(i,1:10)=Q_frame(i,1:10)/M;
    Accel_ex(i,1:10)=Q_ex(i,1:10)/M;
    Accel_M2(i,1:10)=data_T_gravity_all(i,1:10)/M;
    Accel_M3(i,1:10)=data_T_gravity_all(i,11:20)/M;
    Accel_frame(i,1:10)=data_T_gravity_all(i,21:30)/M;
    Accel_fem(i,1:10)=data_T_gravity_all(i,31:40)/M;
    Accel_tib(i,1:10)=data_T_gravity_all(i,41:50)/M;
    Accel_met(i,1:10)=data_T_gravity_all(i,51:60)/M;
    Accel_up(i,1:10)=data_Q_hip_up(i,1:10)/M;
    Accel_pull(i,1:10)=data_Q_hip_pull(i,1:10)/M;
    Accel_heel(i,1:10)=data_Q_heel(i,1:10)/M;
    Accel_toe(i,1:10)=data_Q_toe(i,1:10)/M;


end


%一般化力から股関節にかかる力を計算

if graph_view == true
    %Figure1：theta1のリンク加速度の成分を表示1
    figure(1)
    plot(t(:,1),Accel_Ci(:,5),'LineWidth',2);
    hold on
    plot(t(:,1),Accel_CFLT(:,5),'LineWidth',2);
    plot(t(:,1),Accel_GEo(:,5),'LineWidth',2);
    plot(t(:,1),Accel_GE(:,5),'LineWidth',2);
    plot(t(:,1),Accel_all(:,5),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-5 5]);
    legend('$Ci$','$CFLT$','$GEo$','$GE$','All','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Accel [rad/s^2]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_acceltheta1[1].pdf'], 'ContentType', 'vector');
    end

    %Figure2：theta1のリンク加速度の成分を表示2
    figure(2)
    plot(t(:,1),Accel_up(:,5),'LineWidth',2);
    hold on
    plot(t(:,1),Accel_pull(:,5),'LineWidth',2);
    plot(t(:,1),Accel_heel(:,5),'LineWidth',2);
    plot(t(:,1),Accel_toe(:,5),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-5 5]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('up','pull','heel','toe','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Accel [rad/s^2]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_acceltheta1[2].pdf'], 'ContentType', 'vector');
    end
    
    %Figure3：theta2のリンク加速度の成分を表示1
    figure(3)
    plot(t(:,1),Accel_Ci(:,6),'LineWidth',2);
    hold on
    plot(t(:,1),Accel_CFLT(:,6),'LineWidth',2);
    plot(t(:,1),Accel_GEo(:,6),'LineWidth',2);
    plot(t(:,1),Accel_GE(:,6),'LineWidth',2);
    plot(t(:,1),Accel_all(:,6),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-5 5]);
    legend('$Ci$','$CFLT$','$GEo$','$GE$','All','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Accel [rad/s^2]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_acceltheta2[1].pdf'], 'ContentType', 'vector');
    end

    %Figure4：theta2のリンク加速度の成分を表示2
    figure(4)
    plot(t(:,1),Accel_up(:,6),'LineWidth',2);
    hold on
    plot(t(:,1),Accel_pull(:,6),'LineWidth',2);
    plot(t(:,1),Accel_heel(:,6),'LineWidth',2);
    plot(t(:,1),Accel_toe(:,6),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-5 5]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('up','pull','heel','toe','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Accel [rad/s^2]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_acceltheta2[2].pdf'], 'ContentType', 'vector');
    end
    %Figure5：theta3のリンク加速度の成分を表示1
    figure(5)
    plot(t(:,1),Accel_Ci(:,7),'LineWidth',2);
    hold on
    plot(t(:,1),Accel_CFLT(:,7),'LineWidth',2);
    plot(t(:,1),Accel_GEo(:,7),'LineWidth',2);
    plot(t(:,1),Accel_GE(:,7),'LineWidth',2);
    plot(t(:,1),Accel_all(:,7),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-5 5]);
    legend('$Ci$','$CFLT$','$GEo$','$GE$','All','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Accel [rad/s^2]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_acceltheta3[1].pdf'], 'ContentType', 'vector');
    end

    %Figure6：theta3のリンク加速度の成分を表示2
    figure(6)
    plot(t(:,1),Accel_up(:,7),'LineWidth',2);
    hold on
    plot(t(:,1),Accel_pull(:,7),'LineWidth',2);
    plot(t(:,1),Accel_heel(:,7),'LineWidth',2);
    plot(t(:,1),Accel_toe(:,7),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-5 5]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('up','pull','heel','toe','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Accel [rad/s^2]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_acceltheta3[2].pdf'], 'ContentType', 'vector');
    end
    %Figure7：theta4のリンク加速度の成分を表示1
    figure(7)
    plot(t(:,1),Accel_Ci(:,8),'LineWidth',2);
    hold on
    plot(t(:,1),Accel_CFLT(:,8),'LineWidth',2);
    plot(t(:,1),Accel_GEo(:,8),'LineWidth',2);
    plot(t(:,1),Accel_GE(:,8),'LineWidth',2);
    plot(t(:,1),Accel_all(:,8),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-5 5]);
    legend('$Ci$','$CFLT$','$GEo$','$GE$','All','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Accel [rad/s^2]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_acceltheta4[1].pdf'], 'ContentType', 'vector');
    end

    %Figure8：theta4のリンク加速度の成分を表示2
    figure(8)
    plot(t(:,1),Accel_up(:,8),'LineWidth',2);
    hold on
    plot(t(:,1),Accel_pull(:,8),'LineWidth',2);
    plot(t(:,1),Accel_heel(:,8),'LineWidth',2);
    plot(t(:,1),Accel_toe(:,8),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-5 5]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('up','pull','heel','toe','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Accel [rad/s^2]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_acceltheta4[2].pdf'], 'ContentType', 'vector');
    end
    
    %Figure10：リンク角の加速度をそれぞれ表示2
    figure(10)
    plot(t(:,1),Accel_all(:,5)-Accel_all(:,6),'LineWidth',2);
    hold on
    plot(t(:,1),Accel_all(:,6)-Accel_all(:,7),'LineWidth',2);
    plot(t(:,1),Accel_all(:,7)-Accel_all(:,8),'LineWidth',2);
    plot(t(:,1),Accel_all(:,8),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    %ylim([-1 1]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$\phi_{1}$','$\phi_{2}$','$\phi_{3}$','$\phi_{4}$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Accel [rad/s^2]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_accelall[2].pdf'], 'ContentType', 'vector');
    end
    

end
