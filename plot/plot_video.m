% clear;
 close all;
% clearvars

addpath('/Users/goto/Documents/Matlab/crocodile_sim_PID-main');
load('/Users/goto/Documents/MATLAB/crocodile_sim_PID-main/plot/results/m4/exp20251028_CFL350_Ci44_CFLT100_GEo35_GE185.mat');
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
movie_save = true;

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
    k_wire = data_k_c_wire(i,2:8); %k_Ci, k_CFLT, k_GEo, k_GE
    c_wire = data_k_c_wire(i,9:16); %c_Ci, c_CFLT, c_GEo, c_GE
    % L_wire(1) = data_acce_GRF(i,18);

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

if graph_view == true

    %Figure1：筋腱の張力
    figure(1)
    plot(t(:,1),-data_Q(:,11),'-','LineWidth',2,'Color',CFL_Color); %CFLの張力
    hold on
    plot(t(:,1),muscle_tension(:,1),'--','LineWidth',2,'Color', Ci_Color); %Ciの張力
    if noCFLT == false
        plot(t(:,1),muscle_tension(:,2),'-.','LineWidth',2,'Color', CFLT_Color); %CFLTの張力
    end
    plot(t(:,1),muscle_tension(:,3),':','LineWidth',2,'Color',GEo_Color); %GEoの張力
    if noGE == false
        plot(t(:,1),muscle_tension(:,4),'LineWidth',2,'Color',GE_Color); %GEの張力
    end
    hold off
    xlim([0 time_lim]);
    %ylim([-10 500]);
    h_axes = gca;
    h_axes.XAxis.FontSize = 20;
    h_axes.YAxis.FontSize = 20;
    % ylim([-1 1]);
    xlabel('Time [s]','FontSize',25);
    ylabel('Tension [N]','FontSize',25);
    if noCFLT == false && noGE == false
        legend('CFL','Ci','CFLT','GEo','GE','FontSize',20,'Location','best');
    elseif noCFLT == false && noGE == true
        legend('CFL','Ci','CFLT','GEo','FontSize',20,'Location','best');
    elseif noCFLT == true && noGE == false
        legend('CFL','Ci','GEo','GE','FontSize',20,'Location','best');
    end
    if graph_save == true
        exportgraphics(gca, [save_path name '_tension.pdf'], 'ContentType', 'vector');
    end


    % Figure2：筋腱によるトルク
    figure(2)
    plot(t(:,1),torque_muscle(:,6),'LineWidth',2); %hipのトルク
    hold on
    plot(t(:,1),torque_muscle(:,7),'LineWidth',2); %kneeのトルク
    plot(t(:,1),torque_muscle(:,8),'LineWidth',2); %ankleのトルク
    hold 
    xlim([0 time_lim]);
    ylim([-25 2]);
    h_axes = gca;
    h_axes.XAxis.FontSize = 20;
    h_axes.YAxis.FontSize = 20;
    xlabel('Time [s]','FontSize',25);
    ylabel('Torque [Nm]','FontSize',25);
    legend('hip','knee','ankle','FontSize',25, 'Location','best');
    if graph_save == true
        exportgraphics(gca, [save_path name '_torque.pdf'], 'ContentType', 'vector');
    end


    %Figure3：関節角度
    figure(3)
    plot(t(:,1),rad2deg(q(:,6)),'LineWidth',2); %hipの角度
    hold on
    plot(t(:,1),rad2deg(q(:,7)),'LineWidth',2); %kneeの角度
    plot(t(:,1),rad2deg(q(:,8)),'LineWidth',2); %ankleの角度
    hold off
    xlim([0 time_lim]);
    ylim([-200 150]);
    h_axes = gca;
    h_axes.XAxis.FontSize = 20;
    h_axes.YAxis.FontSize = 20;
    legend('hip','knee','ankle','FontSize',20,'Location','northeast');
    % 判例を取得
    lgd = legend;
    % 判例の現在の位置を取得
    currentPosition = lgd.Position;
    % Y位置を少し下げる
    newPosition = currentPosition;
    newPosition(2) = newPosition(2) - 0.2;  % Y座標を減らすことで下に移動
    % 新しい位置を設定
    lgd.Position = newPosition;
    xlabel('Time [s]','FontSize',25);
    ylabel('Angle [deg]','FontSize',25);
    if graph_save == true
        exportgraphics(gca, [save_path name '_angle.pdf'], 'ContentType', 'vector');
    end


    %Figure4：地面反力
    GRF = data_accel_GRF(:,14:17);
    figure(4)
    plot(t(:,1),GRF(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),GRF(:,2),'LineWidth',2);
    plot(t(:,1),GRF(:,3),'LineWidth',2);
    plot(t(:,1),GRF(:,4),'LineWidth',2);
    plot(t(:,1),(GRF(:,1)+GRF(:,3)),'LineWidth',2);
    plot(t(:,1),(GRF(:,2)+GRF(:,4)),'LineWidth',2);
    % GRF = 10のラインで黒色で点線を引く
    plot([0 time_lim],[10 10],'--','Color','k','LineWidth',1);
    hold off
    xlim([0 time_lim]);
    ylim([-15 55]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$','total-$F_x$','total-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_GRF.pdf'], 'ContentType', 'vector');
    end

    % % Figure5：仕事率
    % figure(5)
    % plot(t(50:end,1),power,'LineWidth',2);
    % xlabel('Time [s]');
    % ylabel('Power [W]');
    % xlim([0 time_lim]);
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_power.pdf'], 'ContentType', 'vector');
    % end

    % Figure6：CFLの目標値との誤差
    % figure(6)
    % plot(t(:,1), data_error_CFL(:,1),'LineWidth',2);
    % xlabel('Time [s]');
    % ylabel('Error');
    % xlim([0 time_lim]);
    % ylim([-0.01 0.01]);
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_error.pdf'], 'ContentType', 'vector');
    % end

    % % Figure7：GEoとfemurの間の角度
    % angle_femur_GEo = calc_angle_femur_GEo(coordinates_x_for_angle, coordinates_y_for_angle);
    % figure(7)
    % yyaxis left
    % plot(t(tension_start_index:end,1)-t(tension_start_index,1),rad2deg(q(tension_start_index:end,6)),'LineWidth',2);
    % ylabel('Hip angle [deg]');
    % ylim([-200 150]);
    % yyaxis right
    % plot(t(tension_start_index:end,1)-t(tension_start_index,1),angle_femur_GEo,'LineWidth',2);
    % xlabel('Time [s]');
    % ylabel('GEo angle [deg]');
    % xlim([0 time_lim-t(tension_start_index,1)]);
    % ylim([0 180]);
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_angle_femur_GEo.pdf'], 'ContentType', 'vector');
    % end

    % Figure8：横軸を時間，縦軸を関節角度と股関節トルク
    figure(8)
    yyaxis left
    plot(t(:,1),rad2deg(q(:,6)),'LineWidth',2);
    ylabel('Angle [deg]');
    ylim([-130 0]);
    yyaxis right
    plot(t(:,1),torque_muscle(:,6),'LineWidth',2);
    ylabel('Torque [Nm]');
    ylim([-25 1]);
    xlabel('Time [s]');
    xlim([0 time_lim]);
    if graph_save == true
        exportgraphics(gca, [save_path name '_torque_angle_hip.pdf'], 'ContentType', 'vector');
    end
    
    % Figure9：横軸を時間，縦軸を関節角度と膝関節トルク
    figure(9)
    yyaxis left
    plot(t(:,1),rad2deg(q(:,7)),'LineWidth',2);
    ylabel('Angle [deg]');
    ylim([-130 0]);
    yyaxis right
    plot(t(:,1),torque_muscle(:,7),'LineWidth',2);
    ylabel('Torque [Nm]');
    ylim([-25 1]);
    xlabel('Time [s]');
    xlim([0 time_lim]);
    if graph_save == true
        exportgraphics(gca, [save_path name '_torque_angle_knee.pdf'], 'ContentType', 'vector');
    end

    % % Figure10：重心の並進の運動量
    % figure(10)
    % % plot(t(:,1),momentum_list(:,1),'LineWidth',2);
    % % hold on
    % plot(t(:,1),momentum_list(:,2),'LineWidth',2);
    % hold on
    % plot(t(:,1),momentum_list(:,3),'LineWidth',2);
    % plot(t(:,1),momentum_list(:,4),'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % ylim([-0.002 0.002]);
    % xlabel('Time [s]');
    % ylabel('Momentum [J・s]');
    % % legend('frame','femur','tibia','metatarsal');
    % legend('femur','tibia','metatarsal');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_momentum.pdf'], 'ContentType', 'vector');
    % end

    % % Figure11：原点から見た角運動量
    % figure(11)
    % % plot(t(:,1),angular_moment_list(:,1),'LineWidth',2);
    % % hold on
    % plot(t(:,1),angular_moment_list(:,2),'LineWidth',2);
    % hold on
    % plot(t(:,1),angular_moment_list(:,3),'LineWidth',2);
    % plot(t(:,1),angular_moment_list(:,4),'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % ylim([-0.0007 0.0007]);
    % xlabel('Time [s]');
    % ylabel('Angular momentum [J・s]');
    % % legend('frame','femur','tibia','metatarsal');
    % legend('femur','tibia','metatarsal');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_angular_momentum.pdf'], 'ContentType', 'vector');
    % end

    % % Figure12：重心の並進運動量
    % figure(12)
    % % plot(t(:,1),momentum_COM_list(:,1),'LineWidth',2);
    % % hold on
    % % plot(t(:,1),momentum_COM_list(:,2),'LineWidth',2);
    % plot(t(:,1),momentum_COM_list(:,3),'LineWidth',2);
    % hold on
    % plot(t(:,1),momentum_COM_list(:,4),'LineWidth',2);
    % plot(t(:,1),momentum_COM_list(:,5),'LineWidth',2);
    % plot(t(:,1),momentum_COM_list(:,6),'LineWidth',2);
    % plot(t(:,1),momentum_COM_list(:,7),'LineWidth',2);
    % plot(t(:,1),momentum_COM_list(:,8),'LineWidth',2);
    % hold off
    % xlim([0 time_lim]);
    % ylim([-0.04 0.04]);
    % xlabel('Time [s]');
    % ylabel('Momentum [kg・m/s]');
    % % legend('caudal x','caudal y','femur x','femur y','tibia x','tibia y','metatarsal x','metatarsal y');
    % legend('femur x','femur y','tibia x','tibia y','metatarsal x','metatarsal y');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_momentum_COM.pdf'], 'ContentType', 'vector');
    % end

    % Figure13：angle_femur_GEo = calc_angle_femur_GEo(coordinates_x, coordinates_y);
    % figure(13)
    % plot(angular_moment_list(tension_start_index:end,2), angular_moment_list(tension_start_index:end,3),'LineWidth',2);
    % xlabel('Femur angular momentum [J・s]');
    % xlim([-0.0007 0.0007]);
    % ylabel('Tibia angular momentum [J・s]');
    % ylim([-0.0007 0.0007]);
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_angular_momentum_femur_tibia.pdf'], 'ContentType', 'vector');
    % end
    
    % % Figure14：angle_femur_GEo = calc_angle_femur_GEo(coordinates_x_for_angle, coordinates_y_for_angle);
    % figure(14)
    % plot( angle_femur_GEo, rad2deg(q(tension_start_index:end,6)), 'LineWidth',2);
    % xlabel('GEo angle [deg]');
    % xlim([0 180]);
    % ylabel('Hip angle [deg]');
    % ylim([-130 0]);
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_angle_femur_GEo_rev.pdf'], 'ContentType', 'vector');
    % end

    % % Figure15：横軸股関節角度，縦軸筋腱の張力
    % figure(15)
    % plot(rad2deg(q(:,6)),-data_Q(:,11),'-','LineWidth',2,'Color','#0072BD') %CFL
    % hold on
    % plot(rad2deg(q(:,6)),muscle_tension(:,1),'--','LineWidth',2,'Color', '#D95319') %Ci
    % if noCFLT == false
    %     plot(rad2deg(q(:,6)),muscle_tension(:,2),'-.','LineWidth',2,'Color', '#EDB120') %CFLT
    % end
    % plot(rad2deg(q(:,6)),muscle_tension(:,3),':','LineWidth',2,'Color','#7E2F8E') %GEo
    % if noGE == false
    %     plot(rad2deg(q(:,6)),muscle_tension(:,4),'LineWidth',2,'Color','#77AC30') %GE
    % end
    % hold off
    % xlabel('Hip angle [deg]','FontSize',25);
    % ylabel('Tension [N]','FontSize',25);
    % xlim([-95 -5]);
    % ylim([-10 500]);
    % h_axes = gca;
    % h_axes.XAxis.FontSize = 20;
    % h_axes.YAxis.FontSize = 20;
    % if noCFLT == false && noGE == false
    %     legend('CFL','Ci','CFLT','GEo','GE','FontSize',20,'Location','best');
    % elseif noCFLT == false && noGE == true
    %     legend('CFL','Ci','CFLT','GEo','FontSize',20,'Location','best');
    % elseif noCFLT == true && noGE == false
    %     legend('CFL','Ci','GEo','GE','FontSize',20,'Location','best');
    % end
    % lgd = legend;
    % % 判例の現在の位置を取得
    % currentPosition = lgd.Position;
    % % Y位置を少し下げる
    % newPosition = currentPosition;
    % newPosition(1) = newPosition(1) + 0.3;  % x座標
    % newPosition(2) = newPosition(2) - 0.1;  % y座標
    % lgd.Position = newPosition;
    % set(gca, 'XDir', 'reverse');
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_tension_hip.pdf'], 'ContentType', 'vector');
    % end

    % % Figure16：横軸股関節角度，縦軸筋腱のトルク
    % figure(16)
    % plot(rad2deg(q(:,6)),torque_muscle(:,6),'LineWidth',2);
    % hold on
    % plot(rad2deg(q(:,6)),torque_muscle(:,7),'LineWidth',2);
    % plot(rad2deg(q(:,6)),torque_muscle(:,8),'LineWidth',2);
    % hold off
    % xlabel('Hip angle [deg]');
    % ylabel('Torque [Nm]');
    % xlim([-95 -5]);
    % ylim([-25 2]);
    % h_axes = gca;
    % h_axes.XAxis.FontSize = 20;
    % h_axes.YAxis.FontSize = 20;
    % legend('hip','knee','ankle','FontSize',25, 'Location','best');
    % set(gca, 'XDir', 'reverse');
    % 
    % if graph_save == true
    %     exportgraphics(gca, [save_path name '_torque_hip_max25.pdf'], 'ContentType', 'vector');
    % end

    %Figure17：外力をそれぞれ表示
    GRF = data_force(:,2:9);
    figure(17)
    plot(t(:,1),GRF(:,3),'LineWidth',2);
    hold on
    plot(t(:,1),GRF(:,4),'LineWidth',2);
    plot(t(:,1),GRF(:,5),'LineWidth',2);
    plot(t(:,1),GRF(:,6),'LineWidth',2);
    plot(t(:,1),GRF(:,7),'LineWidth',2);
    plot(t(:,1),GRF(:,8),'LineWidth',2);
    plot(t(:,1),GRF(:,1),'LineWidth',2);
    plot(t(:,1),GRF(:,2),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    ylim([-15 55]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$F_{hip,x}$','$F_{hip,y}$','$F_{heel,x}$','$F_{heel,y}$','$F_{toe,x}$','$F_{toe,y}$','$F_{pull,x}$','$F_{pull,y}$', 'Interpreter', 'latex', 'FontSize',12,'Location','northeast', 'NumColumns',4);
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_externalforce.pdf'], 'ContentType', 'vector');
    end

    %Figure18：目標CFL長さと実際のCFL長さをそれぞれ表示
    figure(18)
    plot(t(:,1),data_l_CFL(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),data_l_CFL(:,1)+data_error_CFL(:,1),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$target_{CFL}$','$l_{CFL}$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    xlabel('Time [s]');
    ylabel('Length [m]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_CFLlength.pdf'], 'ContentType', 'vector');
    end
    %Figure17：外力をそれぞれ表示
    figure(19)
    plot(t(:,1),GRF(:,3),'LineWidth',2);
    hold on
    plot(t(:,1),GRF(:,4),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    ylim([-15 55]);
    % legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    legend('$F_{hip,x}$','$F_{hip,y}$','Interpreter', 'latex', 'FontSize',12,'Location','northeast');
    xlabel('Time [s]');
    ylabel('Force [N]');
    if graph_save == true
        exportgraphics(gca, [save_path name '_externalforce.pdf'], 'ContentType', 'vector');
    end

    
end

%アニメーションの描写
if movie_save == true
    disp('Writing movie...')

    % ====== 動画出力設定 ======
    v = VideoWriter(new_filename, 'MPEG-4');
    v.Quality = 100;
    v.FrameRate = 100;
    open(v);

    % ====== 表示範囲 ======
    x_lim_neg = -0.4;
    x_lim_pos =  0.4;
    y_lim_neg = -0.4;
    y_lim_pos =  0.4;

    % ====== Figure 設定 ======
    fig = figure('Color','w');
    set(fig, 'Units','pixels','Position',[200,200,600,600]); % 正方形の固定ウィンドウ
    ax = axes('Parent', fig);
    hold(ax, 'on');
    grid(ax, 'on');

    % 固定軸設定（重要）
    xlim(ax, [x_lim_neg x_lim_pos]);
    ylim(ax, [y_lim_neg y_lim_pos]);
    axis(ax, 'manual');
    daspect(ax, [1 1 1]);        % データ単位の縦横比を固定
    pbaspect(ax, [1 1 1]);       % 表示上のアスペクトも固定
    axis(ax, 'tight');
    box(ax, 'on');

    % % 軸目盛を非表示
    % ax.XTick = [];
    % ax.YTick = [];
    % 
    % % ====== ground line ======
    % plot(ax, [x_lim_neg x_lim_pos], [y_fixed(1) y_fixed(1)], ':k', 'LineWidth', 1);
    % plot(ax, [x_lim_neg x_lim_pos], [y_fixed(2) y_fixed(2)], ':k', 'LineWidth', 1);

    nFrames = size(q,1);
    % ====== メインループ ======
    for i = 1:size(q,1)
        % 前フレームの描画を消去（地面は残す）
        delete(findall(ax, 'Type', 'line', '-not', 'LineStyle', ':'));

        % ---- リンク描画 ----
        plot(ax, [coordinates_x(i,1) coordinates_x(i,2)], [coordinates_y(i,1) coordinates_y(i,2)], '-k', 'LineWidth', 2);
        plot(ax, [coordinates_x(i,2) coordinates_x(i,3)], [coordinates_y(i,2) coordinates_y(i,3)], '-k', 'LineWidth', 2);
        plot(ax, [coordinates_x(i,3) coordinates_x(i,4)], [coordinates_y(i,3) coordinates_y(i,4)], '-k', 'LineWidth', 2);
        plot(ax, [coordinates_x(i,4) coordinates_x(i,5)], [coordinates_y(i,4) coordinates_y(i,5)], '-k', 'LineWidth', 2);

        % ---- 筋線描画 ----
        plot(ax, [coordinates_x(i,1) coordinates_x(i,6)], [coordinates_y(i,1) coordinates_y(i,6)], '-r', 'LineWidth', 2, 'Color', CFL_Color);
        plot(ax, [coordinates_x(i,6) coordinates_x(i,10)], [coordinates_y(i,6) coordinates_y(i,10)], '-r', 'LineWidth', 2, 'Color', Ci_Color);
        if noCFLT == false
            plot(ax, [coordinates_x(i,6) coordinates_x(i,7)], [coordinates_y(i,6) coordinates_y(i,7)], '-r', 'LineWidth', 2, 'Color', CFLT_Color);
        end
        plot(ax, [coordinates_x(i,7) coordinates_x(i,11)], [coordinates_y(i,7) coordinates_y(i,11)], '-r', 'LineWidth', 2, 'Color', GEo_Color);

        % ---- プーリ描画 ----
        plot(ax, r*cos(phi)+coordinates_x(i,4), r*sin(phi)+coordinates_y(i,4), '-k', 'LineWidth', 2);
        if noGE == false
            plot(ax, [coordinates_x(i,7) coordinates_x(i,8)], [coordinates_y(i,7) coordinates_y(i,8)], '-r', 'LineWidth', 2, 'Color', GE_Color);
        end

        % ---- 重心描写 ----
        plot(ax, COM(i,1), COM(i,2), "o",'MarkerSize', 10,'MarkerFaceColor', 'k','MarkerEdgeColor', 'k');


        % ---- 縦横比・表示範囲を固定（フレームごと確認）----
        xlim(ax, [x_lim_neg x_lim_pos]);
        ylim(ax, [y_lim_neg y_lim_pos]);
        daspect(ax, [1 1 1]);
        pbaspect(ax, [1 1 1]);
        axis(ax, 'manual');

        % ---- フレーム番号を動画内に表示 ----
        txt_x = x_lim_neg + 0.02*(x_lim_pos - x_lim_neg);  % 左下に配置
        txt_y = y_lim_neg + 0.05*(y_lim_pos - y_lim_neg);
        text(txt_x, txt_y, sprintf('Frame %d / %d', i, nFrames), ...
            'FontSize', 12, 'Color', [0 0 0], 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.6]);

        % ---- 描画とフレーム取得 ----
        drawnow limitrate;
        frame = getframe(gcf); % Figure全体をキャプチャ
        writeVideo(v, frame);
    end

    close(v);
    disp('Movie saved successfully.');
end


% figure()
% plot(t(:,1),data_Q(:,8),'LineWidth',2);
% xlabel('Time [s]');
% ylabel('Q7');
% xlim([0 time_lim]);

% CFLとCiの角度
function alpha_Ci = calc_alpha_Ci(coordinates_x, coordinates_y)
    alpha_Ci = zeros(size(coordinates_x,1),1);
    for i = 1:size(coordinates_x,1)
        vec_m2_to_m1 = [coordinates_x(i,1)-coordinates_x(i,6); coordinates_y(i,1)-coordinates_y(i,6)]; %CFL
        vec_m2_to_4th_troch = [coordinates_x(i,10)-coordinates_x(i,6); coordinates_y(i,10)-coordinates_y(i,6)]; %Ci
        alpha_Ci(i) = rad2deg(acos(dot(vec_m2_to_m1,vec_m2_to_4th_troch)/(norm(vec_m2_to_m1)*norm(vec_m2_to_4th_troch))));
    end
end

% CFLとCFLTの角度
function alpha_CFLT = calc_alpha_CFLT(coordinates_x, coordinates_y)
    alpha_CFLT = zeros(size(coordinates_x,1),1);
    for i = 1:size(coordinates_x,1)
        vec_m2_to_m1 = [coordinates_x(i,1)-coordinates_x(i,6); coordinates_y(i,1)-coordinates_y(i,6)]; %CFL
        vec_m2_to_m3 = [coordinates_x(i,7)-coordinates_x(i,6); coordinates_y(i,7)-coordinates_y(i,6)]; %CFLT
        alpha_CFLT(i) = rad2deg(acos(dot(vec_m2_to_m1,vec_m2_to_m3)/(norm(vec_m2_to_m1)*norm(vec_m2_to_m3))));
    end
end

% CFLTとGEoの角度
function alpha_GEo = calc_alpha_GEo(coordinates_x, coordinates_y)
    alpha_GEo = zeros(size(coordinates_x,1),1);
    for i = 1:size(coordinates_x,1)
        vec_m3_to_m2 = [coordinates_x(i,6)-coordinates_x(i,7); coordinates_y(i,6)-coordinates_y(i,7)]; %CFLT
        vec_m3_to_GE_origin = [coordinates_x(i,11)-coordinates_x(i,7); coordinates_y(i,11)-coordinates_y(i,7)]; %GEo
        alpha_GEo(i) = rad2deg(acos(dot(vec_m3_to_m2,vec_m3_to_GE_origin)/(norm(vec_m3_to_m2)*norm(vec_m3_to_GE_origin))));
    end
end

% CFLTとGEの角度
function alpha_GE = calc_alpha_GE(coordinates_x, coordinates_y)
    alpha_GE = zeros(size(coordinates_x,1),1);
    for i = 1:size(coordinates_x,1)
        vec_m3_to_m2 = [coordinates_x(i,6)-coordinates_x(i,7); coordinates_y(i,6)-coordinates_y(i,7)]; %CFLT
        vec_m3_to_pulley = [coordinates_x(i,8)-coordinates_x(i,7); coordinates_x(i,8)-coordinates_y(i,7)]; %GE
        alpha_GE(i) = rad2deg(acos(dot(vec_m3_to_m2,vec_m3_to_pulley)/(norm(vec_m3_to_m2)*norm(vec_m3_to_pulley))));
    end
end

% GEoとfemurの角度
function angle_femur_GEo = calc_angle_femur_GEo(coordinates_x, coordinates_y)
    angle_femur_GEo = zeros(size(coordinates_x,1),1);
    for i = 1:size(coordinates_x,1)
        vec_GEorigin_to_hip = [coordinates_x(i,2)-coordinates_x(i,11); coordinates_y(i,2)-coordinates_y(i,11)]; %femur
        vec_GEorigin_to_m3 = [coordinates_x(i,7)-coordinates_x(i,11); coordinates_y(i,7)-coordinates_y(i,11)]; %GEo
        angle_femur_GEo(i) = rad2deg(acos(dot(vec_GEorigin_to_hip,vec_GEorigin_to_m3)/(norm(vec_GEorigin_to_hip)*norm(vec_GEorigin_to_m3))));
    end
end