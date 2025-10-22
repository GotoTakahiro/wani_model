% clear;
% close all;
% clearvars
% load('results/20240712_MuscleLengthTest_PID/exp20240712_MuscleLengthTest_PID_81_P50000_I50_D550_CFL350_Ci36_CFLT98_GEo41_GE176.mat');
% load('results/20240726_init_condition_test_per2mm/exp20240726_init_condition_test_per2mm_223_Hip20_Knee44_CFL350_Ci44_CFLT107_GEo37_GE188.mat');
load('results/20240822_for_nolta_paper_rev/20241024/exp20240822_for_nolta_paper_knee83_CFL350_Ci44_CFLT100_GEo35_GE185.mat')
[filepath,name,ext] = fileparts(filename);
% save_path = 'results/20240712_MuscleLengthTest_PID/';
save_path = 'results/20240822_for_nolta_paper_rev/20241024/';
% save_path = '';
new_filename = fullfile([save_path name '.mp4']);
graph_save = true;
graph_view = false;
movie_save = true;
time_lim = max(t(:,1));

noCFLT = false;
noGE = false;
if contains(name,'noCFLT')
    noCFLT = true;
end
if contains(name,'noGE')
    noGE = true;
end

phi = linspace(0,2*pi,100);
r = l_link_list(7);

[coordinates_x, coordinates_y, angle_wire_pulley] = calc_coordinate_for_plot(q,r,l_link_list); %1:hip, 2:4th trochanter, 3:GE origin, 4:knee, 5:ankle, 6:toem 7:CFTLT branch, 8:Y shaped branch
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
    general_q = q(j,1:10).';
    COM(j,:) = (calc_COM(m_list,l_link_list,general_q))';
end

for i = 1:size(q,1)
    k_wire = data_k_c_wire(i,2:5);
    c_wire = data_k_c_wire(i,6:10);
    % L_wire(1) = data_acce_GRF(i,18);

    general_q = q(i,1:10).';
    general_dq = q(i,11:20).';

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

coordinates_x_for_angle = coordinates_x(tension_start_index:end,:);
coordinates_y_for_angle = coordinates_y(tension_start_index:end,:);

power = q(50:end,20).*data_Q(50:end,11);
work = trapz(t(50:end,1),power);
disp(['Work: ', num2str(work)]);

% プロットの色を固定
CFL_Color = '#0072BD';
Ci_Color = '#D95319';
CFLT_Color = '#EDB120';
GEo_Color = '#7E2F8E';
GE_Color = '#77AC30';

phase1_DB_t = 4.35;
phase2_DB_t = 5.2;
phase3_DB_t = 5.75;
phase4_DB_t = 8.5;
phase1_DB_index = find(t(:,1) > phase1_DB_t,1);
phase2_DB_index = find(t(:,1) > phase2_DB_t,1);
phase3_DB_index = find(t(:,1) > phase3_DB_t,1);
phase4_DB_index = find(t(:,1) > phase4_DB_t,1);

phase1_CFLTE_t = 5.5;
phase2_CFLTE_t = 7;
phase1_CFLTE_index = find(t(:,1) > phase1_CFLTE_t,1);
phase2_CFLTE_index = find(t(:,1) > phase2_CFLTE_t,1);

phase1_GEE_t = 4.5;
phase2_GEE_t = 5.35;
phase1_GEE_index = find(t(:,1) > phase1_GEE_t,1);
phase2_GEE_index = find(t(:,1) > phase2_GEE_t,1);

if graph_view == true

    % 張力
    figure(1)
    plot(t(:,1),-data_Q(:,11),'-','LineWidth',2,'Color',CFL_Color); %CFL
    hold on
    plot(t(:,1),muscle_tension(:,1),'--','LineWidth',2,'Color', Ci_Color); %Ci
    if noCFLT == false
        plot(t(:,1),muscle_tension(:,2),'-.','LineWidth',2,'Color', CFLT_Color); %CFLT
    end
    plot(t(:,1),muscle_tension(:,3),':','LineWidth',2,'Color',GEo_Color); %GEo
    if noGE == false
        plot(t(:,1),muscle_tension(:,4),'LineWidth',2,'Color',GE_Color); %GE
    end
    xlim([0 time_lim]);
    ylim([-10 500]);
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
    %legendを少し右にずらす
    % lgd = legend;
    % newPosition = [0.7 0.7 0.1 0.1];
    % lgd.Position = newPosition;
    % legend('CFL','Ci','CFLT','GEo','GE','FontSize',20,'Location','best');
    % 各領域に色を指定して背景色をつける (RGB指定)
    if noCFLT == false && noGE == false
        x_region0 = [0 t(phase1_DB_index,1) t(phase1_DB_index,1) 0];
        x_region1 = [t(phase1_DB_index,1) t(phase2_DB_index,1) t(phase2_DB_index,1) t(phase1_DB_index,1)];
        x_region2 = [t(phase2_DB_index,1) t(phase3_DB_index,1) t(phase3_DB_index,1) t(phase2_DB_index,1)];
        x_region3 = [t(phase3_DB_index,1) t(phase4_DB_index,1) t(phase4_DB_index,1) t(phase3_DB_index,1)];
        x_region4 = [t(phase4_DB_index,1) t(1000,1) t(1000,1) t(phase4_DB_index,1)];
        y_region = [-10 -10 500 500];
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        p3 = patch(x_region3, y_region, [0.7, 0.7, 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い青
        p4 = patch(x_region4, y_region, [0.9, 0.7, 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  %薄いオレンジ
        hold off
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        uistack(p2, 'bottom');
        uistack(p3, 'bottom');
        uistack(p4, 'bottom');
    elseif noCFLT == true && noGE == false
        x_region0 = [0 t(phase1_CFLTE_index,1) t(phase1_CFLTE_index,1) 0];
        x_region1 = [t(phase1_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase1_CFLTE_index,1)];
        % x_region2 = [t(phase2_CFLTE_index,1) t(1000,1) t(1000,1) t(phase2_CFLTE_index,1)];
        y_region = [-10 -10 500 500];
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        hold off
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        % uistack(p2, 'bottom');
    elseif noCFLT == false && noGE == true
        x_region0 = [0 t(phase1_GEE_index,1) t(phase1_GEE_index,1) 0];
        x_region1 = [t(phase1_GEE_index,1) t(phase2_GEE_index,1) t(phase2_GEE_index,1) t(phase1_GEE_index,1)];
        y_region = [-10 -10 500 500];
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        hold off
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
    end
    if graph_save == true
        exportgraphics(gca, [save_path name '_tension.pdf'], 'ContentType', 'vector');
    end

    % トルク
    figure(2)
    plot(t(:,1),torque_muscle(:,6),'LineWidth',2);
    hold on
    plot(t(:,1),torque_muscle(:,7),'LineWidth',2);
    plot(t(:,1),torque_muscle(:,8),'LineWidth',2);

    if noCFLT == false && noGE == false
        x_region0 = [0 t(phase1_DB_index,1) t(phase1_DB_index,1) 0];
        x_region1 = [t(phase1_DB_index,1) t(phase2_DB_index,1) t(phase2_DB_index,1) t(phase1_DB_index,1)];
        x_region2 = [t(phase2_DB_index,1) t(phase3_DB_index,1) t(phase3_DB_index,1) t(phase2_DB_index,1)];
        x_region3 = [t(phase3_DB_index,1) t(phase4_DB_index,1) t(phase4_DB_index,1) t(phase3_DB_index,1)];
        x_region4 = [t(phase4_DB_index,1) t(1000,1) t(1000,1) t(phase4_DB_index,1)];
        y_region = [-15 -15 2 2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        p3 = patch(x_region3, y_region, [0.7, 0.7, 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い青
        p4 = patch(x_region4, y_region, [0.9, 0.7, 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  %薄いオレンジ
        hold off
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        uistack(p2, 'bottom');
        uistack(p3, 'bottom');
        uistack(p4, 'bottom');
    elseif noCFLT == true && noGE == false
        x_region0 = [0 t(phase1_CFLTE_index,1) t(phase1_CFLTE_index,1) 0];
        x_region1 = [t(phase1_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase1_CFLTE_index,1)];
        % x_region2 = [t(phase2_CFLTE_index,1) t(1000,1) t(1000,1) t(phase2_CFLTE_index,1)];
        y_region = [-15 -15 2 2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        hold off
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        % uistack(p2, 'bottom');
    elseif noCFLT == false && noGE == true
        x_region0 = [0 t(phase1_GEE_index,1) t(phase1_GEE_index,1) 0];
        x_region1 = [t(phase1_GEE_index,1) t(phase2_GEE_index,1) t(phase2_GEE_index,1) t(phase1_GEE_index,1)];
        y_region = [-15 -15 2 2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        hold off
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
    end

    hold off
    xlim([0 time_lim]);
    ylim([-15 2]);
    h_axes = gca;
    h_axes.XAxis.FontSize = 20;
    h_axes.YAxis.FontSize = 20;
    xlabel('Time [s]','FontSize',25);
    ylabel('Torque [Nm]','FontSize',25);
    legend('hip','knee','ankle','FontSize',25, 'Location','best');
    if graph_save == true
        exportgraphics(gca, [save_path name '_torque.pdf'], 'ContentType', 'vector');
    end

    
    %関節角度
    figure(3)
    plot(t(:,1),rad2deg(q(:,6)),'LineWidth',2);
    hold on
    plot(t(:,1),rad2deg(q(:,7)),'LineWidth',2);
    plot(t(:,1),rad2deg(q(:,8)),'LineWidth',2);
    
    if noCFLT == false && noGE == false
        x_region0 = [0 t(phase1_DB_index,1) t(phase1_DB_index,1) 0];
        x_region1 = [t(phase1_DB_index,1) t(phase2_DB_index,1) t(phase2_DB_index,1) t(phase1_DB_index,1)];
        x_region2 = [t(phase2_DB_index,1) t(phase3_DB_index,1) t(phase3_DB_index,1) t(phase2_DB_index,1)];
        x_region3 = [t(phase3_DB_index,1) t(phase4_DB_index,1) t(phase4_DB_index,1) t(phase3_DB_index,1)];
        x_region4 = [t(phase4_DB_index,1) t(1000,1) t(1000,1) t(phase4_DB_index,1)];
        y_region = [-200 -200 150 150];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        p3 = patch(x_region3, y_region, [0.7, 0.7, 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い青
        p4 = patch(x_region4, y_region, [0.9, 0.7, 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  %薄いオレンジ
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        uistack(p2, 'bottom');
        uistack(p3, 'bottom');
        uistack(p4, 'bottom');
    elseif noCFLT == true && noGE == false
        x_region0 = [0 t(phase1_CFLTE_index,1) t(phase1_CFLTE_index,1) 0];
        x_region1 = [t(phase1_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase1_CFLTE_index,1)];
        % x_region2 = [t(phase2_CFLTE_index,1) t(1000,1) t(1000,1) t(phase2_CFLTE_index,1)];
        y_region = [-200 -200 150 150];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        % uistack(p2, 'bottom');
    elseif noCFLT == false && noGE == true
        x_region0 = [0 t(phase1_GEE_index,1) t(phase1_GEE_index,1) 0];
        x_region1 = [t(phase1_GEE_index,1) t(phase2_GEE_index,1) t(phase2_GEE_index,1) t(phase1_GEE_index,1)];
        y_region = [-200 -200 150 150];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
    end
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

    %地面反力
    GRF = data_acce_GRF(:,14:17);
    figure(4)
    plot(t(:,1),GRF(:,1),'LineWidth',2);
    hold on
    plot(t(:,1),GRF(:,2),'LineWidth',2);
    plot(t(:,1),GRF(:,3),'LineWidth',2);
    plot(t(:,1),GRF(:,4),'LineWidth',2);
    % plot(t(:,1),(GRF(:,1)+GRF(:,3)),'LineWidth',2);
    % plot(data_acce_GRF(:,1),(GRF(:,2)+GRF(:,4)),'LineWidth',2);

    
    if noCFLT == false && noGE == false
        x_region0 = [0 t(phase1_DB_index,1) t(phase1_DB_index,1) 0];
        x_region1 = [t(phase1_DB_index,1) t(phase2_DB_index,1) t(phase2_DB_index,1) t(phase1_DB_index,1)];
        x_region2 = [t(phase2_DB_index,1) t(phase3_DB_index,1) t(phase3_DB_index,1) t(phase2_DB_index,1)];
        x_region3 = [t(phase3_DB_index,1) t(phase4_DB_index,1) t(phase4_DB_index,1) t(phase3_DB_index,1)];
        x_region4 = [t(phase4_DB_index,1) t(1000,1) t(1000,1) t(phase4_DB_index,1)];
        y_region = [-15 -15 55 55];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        p3 = patch(x_region3, y_region, [0.7, 0.7, 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い青
        p4 = patch(x_region4, y_region, [0.9, 0.7, 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  %薄いオレンジ
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        uistack(p2, 'bottom');
        uistack(p3, 'bottom');
        uistack(p4, 'bottom');
    elseif noCFLT == true && noGE == false
        x_region0 = [0 t(phase1_CFLTE_index,1) t(phase1_CFLTE_index,1) 0];
        x_region1 = [t(phase1_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase1_CFLTE_index,1)];
        % x_region2 = [t(phase2_CFLTE_index,1) t(1000,1) t(1000,1) t(phase2_CFLTE_index,1)];
        y_region = [-15 -15 55 55];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        % uistack(p2, 'bottom');
    elseif noCFLT == false && noGE == true
        x_region0 = [0 t(phase1_GEE_index,1) t(phase1_GEE_index,1) 0];
        x_region1 = [t(phase1_GEE_index,1) t(phase2_GEE_index,1) t(phase2_GEE_index,1) t(phase1_GEE_index,1)];
        y_region = [-15 -15 55 55];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
    end

    hold off
    xlim([0 time_lim]);
    ylim([-15 55]);
    h_axes = gca;
    h_axes.XAxis.FontSize = 20;
    h_axes.YAxis.FontSize = 20;
    % legend('Heel-Fx','Heel-Fy','Toe-Fx','Toe-y','Total-Fx','Total-Fy','FontSize',20,'Location','best');
    legend('Heel-$F_x$','Heel-$F_y$','Toe-$F_x$','Toe-$F_y$', 'Interpreter', 'latex', 'FontSize',20,'Location','best');
    xlabel('Time [s]','FontSize',25);
    ylabel('Force [N]', 'FontSize',25);
    if graph_save == true
        exportgraphics(gca, [save_path name '_GRF.pdf'], 'ContentType', 'vector');
    end

    figure(5)
    plot(t(50:end,1),power,'LineWidth',2);
    xlabel('Time [s]');
    ylabel('Power [W]');
    xlim([0 time_lim]);
    if graph_save == true
        exportgraphics(gca, [save_path name '_power.pdf'], 'ContentType', 'vector');
    end

    figure(6)
    plot(t(:,1), stationary_error_CFL(:,1),'LineWidth',2);
    xlabel('Time [s]');
    ylabel('Error');
    xlim([0 time_lim]);
    ylim([-0.01 0.01]);
    if graph_save == true
        exportgraphics(gca, [save_path name '_error.pdf'], 'ContentType', 'vector');
    end

    angle_femur_GEo = calc_angle_femur_GEo(coordinates_x_for_angle, coordinates_y_for_angle);
    figure(7)
    yyaxis left
    plot(t(tension_start_index:end,1)-t(tension_start_index,1),rad2deg(q(tension_start_index:end,6)),'LineWidth',2);
    ylabel('Hip angle [deg]');
    ylim([-200 150]);
    yyaxis right
    plot(t(tension_start_index:end,1)-t(tension_start_index,1),angle_femur_GEo,'LineWidth',2);
    xlabel('Time [s]');
    ylabel('GEo angle [deg]');
    xlim([0 time_lim-t(tension_start_index,1)]);
    ylim([0 180]);
    if graph_save == true
        exportgraphics(gca, [save_path name '_angle_femur_GEo.pdf'], 'ContentType', 'vector');
    end

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

    figure(10)
    % plot(t(:,1),momentum_list(:,1),'LineWidth',2);
    % hold on
    plot(t(:,1),momentum_list(:,2),'LineWidth',2);
    hold on
    plot(t(:,1),momentum_list(:,3),'LineWidth',2);
    plot(t(:,1),momentum_list(:,4),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    ylim([-0.002 0.002]);
    xlabel('Time [s]');
    ylabel('Momentum [J・s]');
    % legend('frame','femur','tibia','metatarsal');
    legend('femur','tibia','metatarsal');
    if graph_save == true
        exportgraphics(gca, [save_path name '_momentum.pdf'], 'ContentType', 'vector');
    end

    figure(11)
    % plot(t(:,1),angular_moment_list(:,1),'LineWidth',2);
    % hold on
    plot(t(:,1),angular_moment_list(:,2),'LineWidth',2);
    hold on
    plot(t(:,1),angular_moment_list(:,3),'LineWidth',2);
    plot(t(:,1),angular_moment_list(:,4),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    ylim([-0.0007 0.0007]);
    xlabel('Time [s]');
    ylabel('Angular momentum [J・s]');
    % legend('frame','femur','tibia','metatarsal');
    legend('femur','tibia','metatarsal');
    if graph_save == true
        exportgraphics(gca, [save_path name '_angular_momentum.pdf'], 'ContentType', 'vector');
    end

    figure(12)
    % plot(t(:,1),momentum_COM_list(:,1),'LineWidth',2);
    % hold on
    % plot(t(:,1),momentum_COM_list(:,2),'LineWidth',2);
    plot(t(:,1),momentum_COM_list(:,3),'LineWidth',2);
    hold on
    plot(t(:,1),momentum_COM_list(:,4),'LineWidth',2);
    plot(t(:,1),momentum_COM_list(:,5),'LineWidth',2);
    plot(t(:,1),momentum_COM_list(:,6),'LineWidth',2);
    plot(t(:,1),momentum_COM_list(:,7),'LineWidth',2);
    plot(t(:,1),momentum_COM_list(:,8),'LineWidth',2);
    hold off
    xlim([0 time_lim]);
    ylim([-0.04 0.04]);
    xlabel('Time [s]');
    ylabel('Momentum [kg・m/s]');
    % legend('caudal x','caudal y','femur x','femur y','tibia x','tibia y','metatarsal x','metatarsal y');
    legend('femur x','femur y','tibia x','tibia y','metatarsal x','metatarsal y');
    if graph_save == true
        exportgraphics(gca, [save_path name '_momentum_COM.pdf'], 'ContentType', 'vector');
    end

    % angle_femur_GEo = calc_angle_femur_GEo(coordinates_x, coordinates_y);
    figure(13)
    plot(angular_moment_list(tension_start_index:end,2), angular_moment_list(tension_start_index:end,3),'LineWidth',2);
    xlabel('Femur angular momentum [J・s]');
    xlim([-0.0007 0.0007]);
    ylabel('Tibia angular momentum [J・s]');
    ylim([-0.0007 0.0007]);
    if graph_save == true
        exportgraphics(gca, [save_path name '_angular_momentum_femur_tibia.pdf'], 'ContentType', 'vector');
    end
    
    % angle_femur_GEo = calc_angle_femur_GEo(coordinates_x_for_angle, coordinates_y_for_angle);
    figure(14)
    plot( angle_femur_GEo, rad2deg(q(tension_start_index:end,6)), 'LineWidth',2);
    xlabel('GEo angle [deg]');
    xlim([0 180]);
    ylabel('Hip angle [deg]');
    ylim([-130 0]);
    if graph_save == true
        exportgraphics(gca, [save_path name '_angle_femur_GEo_rev.pdf'], 'ContentType', 'vector');
    end

    figure(15)
    plot(rad2deg(q(:,6)),-data_Q(:,11),'-','LineWidth',2,'Color','#0072BD') %CFL
    hold on
    plot(rad2deg(q(:,6)),muscle_tension(:,1),'--','LineWidth',2,'Color', '#D95319') %Ci
    if noCFLT == false
        plot(rad2deg(q(:,6)),muscle_tension(:,2),'-.','LineWidth',2,'Color', '#EDB120') %CFLT
    end
    plot(rad2deg(q(:,6)),muscle_tension(:,3),':','LineWidth',2,'Color','#7E2F8E') %GEo
    if noGE == false
        plot(rad2deg(q(:,6)),muscle_tension(:,4),'LineWidth',2,'Color','#77AC30') %GE
    end

    
    if noCFLT == false && noGE == false
        x_region0 = [-5 rad2deg(q(phase1_DB_index,6)) rad2deg(q(phase1_DB_index,6)) -5];
        x_region1 = [rad2deg(q(phase1_DB_index,6)) rad2deg(q(phase2_DB_index,6)) rad2deg(q(phase2_DB_index,6)) rad2deg(q(phase1_DB_index,6))];
        x_region2 = [rad2deg(q(phase2_DB_index,6)) rad2deg(q(phase3_DB_index,6)) rad2deg(q(phase3_DB_index,6)) rad2deg(q(phase2_DB_index,6))];
        x_region3 = [rad2deg(q(phase3_DB_index,6)) rad2deg(q(phase4_DB_index,6)) rad2deg(q(phase4_DB_index,6)) rad2deg(q(phase3_DB_index,6))];
        x_region4 = [rad2deg(q(phase4_DB_index,6)) rad2deg(q(1000,6)) rad2deg(q(1000,6)) rad2deg(q(phase4_DB_index,6))];
        y_region = [-10 -10 500 500];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        p3 = patch(x_region3, y_region, [0.7, 0.7, 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い青
        p4 = patch(x_region4, y_region, [0.9, 0.7, 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  %薄いオレンジ
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        uistack(p2, 'bottom');
        uistack(p3, 'bottom');
        uistack(p4, 'bottom');
    elseif noCFLT == true && noGE == false
        x_region0 = [0 rad2deg(q(phase1_CFLTE_index,6)) rad2deg(q(phase1_CFLTE_index,6)) 0];
        x_region1 = [rad2deg(q(phase1_CFLTE_index,6)) rad2deg(q(phase2_CFLTE_index,6)) rad2deg(q(phase2_CFLTE_index,6)) rad2deg(q(phase1_CFLTE_index,6))];
        % x_region2 = [rad2deg(q(phase2_CFLTE_index,6)) rad2deg(q(1000,6)) rad2deg(q(1000,6)) rad2deg(q(phase2_CFLTE_index,6))];
        y_region = [-10 -10 500 500];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        % uistack(p2, 'bottom');
    elseif noCFLT == false && noGE == true
        x_region0 = [0 rad2deg(q(phase1_GEE_index,6)) rad2deg(q(phase1_GEE_index,6)) 0];
        x_region1 = [rad2deg(q(phase1_GEE_index,6)) rad2deg(q(phase2_GEE_index,6)) rad2deg(q(phase2_GEE_index,6)) rad2deg(q(phase1_GEE_index,6))];
        y_region = [-10 -10 500 500];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
    end

    hold off
    xlabel('Hip angle [deg]','FontSize',25);
    ylabel('Tension [N]','FontSize',25);
    xlim([-95 -5]);
    ylim([-10 500]);
    h_axes = gca;
    h_axes.XAxis.FontSize = 20;
    h_axes.YAxis.FontSize = 20;
    % legend('CFL','Ci','CFLT','GEo','GE','FontSize',20,'Location','best');
    if noCFLT == false && noGE == false
        % disp('CFLT and GE');
        legend('CFL','Ci','CFLT','GEo','GE','FontSize',20,'Location','best');
        % 判例の位置を少し右下に移動
        lgd = legend;
        % 判例の現在の位置を取得
        currentPosition = lgd.Position;
        % Y位置を少し下げる
        newPosition = currentPosition;
        newPosition(1) = newPosition(1) + 0.1;  % x座標
        newPosition(2) = newPosition(2) - 0.1;  % y座標
        lgd.Position = newPosition;

    elseif noCFLT == false && noGE == true
        legend('CFL','Ci','CFLT','GEo','FontSize',20,'Location','best');
    elseif noCFLT == true && noGE == false
        legend('CFL','Ci','GEo','GE','FontSize',20,'Location','best');
    end
    % lgd = legend;
    % % 判例の現在の位置を取得
    % currentPosition = lgd.Position;
    % % Y位置を少し下げる
    % newPosition = currentPosition;
    % newPosition(1) = newPosition(1) + 0.3;  % x座標
    % newPosition(2) = newPosition(2) - 0.1;  % y座標
    % lgd.Position = newPosition;
    set(gca, 'XDir', 'reverse');
    if graph_save == true
        exportgraphics(gca, [save_path name '_tension_hip.pdf'], 'ContentType', 'vector');
    end

    figure(16)
    plot(rad2deg(q(:,6)),torque_muscle(:,6),'LineWidth',2);
    hold on
    plot(rad2deg(q(:,6)),torque_muscle(:,7),'LineWidth',2);
    plot(rad2deg(q(:,6)),torque_muscle(:,8),'LineWidth',2);

    
    if noCFLT == false && noGE == false
        x_region0 = [-5 rad2deg(q(phase1_DB_index,6)) rad2deg(q(phase1_DB_index,6)) -5];
        x_region1 = [rad2deg(q(phase1_DB_index,6)) rad2deg(q(phase2_DB_index,6)) rad2deg(q(phase2_DB_index,6)) rad2deg(q(phase1_DB_index,6))];
        x_region2 = [rad2deg(q(phase2_DB_index,6)) rad2deg(q(phase3_DB_index,6)) rad2deg(q(phase3_DB_index,6)) rad2deg(q(phase2_DB_index,6))];
        x_region3 = [rad2deg(q(phase3_DB_index,6)) rad2deg(q(phase4_DB_index,6)) rad2deg(q(phase4_DB_index,6)) rad2deg(q(phase3_DB_index,6))];
        x_region4 = [rad2deg(q(phase4_DB_index,6)) rad2deg(q(1000,6)) rad2deg(q(1000,6)) rad2deg(q(phase4_DB_index,6))];
        y_region = [-15 -15 2 2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        p3 = patch(x_region3, y_region, [0.7, 0.7, 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い青
        p4 = patch(x_region4, y_region, [0.9, 0.7, 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  %薄いオレンジ
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        uistack(p2, 'bottom');
        uistack(p3, 'bottom');
        uistack(p4, 'bottom');
    elseif noCFLT == true && noGE == false
        x_region0 = [0 rad2deg(q(phase1_CFLTE_index,6)) rad2deg(q(phase1_CFLTE_index,6)) 0];
        x_region1 = [rad2deg(q(phase1_CFLTE_index,6)) rad2deg(q(phase2_CFLTE_index,6)) rad2deg(q(phase2_CFLTE_index,6)) rad2deg(q(phase1_CFLTE_index,6))];
        % x_region2 = [rad2deg(q(phase2_CFLTE_index,6)) rad2deg(q(1000,6)) rad2deg(q(1000,6)) rad2deg(q(phase2_CFLTE_index,6))];
        y_region = [-15 -15 2 2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        % uistack(p2, 'bottom');
    elseif noCFLT == false && noGE == true
        x_region0 = [0 rad2deg(q(phase1_GEE_index,6)) rad2deg(q(phase1_GEE_index,6)) 0];
        x_region1 = [rad2deg(q(phase1_GEE_index,6)) rad2deg(q(phase2_GEE_index,6)) rad2deg(q(phase2_GEE_index,6)) rad2deg(q(phase1_GEE_index,6))];
        y_region = [-15 -15 2 2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
    end

    hold off
    xlabel('Hip angle [deg]');
    ylabel('Torque [Nm]');
    xlim([-95 -5]);
    ylim([-15 2]);
    h_axes = gca;
    h_axes.XAxis.FontSize = 20;
    h_axes.YAxis.FontSize = 20;
    legend('hip','knee','ankle','FontSize',25, 'Location','best');
    set(gca, 'XDir', 'reverse');
    if graph_save == true
        exportgraphics(gca, [save_path name '_torque_hip.pdf'], 'ContentType', 'vector');
    end

    figure(17)
    x_lim_neg = -0.4;
    x_lim_pos = 0.4;
    y_lim_neg = -0.4;
    y_lim_pos = 0.4;
    plot(t(:,1),coordinates_y(:,2),'LineWidth',2);
    hold on
    plot([0 t(end,1)],[y_fixed(1) y_fixed(1)],':k','LineWidth',2); %ground

    
    if noCFLT == false && noGE == false
        x_region0 = [0 t(phase1_DB_index,1) t(phase1_DB_index,1) 0];
        x_region1 = [t(phase1_DB_index,1) t(phase2_DB_index,1) t(phase2_DB_index,1) t(phase1_DB_index,1)];
        x_region2 = [t(phase2_DB_index,1) t(phase3_DB_index,1) t(phase3_DB_index,1) t(phase2_DB_index,1)];
        x_region3 = [t(phase3_DB_index,1) t(phase4_DB_index,1) t(phase4_DB_index,1) t(phase3_DB_index,1)];
        x_region4 = [t(phase4_DB_index,1) t(1000,1) t(1000,1) t(phase4_DB_index,1)];
        y_region = [-0.1 -0.1 0.2 0.2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        p3 = patch(x_region3, y_region, [0.7, 0.7, 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い青
        p4 = patch(x_region4, y_region, [0.9, 0.7, 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  %薄いオレンジ
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        uistack(p2, 'bottom');
        uistack(p3, 'bottom');
        uistack(p4, 'bottom');
    elseif noCFLT == true && noGE == false
        x_region0 = [0 t(phase1_CFLTE_index,1) t(phase1_CFLTE_index,1) 0];
        x_region1 = [t(phase1_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase1_CFLTE_index,1)];
        % x_region2 = [t(phase2_CFLTE_index,1) t(1000,1) t(1000,1) t(phase2_CFLTE_index,1)];
        y_region = [-0.1 -0.1 0.2 0.2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        % uistack(p2, 'bottom');
    elseif noCFLT == false && noGE == true
        x_region0 = [0 t(phase1_GEE_index,1) t(phase1_GEE_index,1) 0];
        x_region1 = [t(phase1_GEE_index,1) t(phase2_GEE_index,1) t(phase2_GEE_index,1) t(phase1_GEE_index,1)];
        y_region = [-0.1 -0.1 0.2 0.2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
    end

    hold off
    xlim([0 time_lim]);
    ylim([-0.05 0.15]);
    h_axes = gca;
    h_axes.XAxis.FontSize = 20;
    h_axes.YAxis.FontSize = 20;
    xlabel('Time [s]','FontSize',25);
    ylabel('$y$ [m]','Interpreter','latex','FontSize',25);
    text(10.1,y_fixed(1)-0.005,'$y_{lim}$','Interpreter','latex','FontSize',20)
    if graph_save == true
        exportgraphics(gca, [save_path name '_hip_y.pdf'], 'ContentType', 'vector');
    end

    figure(18)
    plot(t(:,1),COM(:,1),'LineWidth',2);
    hold on
    if noCFLT == false && noGE == false
        x_region0 = [0 t(phase1_DB_index,1) t(phase1_DB_index,1) 0];
        x_region1 = [t(phase1_DB_index,1) t(phase2_DB_index,1) t(phase2_DB_index,1) t(phase1_DB_index,1)];
        x_region2 = [t(phase2_DB_index,1) t(phase3_DB_index,1) t(phase3_DB_index,1) t(phase2_DB_index,1)];
        x_region3 = [t(phase3_DB_index,1) t(phase4_DB_index,1) t(phase4_DB_index,1) t(phase3_DB_index,1)];
        x_region4 = [t(phase4_DB_index,1) t(1000,1) t(1000,1) t(phase4_DB_index,1)];
        y_region = [-0.1 -0.1 0.2 0.2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        p3 = patch(x_region3, y_region, [0.7, 0.7, 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い青
        p4 = patch(x_region4, y_region, [0.9, 0.7, 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  %薄いオレンジ
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        uistack(p2, 'bottom');
        uistack(p3, 'bottom');
        uistack(p4, 'bottom');
    elseif noCFLT == true && noGE == false
        x_region0 = [0 t(phase1_CFLTE_index,1) t(phase1_CFLTE_index,1) 0];
        x_region1 = [t(phase1_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase2_CFLTE_index,1) t(phase1_CFLTE_index,1)];
        % x_region2 = [t(phase2_CFLTE_index,1) t(1000,1) t(1000,1) t(phase2_CFLTE_index,1)];
        y_region = [-0.1 -0.1 0.2 0.2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % p2 = patch(x_region2, y_region, [0.7, 0.9, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い緑
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
        % uistack(p2, 'bottom');
    elseif noCFLT == false && noGE == true
        x_region0 = [0 t(phase1_GEE_index,1) t(phase1_GEE_index,1) 0];
        x_region1 = [t(phase1_GEE_index,1) t(phase2_GEE_index,1) t(phase2_GEE_index,1) t(phase1_GEE_index,1)];
        y_region = [-0.1 -0.1 0.2 0.2];
        % 各領域に色を指定して背景色をつける (RGB指定)
        p0 = patch(x_region0, y_region, [0.8, 0.8, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い灰色
        p1 = patch(x_region1, y_region, [0.9, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');  % 薄い赤
        % patchを背面に移動
        uistack(p0, 'bottom');
        uistack(p1, 'bottom');
    end

    xlim([0 time_lim]);
    ylim([-0 0.2]);
    h_axes = gca;
    h_axes.XAxis.FontSize = 20;
    h_axes.YAxis.FontSize = 20;
    xlabel('Time [s]', 'FontSize',25);
    ylabel('$x$ [m]', 'Interpreter','latex', 'FontSize',25);
    if graph_save == true
        exportgraphics(gca, [save_path name '_COM_x.pdf'], 'ContentType', 'vector');
    end

end


if movie_save == true
    % % %結果を動画としてプロット
    clf;
    v = VideoWriter(new_filename,'MPEG-4');
    v.Quality = 100;
    v.FrameRate = 100;
    open(v)

    x_lim_neg = -0.4;
    x_lim_pos = 0.4;
    y_lim_neg = -0.4;
    y_lim_pos = 0.4;

    ax = gca();

    for i=1:size(q,1)
        % plot(coordinates_x(i,1),coordinates_y(i,1),'.', 'MarkerSize', 50); %hip
        % hold on
        plot([coordinates_x(i,1) coordinates_x(i,2)],[coordinates_y(i,1) coordinates_y(i,2)],'-','Color','k','MarkerSize',2,'LineWidth',2); %CFL_origin to hip
        hold on
        plot([coordinates_x(i,2) coordinates_x(i,3)],[coordinates_y(i,2) coordinates_y(i,3)],'-','Color','k','MarkerSize',2,'LineWidth',2); %hip to knee
        plot([coordinates_x(i,3) coordinates_x(i,4)],[coordinates_y(i,3) coordinates_y(i,4)],'-','Color','k','MarkerSize',2,'LineWidth',2); %knee to ankle
        plot([coordinates_x(i,4) coordinates_x(i,5)],[coordinates_y(i,4) coordinates_y(i,5)],'-','Color','k','MarkerSize',2,'LineWidth',2); %ankle to toe
        plot([coordinates_x(i,1) coordinates_x(i,6)],[coordinates_y(i,1) coordinates_y(i,6)],'-r','MarkerSize',2,'LineWidth',2,'Color',CFL_Color); %CFL
        plot([coordinates_x(i,6) coordinates_x(i,10)],[coordinates_y(i,6) coordinates_y(i,10)],'-r','MarkerSize',2,'LineWidth',2,'Color',Ci_Color); %Ci
        if noCFLT == false
            plot([coordinates_x(i,6) coordinates_x(i,7)],[coordinates_y(i,6) coordinates_y(i,7)],'-r','MarkerSize',2,'LineWidth',2,'Color',CFLT_Color); %CFLT
        end
        plot([coordinates_x(i,7) coordinates_x(i,11)],[coordinates_y(i,7) coordinates_y(i,11)],'-r','MarkerSize',2,'LineWidth',2,'Color',GEo_Color); %GEo

        plot(r*cos(phi)+coordinates_x(i,4) , r*sin(phi)+coordinates_y(i,4),'-','Color','k','MarkerSize',2,'LineWidth',2); %pulley
        if noGE == false
            plot([coordinates_x(i,7) coordinates_x(i,8)],[coordinates_y(i,7) coordinates_y(i,8)],'-r','MarkerSize',2,'LineWidth',2,'Color',GE_Color); %GE
        end

        plot([x_lim_neg x_lim_pos],[y_fixed(1) y_fixed(1)],':k','LineWidth',1); %ground
        plot([x_lim_neg x_lim_pos],[y_fixed(2) y_fixed(2)],':k','LineWidth',1); %ground

        plot(COM(i,1),COM(i,2),'o','MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','r'); %COM

        drawnow limitrate % drawnowを入れるとアニメーションになる
        hold off
        grid on
        xlim([x_lim_neg x_lim_pos]);
        ylim([y_lim_neg y_lim_pos]);
        axis square
        frame= getframe(ax); % アニメーションのフレームをゲットする 
        writeVideo(v,frame);

        % 0.25秒ごとにpdfとして保存．座標軸はなし
        ax.XTick = [];
        ax.YTick = [];
        % if mod(i,25) == 0
        %     exportgraphics(gca, [save_path name '_movie_' num2str(i) '.pdf'], 'ContentType', 'vector');
        % end
    end
    close(v)
end

function alpha_Ci = calc_alpha_Ci(coordinates_x, coordinates_y)
    alpha_Ci = zeros(size(coordinates_x,1),1);
    for i = 1:size(coordinates_x,1)
        vec_m2_to_m1 = [coordinates_x(i,1)-coordinates_x(i,6); coordinates_y(i,1)-coordinates_y(i,6)];
        vec_m2_to_4th_troch = [coordinates_x(i,10)-coordinates_x(i,6); coordinates_y(i,10)-coordinates_y(i,6)];
        alpha_Ci(i) = rad2deg(acos(dot(vec_m2_to_m1,vec_m2_to_4th_troch)/(norm(vec_m2_to_m1)*norm(vec_m2_to_4th_troch))));
    end
end

function alpha_CFLT = calc_alpha_CFLT(coordinates_x, coordinates_y)
    alpha_CFLT = zeros(size(coordinates_x,1),1);
    for i = 1:size(coordinates_x,1)
        vec_m2_to_m1 = [coordinates_x(i,1)-coordinates_x(i,6); coordinates_y(i,1)-coordinates_y(i,6)];
        vec_m2_to_m3 = [coordinates_x(i,7)-coordinates_x(i,6); coordinates_y(i,7)-coordinates_y(i,6)];
        alpha_CFLT(i) = rad2deg(acos(dot(vec_m2_to_m1,vec_m2_to_m3)/(norm(vec_m2_to_m1)*norm(vec_m2_to_m3))));
    end
end

function alpha_GEo = calc_alpha_GEo(coordinates_x, coordinates_y)
    alpha_GEo = zeros(size(coordinates_x,1),1);
    for i = 1:size(coordinates_x,1)
        vec_m3_to_m2 = [coordinates_x(i,6)-coordinates_x(i,7); coordinates_y(i,6)-coordinates_y(i,7)];
        vec_m3_to_GE_origin = [coordinates_x(i,11)-coordinates_x(i,7); coordinates_y(i,11)-coordinates_y(i,7)];
        alpha_GEo(i) = rad2deg(acos(dot(vec_m3_to_m2,vec_m3_to_GE_origin)/(norm(vec_m3_to_m2)*norm(vec_m3_to_GE_origin))));
    end
end

function alpha_GE = calc_alpha_GE(coordinates_x, coordinates_y)
    alpha_GE = zeros(size(coordinates_x,1),1);
    for i = 1:size(coordinates_x,1)
        vec_m3_to_m2 = [coordinates_x(i,6)-coordinates_x(i,7); coordinates_y(i,6)-coordinates_y(i,7)];
        vec_m3_to_pulley = [coordinates_x(i,8)-coordinates_x(i,7); coordinates_x(i,8)-coordinates_y(i,7)];
        alpha_GE(i) = rad2deg(acos(dot(vec_m3_to_m2,vec_m3_to_pulley)/(norm(vec_m3_to_m2)*norm(vec_m3_to_pulley))));
    end
end

function angle_femur_GEo = calc_angle_femur_GEo(coordinates_x, coordinates_y)
    angle_femur_GEo = zeros(size(coordinates_x,1),1);
    for i = 1:size(coordinates_x,1)
        vec_GEorigin_to_hip = [coordinates_x(i,2)-coordinates_x(i,11); coordinates_y(i,2)-coordinates_y(i,11)];
        vec_GEorigin_to_m3 = [coordinates_x(i,7)-coordinates_x(i,11); coordinates_y(i,7)-coordinates_y(i,11)];
        angle_femur_GEo(i) = rad2deg(acos(dot(vec_GEorigin_to_hip,vec_GEorigin_to_m3)/(norm(vec_GEorigin_to_hip)*norm(vec_GEorigin_to_m3))));
    end
end