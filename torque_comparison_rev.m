clear
% load('results/20240718_HighWalk_init/exp20240718_HighWalk_init_PID_length_and_gain_combination.mat');
% load('results/20240718_HighWalk_init/top10_standing_not_standing_HighWalk_init.mat');
% save_path = 'results/20240718_HighWalk_init/';

load('results/20240712_MuscleLengthTest_PID/exp20240712_MuscleLengthTest_PID_length_and_gain_combination.mat');
load('results/20240712_MuscleLengthTest_PID/top10_standing_not_standing_MuscleLengthTest_PID.mat');
save_path = 'results/20240712_MuscleLengthTest_PID/';

graph_save = true;

exp_num = size(length_and_gain_combination,1);
muscle_lengths = zeros(size(exp_num,1),6); 

for i = 1:exp_num
    L_CFL = length_and_gain_combination(i,1);
    L_Ci = length_and_gain_combination(i,2);
    L_CFLT = length_and_gain_combination(i,3);
    L_GEo = length_and_gain_combination(i,4);
    L_GE = length_and_gain_combination(i,5);
    Pgain = length_and_gain_combination(i,6);
    Igain = length_and_gain_combination(i,7);
    Dgain = length_and_gain_combination(i,8);
    muscle_lengths(i,:) = [i, L_CFL, L_Ci, L_CFLT, L_GEo, L_GE];

    clear t q muscle_tension torque_muscle L_wire k_wire c_wire general_q general_dq l_muscle_list l_link_list data_Q power ankle_lim_index
    % filename = sprintf('results/20240718_HighWalk_init/exp20240718_HighWalk_init_PID_%d_P%d_I%d_D%d_CFL%d_Ci%d_CFLT%d_GEo%d_GE%d.mat',i,Pgain,Igain,Dgain,L_CFL*1000,L_Ci*1000,L_CFLT*1000,L_GEo*1000,L_GE*1000);

    filename = sprintf('results/20240712_MuscleLengthTest_PID/exp20240712_MuscleLengthTest_PID_%d_P%d_I%d_D%d_CFL%d_Ci%d_CFLT%d_GEo%d_GE%d.mat',i,Pgain,Igain,Dgain,L_CFL*1000,L_Ci*1000,L_CFLT*1000,L_GEo*1000,L_GE*1000);
    load(filename);

    muscle_tension = zeros(size(q,1),4);
    torque_muscle = zeros(size(q,1),10);

    for j = 1:size(q,1)
        k_wire = data_k_c_wire(j,2:5).';
        c_wire = data_k_c_wire(j,6:10).';

        general_q = q(j,1:10).';
        general_dq = q(j,11:20).';

        muscle_tension(j,:) = (calc_muscle_tension(l_link_list,l_muscle_list,k_wire,c_wire,general_q, general_dq))';
        torque_muscle(j,:) = (calc_torque_muscle(l_link_list,l_muscle_list,k_wire,c_wire,general_q, general_dq))';
    end

    hip_torque_all(i,:) = torque_muscle(:,6);
    knee_torque_all(i,:) = torque_muscle(:,7);
    ankle_torque_all(i,:) = torque_muscle(:,8);
    hip_angle_all(i,:) = rad2deg(q(:,6));
    knee_angle_all(i,:) = rad2deg(q(:,7));
    ankle_angle_all(i,:) = rad2deg(q(:,8));
    GE_tension_all(i,:) = muscle_tension(:,4);

    % if max(rad2deg(q(:,6))) > -34.9
    %     disp(i);
    % end
end

start_time = 1;
start_num = start_time/0.01;

% start_numからendまでのデータを取得
hip_torque_all = hip_torque_all(:,start_num:end);
knee_torque_all = knee_torque_all(:,start_num:end);
ankle_torque_all = ankle_torque_all(:,start_num:end);
hip_angle_all = hip_angle_all(:,start_num:end);
knee_angle_all = knee_angle_all(:,start_num:end);
ankle_angle_all = ankle_angle_all(:,start_num:end);
GE_tension_all = GE_tension_all(:,start_num:end);

% hip_angle_allの最大値と最小値を取得
hip_angle_max = max(max(hip_angle_all));
hip_angle_min = min(min(hip_angle_all));

% torqueを一次内挿．外挿になる場合はデータを入れない
angle_interp = linspace(hip_angle_min, hip_angle_max, 10000);
for i = 1:exp_num
    hip_torque_interp(i,:) = interp1(hip_angle_all(i,:), hip_torque_all(i,:), angle_interp, 'linear'); %, 'extrap'
    knee_torque_interp(i,:) = interp1(hip_angle_all(i,:), knee_torque_all(i,:), angle_interp, 'linear'); % 'extrap'
    ankle_torque_interp(i,:) = interp1(hip_angle_all(i,:), ankle_torque_all(i,:), angle_interp, 'linear'); % 'extrap'
end

% top10_standing_tension_indexにあるデータのみ取得
hip_torque_all_top10_interp = hip_torque_interp(top10_standing_tension_index,:);
knee_torque_all_top10_interp = knee_torque_interp(top10_standing_tension_index,:);
ankle_torque_all_top10_interp = ankle_torque_interp(top10_standing_tension_index,:);
hip_angle_all_top10_interp = hip_angle_all(top10_standing_tension_index,:);
knee_angle_all_top10_interp = knee_angle_all(top10_standing_tension_index,:);
ankle_angle_all_top10_interp = ankle_angle_all(top10_standing_tension_index,:);
GE_tension_all_top10_interp = GE_tension_all(top10_standing_tension_index,:);
muscle_lengths_top10 = muscle_lengths(top10_standing_tension_index,:);
hip_torque_all_top10_interp_min = min(hip_torque_all_top10_interp);
hip_torque_all_top10_interp_max = max(hip_torque_all_top10_interp);
knee_torque_all_top10_interp_min = min(knee_torque_all_top10_interp);
knee_torque_all_top10_interp_max = max(knee_torque_all_top10_interp);
ankle_torque_all_top10_interp_min = min(ankle_torque_all_top10_interp);
ankle_torque_all_top10_interp_max = max(ankle_torque_all_top10_interp);


% standing_indexのデータのみ取得
hip_torque_all_standing_interp = hip_torque_interp(standing_index,:);
knee_torque_all_standing_interp = knee_torque_interp(standing_index,:);
ankle_torque_all_standing_interp = ankle_torque_interp(standing_index,:);
hip_angle_all_standing_interp = hip_angle_all(standing_index,:);
knee_angle_all_standing_interp = knee_angle_all(standing_index,:);
ankle_angle_all_standing_interp = ankle_angle_all(standing_index,:);
GE_tension_all_standing_interp = GE_tension_all(standing_index,:);
muscle_lengths_standing = muscle_lengths(standing_index,:);
hip_torque_all_standing_interp_min = min(hip_torque_all_standing_interp);
hip_torque_all_standing_interp_max = max(hip_torque_all_standing_interp);
knee_torque_all_standing_interp_min = min(knee_torque_all_standing_interp);
knee_torque_all_standing_interp_max = max(knee_torque_all_standing_interp);
ankle_torque_all_standing_interp_min = min(ankle_torque_all_standing_interp);
ankle_torque_all_standing_interp_max = max(ankle_torque_all_standing_interp);

% not_standing_indexのデータのみ取得
hip_torque_all_not_standing_interp = hip_torque_interp(not_standing_index,:);
knee_torque_all_not_standing_interp = knee_torque_interp(not_standing_index,:);
ankle_torque_all_not_standing_interp = ankle_torque_interp(not_standing_index,:);
hip_angle_all_not_standing_interp = hip_angle_all(not_standing_index,:);
knee_angle_all_not_standing_interp = knee_angle_all(not_standing_index,:);
ankle_angle_all_not_standing_interp = ankle_angle_all(not_standing_index,:);
GE_tension_all_not_standing_interp = GE_tension_all(not_standing_index,:);
muscle_lengths_not_standing = muscle_lengths(not_standing_index,:);
hip_torque_all_not_standing_interp_min = min(hip_torque_all_not_standing_interp);
hip_torque_all_not_standing_interp_max = max(hip_torque_all_not_standing_interp);
knee_torque_all_not_standing_interp_min = min(knee_torque_all_not_standing_interp);
knee_torque_all_not_standing_interp_max = max(knee_torque_all_not_standing_interp);
ankle_torque_all_not_standing_interp_min = min(ankle_torque_all_not_standing_interp);
ankle_torque_all_not_standing_interp_max = max(ankle_torque_all_not_standing_interp);

% hyperextension_index
hip_torque_all_hyperextension_interp = hip_torque_interp(hyperextension_index,:);
knee_torque_all_hyperextension_interp = knee_torque_interp(hyperextension_index,:);
ankle_torque_all_hyperextension_interp = ankle_torque_interp(hyperextension_index,:);
hip_angle_all_hyperextension_interp = hip_angle_all(hyperextension_index,:);
knee_angle_all_hyperextension_interp = knee_angle_all(hyperextension_index,:);
ankle_angle_all_hyperextension_interp = ankle_angle_all(hyperextension_index,:);
GE_tension_all_hyperextension_interp = GE_tension_all(hyperextension_index,:);
muscle_lengths_hyperextension = muscle_lengths(hyperextension_index,:);
hip_torque_all_hyperextension_interp_min = min(hip_torque_all_hyperextension_interp);
hip_torque_all_hyperextension_interp_max = max(hip_torque_all_hyperextension_interp);
knee_torque_all_hyperextension_interp_min = min(knee_torque_all_hyperextension_interp);
knee_torque_all_hyperextension_interp_max = max(knee_torque_all_hyperextension_interp);
ankle_torque_all_hyperextension_interp_min = min(ankle_torque_all_hyperextension_interp);
ankle_torque_all_hyperextension_interp_max = max(ankle_torque_all_hyperextension_interp);


% hip_angle_all_interpに対する，各関節角度におけるtorqueの最大と最小を取得
hip_torque_max = max(hip_torque_interp);
hip_torque_min = min(hip_torque_interp);
knee_torque_max = max(knee_torque_interp);
knee_torque_min = min(knee_torque_interp);
ankle_torque_max = max(ankle_torque_interp);
ankle_torque_min = min(ankle_torque_interp);


% 平均値を取得
hip_torque_all_top10_interp_mean = mean(hip_torque_all_top10_interp);
hip_torque_all_standing_interp_mean = mean(hip_torque_all_standing_interp);
hip_torque_all_not_standing_interp_mean = mean(hip_torque_all_not_standing_interp);
hip_torque_all_hyperextension_interp_mean = mean(hip_torque_all_hyperextension_interp);
knee_torque_all_top10_interp_mean = mean(knee_torque_all_top10_interp);
knee_torque_all_standing_interp_mean = mean(knee_torque_all_standing_interp);
knee_torque_all_not_standing_interp_mean = mean(knee_torque_all_not_standing_interp);
knee_torque_all_hyperextension_interp_mean = mean(knee_torque_all_hyperextension_interp);


green = [0.4660 0.6740 0.1880];
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
purple = [0.4940 0.1840 0.5560];

validMask_top10 = ~isnan(hip_torque_all_top10_interp_min) & ~isnan(hip_torque_all_top10_interp_max);
validMask_standing = ~isnan(hip_torque_all_standing_interp_min) & ~isnan(hip_torque_all_standing_interp_max);
validMask_not_standing = ~isnan(hip_torque_all_not_standing_interp_min) & ~isnan(hip_torque_all_not_standing_interp_max);
validMask_hyperextension = ~isnan(hip_torque_all_hyperextension_interp_min) & ~isnan(hip_torque_all_hyperextension_interp_max);

figure();
hold on;
% fill([angle_interp(validMask_top10), fliplr(angle_interp(validMask_top10))],[hip_torque_all_top10_interp_max(validMask_top10), fliplr(hip_torque_all_top10_interp_min(validMask_top10))], green, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
fill([angle_interp(validMask_standing), fliplr(angle_interp(validMask_standing))],[hip_torque_all_standing_interp_max(validMask_standing), fliplr(hip_torque_all_standing_interp_min(validMask_standing))], blue, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.3);
fill([angle_interp(validMask_not_standing), fliplr(angle_interp(validMask_not_standing))],[hip_torque_all_not_standing_interp_max(validMask_not_standing), fliplr(hip_torque_all_not_standing_interp_min(validMask_not_standing))], red, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.3);
if ~isempty(hyperextension_index)
    fill([angle_interp(validMask_hyperextension), fliplr(angle_interp(validMask_hyperextension))],[hip_torque_all_hyperextension_interp_max(validMask_hyperextension), fliplr(hip_torque_all_hyperextension_interp_min(validMask_hyperextension))], purple, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.3);
end
fill([angle_interp(validMask_top10), fliplr(angle_interp(validMask_top10))],[hip_torque_all_top10_interp_max(validMask_top10), fliplr(hip_torque_all_top10_interp_min(validMask_top10))], [30/255 230/255 30/255], 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
% p1 = plot(angle_interp, hip_torque_all_top10_interp_mean, 'Color', green, 'LineWidth', 2);
p2 = plot(angle_interp, hip_torque_all_standing_interp_mean, 'Color', blue, 'LineWidth', 2);
p3 = plot(angle_interp, hip_torque_all_not_standing_interp_mean, 'Color', red, 'LineWidth', 2);
if ~isempty(hyperextension_index)
    p4 = plot(angle_interp, hip_torque_all_hyperextension_interp_mean, 'Color', purple, 'LineWidth', 2);
end
p1 = plot(angle_interp, hip_torque_all_top10_interp_mean, 'Color', [30/255 170/255 70/255], 'LineWidth', 2);
hold off;
set(gca, 'XDir', 'reverse');
xlim([-100, -10]);
ylim([-40, 0]);
h_axes = gca;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;
xlabel('Hip angle [deg]', 'FontSize', 25);
ylabel('Hip torque [Nm]','FontSize', 25);
if ~isempty(hyperextension_index)
    legend([p1 p2 p3 p4], '(i)', '(ii)', '(iii)', '(iv)', 'FontSize', 25, 'Location', 'best');
else
    legend([p1 p2 p3], '(i)', '(ii)', '(iii)', 'FontSize', 25, 'Location', 'best');
end
if graph_save == true
    exportgraphics(gca, [save_path 'hip_torque_comparison_rev.pdf'], 'ContentType', 'vector');
    saveas(gcf, [save_path 'hip_torque_comparison_rev.png']);
end


validMask_top10 = ~isnan(knee_torque_all_top10_interp_min) & ~isnan(knee_torque_all_top10_interp_max);
validMask_standing = ~isnan(knee_torque_all_standing_interp_min) & ~isnan(knee_torque_all_standing_interp_max);
validMask_not_standing = ~isnan(knee_torque_all_not_standing_interp_min) & ~isnan(knee_torque_all_not_standing_interp_max);
validMask_hyperextension = ~isnan(knee_torque_all_hyperextension_interp_min) & ~isnan(knee_torque_all_hyperextension_interp_max);


figure();
hold on;
% fill([angle_interp(validMask_top10), fliplr(angle_interp(validMask_top10))],[knee_torque_all_top10_interp_max(validMask_top10), fliplr(knee_torque_all_top10_interp_min(validMask_top10))], green, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
fill([angle_interp(validMask_standing), fliplr(angle_interp(validMask_standing))],[knee_torque_all_standing_interp_max(validMask_standing), fliplr(knee_torque_all_standing_interp_min(validMask_standing))], blue, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.3);
fill([angle_interp(validMask_not_standing), fliplr(angle_interp(validMask_not_standing))],[knee_torque_all_not_standing_interp_max(validMask_not_standing), fliplr(knee_torque_all_not_standing_interp_min(validMask_not_standing))], red, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.3);
if ~isempty(hyperextension_index)
    fill([angle_interp(validMask_hyperextension), fliplr(angle_interp(validMask_hyperextension))],[knee_torque_all_hyperextension_interp_max(validMask_hyperextension), fliplr(knee_torque_all_hyperextension_interp_min(validMask_hyperextension))], purple, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.3);
end
fill([angle_interp(validMask_top10), fliplr(angle_interp(validMask_top10))],[knee_torque_all_top10_interp_max(validMask_top10), fliplr(knee_torque_all_top10_interp_min(validMask_top10))], [30/255 230/255 30/255], 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
% p1 = plot(angle_interp, knee_torque_all_top10_interp_mean, 'Color', green, 'LineWidth', 2);
p2 = plot(angle_interp, knee_torque_all_standing_interp_mean, 'Color', blue, 'LineWidth', 2);
p3 = plot(angle_interp, knee_torque_all_not_standing_interp_mean, 'Color', red, 'LineWidth', 2);
if ~isempty(hyperextension_index)
    p4 = plot(angle_interp, knee_torque_all_hyperextension_interp_mean, 'Color', purple, 'LineWidth', 2);
end
p1 = plot(angle_interp, knee_torque_all_top10_interp_mean, 'Color', [30/255 170/255 70/255], 'LineWidth', 2);
hold off;
set(gca, 'XDir', 'reverse');
xlim([-100, -10]);
ylim([-18, 0]);
h_axes = gca;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;
xlabel('Hip angle [deg]','FontSize', 25);
ylabel('Knee torque [Nm]','FontSize', 25);
if ~isempty(hyperextension_index)
    legend([p1 p2 p3 p4], '(i)', '(ii)', '(iii)', '(iv)', 'FontSize', 25, 'Location', 'best');
else
    legend([p1 p2 p3], '(i)', '(ii)', '(iii)', 'FontSize', 25, 'Location', 'best');
end
if graph_save == true
    exportgraphics(gca, [save_path 'knee_torque_comparison_rev.pdf'], 'ContentType', 'vector');
    saveas(gcf, [save_path 'knee_torque_comparison_rev.png']);
end



validMask_top10 = ~isnan(knee_torque_all_top10_interp_min) & ~isnan(knee_torque_all_top10_interp_max) & ~isnan(hip_torque_all_top10_interp_min) & ~isnan(hip_torque_all_top10_interp_max);
validMask_standing = ~isnan(knee_torque_all_standing_interp_min) & ~isnan(knee_torque_all_standing_interp_max) & ~isnan(hip_torque_all_standing_interp_min) & ~isnan(hip_torque_all_standing_interp_max);
validMask_not_standing = ~isnan(knee_torque_all_not_standing_interp_min) & ~isnan(knee_torque_all_not_standing_interp_max) & ~isnan(hip_torque_all_not_standing_interp_min) & ~isnan(hip_torque_all_not_standing_interp_max);
validMask_hyperextension = ~isnan(knee_torque_all_hyperextension_interp_min) & ~isnan(knee_torque_all_hyperextension_interp_max) & ~isnan(hip_torque_all_hyperextension_interp_min) & ~isnan(hip_torque_all_hyperextension_interp_max);

knee_hip_ratio_top10 = knee_torque_all_top10_interp./hip_torque_all_top10_interp;
knee_hip_ratio_standing = knee_torque_all_standing_interp./hip_torque_all_standing_interp;
knee_hip_ratio_not_standing = knee_torque_all_not_standing_interp./hip_torque_all_not_standing_interp;
knee_hip_ratio_hyperextension = knee_torque_all_hyperextension_interp./hip_torque_all_hyperextension_interp;

% 最大値と最小値を取得
knee_hip_ratio_top10_max = max(knee_hip_ratio_top10);
knee_hip_ratio_top10_min = min(knee_hip_ratio_top10);
knee_hip_ratio_standing_max = max(knee_hip_ratio_standing);
knee_hip_ratio_standing_min = min(knee_hip_ratio_standing);
knee_hip_ratio_not_standing_max = max(knee_hip_ratio_not_standing);
knee_hip_ratio_not_standing_min = min(knee_hip_ratio_not_standing);
knee_hip_ratio_hyperextension_max = max(knee_hip_ratio_hyperextension);
knee_hip_ratio_hyperextension_min = min(knee_hip_ratio_hyperextension);


figure();
hold on;
% fill([angle_interp(validMask_top10), fliplr(angle_interp(validMask_top10))],[knee_hip_ratio_top10_max(validMask_top10), fliplr(knee_hip_ratio_top10_min(validMask_top10))], green, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
fill([angle_interp(validMask_standing), fliplr(angle_interp(validMask_standing))],[knee_hip_ratio_standing_max(validMask_standing), fliplr(knee_hip_ratio_standing_min(validMask_standing))], blue, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.3);
fill([angle_interp(validMask_not_standing), fliplr(angle_interp(validMask_not_standing))],[knee_hip_ratio_not_standing_max(validMask_not_standing), fliplr(knee_hip_ratio_not_standing_min(validMask_not_standing))], red, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.3);
if ~isempty(hyperextension_index)
    fill([angle_interp(validMask_hyperextension), fliplr(angle_interp(validMask_hyperextension))],[knee_hip_ratio_hyperextension_max(validMask_hyperextension), fliplr(knee_hip_ratio_hyperextension_min(validMask_hyperextension))], purple, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.3);
end
fill([angle_interp(validMask_top10), fliplr(angle_interp(validMask_top10))],[knee_hip_ratio_top10_max(validMask_top10), fliplr(knee_hip_ratio_top10_min(validMask_top10))], [30/255 230/255 30/255], 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
% p1 = plot(angle_interp, knee_torque_all_top10_interp_mean./hip_torque_all_top10_interp_mean, 'Color', green, 'LineWidth', 2);
p2 = plot(angle_interp, knee_torque_all_standing_interp_mean./hip_torque_all_standing_interp_mean, 'Color', blue, 'LineWidth', 2);
p3 = plot(angle_interp, knee_torque_all_not_standing_interp_mean./hip_torque_all_not_standing_interp_mean, 'Color', red, 'LineWidth', 2);
if ~isempty(hyperextension_index)
    p4 = plot(angle_interp, knee_torque_all_hyperextension_interp_mean./hip_torque_all_hyperextension_interp_mean, 'Color', purple, 'LineWidth', 2);
end
p1 = plot(angle_interp, knee_torque_all_top10_interp_mean./hip_torque_all_top10_interp_mean, 'Color', [30/255 170/255 70/255], 'LineWidth', 2);
hold off;
set(gca, 'XDir', 'reverse');
xlim([-100, -10]);
ylim([0, 0.7]);
h_axes = gca;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;
xlabel('Hip angle [deg]','FontSize', 25);
ylabel('Knee/Hip torque ratio [-]', 'FontSize', 25);
if ~isempty(hyperextension_index)
    legend([p1 p2 p3 p4], '(i)', '(ii)', '(iii)', '(iv)', 'FontSize', 25, 'Location', 'best');
else
    legend([p1 p2 p3], '(i)', '(ii)', '(iii)', 'FontSize', 25, 'Location', 'best');
end
if graph_save == true
    exportgraphics(gca, [save_path 'knee_hip_torque_ratio_comparison_rev.pdf'], 'ContentType', 'vector');
    saveas(gcf, [save_path 'knee_hip_torque_ratio_comparison_rev.png']);
end