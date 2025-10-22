clear

load('results/20240726_init_condition_test_MinWork_per2mm/exp20240726_init_conditions.mat');
load('results/20240726_init_condition_test_MinWork_per2mm/exp20240726_init_condition_test_MinWork_per2mm_length_and_gain_combination.mat');
save_path = 'results/20240726_init_condition_test_MinWork_per2mm/';

% load('results/20240726_init_condition_test_per2mm/exp20240726_init_conditions.mat');
% load('results/20240726_init_condition_test_per2mm/exp20240726_init_condition_test_per2mm_length_and_gain_combination.mat');
% save_path = 'results/20240726_init_condition_test_per2mm/';

MinWork_flag = false;
if contains(save_path, '20240726_init_condition_test_MinWork_per2mm')
    MinWork_flag = true;
end

L_CFL = length_and_gain_combination(1,1);
L_Ci = length_and_gain_combination(1,2);
L_CFLT = length_and_gain_combination(1,3);
L_GEo = length_and_gain_combination(1,4);
L_GE = length_and_gain_combination(1,5);

init_theta2_ = 10:2:70;
init_theta2_ = -init_theta2_/180*pi;
init_theta3_ = 30:2:110;
init_theta3_ = -init_theta3_/180*pi;
[init_theta2_list, init_theta3_list] = ndgrid(init_theta2_, init_theta3_);
init_condition_list_for_lim = [init_theta2_list(:), init_theta3_list(:)];
r = 0.015;
L_met = 90.739/1000;
limit_ankle_ex = 120/180*pi;

evaluate = []; % 1列目：theta2, 2列目：theta3, 3列目：work (立てなかった場合は-1, 除外された場合は-2)

% Evaluate initial conditions
for i = size(init_condition_list_for_lim,1):-1:1
    theta2_ = init_condition_list_for_lim(i,1);
    theta3_ = init_condition_list_for_lim(i,2);
    theta4_ = 90/180*pi-(90/180*pi+theta2_+theta3_)-(asin(r/L_met)-0.03);
    if theta4_ > limit_ankle_ex %初期の足関節角度がlimit_ankle_exよりも大きくなる時は除外
        evaluate(end+1,1) = rad2deg(theta2_);
        evaluate(end,2) = rad2deg(theta3_);
        evaluate(end,3) = -2;
    end
end

for i = 1:size(init_condition_list,1)
    theta2 = init_condition_list(i,1); %hip
    theta3 = init_condition_list(i,2); %knee
    
    clear t q muscle_tension torque_muscle L_wire k_wire c_wire general_q general_dq l_muscle_list l_link_list data_Q power ankle_lim_index
    filename = sprintf('results/20240726_init_condition_test_MinWork_per2mm/exp20240726_init_condition_test_MinWork_per2mm_%d_Hip%d_Knee%d_CFL%d_Ci%d_CFLT%d_GEo%d_GE%d.mat', i, -rad2deg(theta2), -rad2deg(theta3), L_CFL*1000, L_Ci*1000, L_CFLT*1000, L_GEo*1000, L_GE*1000);
    if MinWork_flag == true && i > 699
        filename = sprintf('results/20240726_init_condition_test_MinWork_per2mm/exp20240726_init_condition_test_per2mm_MinWork_%d_Hip%d_Knee%d_CFL%d_Ci%d_CFLT%d_GEo%d_GE%d.mat', i, -rad2deg(theta2), -rad2deg(theta3), L_CFL*1000, L_Ci*1000, L_CFLT*1000, L_GEo*1000, L_GE*1000);
    end
    load(filename);

    % t = 10s のときのy座標
    y_hip = q(1200,2)-l_link_list(6)*cos(q(1200,5));

    if y_hip > 0
        power = q(50:end,20).*data_Q(50:end,11);
        work= trapz(t(50:end,1),power);
        evaluate(end+1,1) = rad2deg(theta2);
        evaluate(end,2) = rad2deg(theta3);
        evaluate(end,3) = work;
    elseif y_hip <= 0
        evaluate(end+1,1) = rad2deg(theta2);
        evaluate(end,2) = rad2deg(theta3);
        evaluate(end,3) = -1;
    end
end

% グリッドデータの準備
theta2_vals = unique(evaluate(:,1));
theta3_vals = unique(evaluate(:,2));
work_grid = nan(length(theta2_vals), length(theta3_vals));

for i = 1:size(evaluate,1)
    row = find(theta2_vals == evaluate(i,1));
    col = find(theta3_vals == evaluate(i,2));
    work_grid(row, col) = evaluate(i,3);
end

% メッシュデータの作成
[X, Y] = meshgrid(theta2_vals, theta3_vals);
Z = work_grid';

% プロット
figure
hold on
for i = 1:size(Z, 1)-1
    for j = 1:size(Z, 2)-1
        if ~isnan(Z(i, j))
            if Z(i, j) == -1
                fill3([X(i,j) X(i,j+1) X(i+1,j+1) X(i+1,j)], [Y(i,j) Y(i,j+1) Y(i+1,j+1) Y(i+1,j)], [0 0 0 0], 'r');
            elseif Z(i, j) == -2
                fill3([X(i,j) X(i,j+1) X(i+1,j+1) X(i+1,j)], [Y(i,j) Y(i,j+1) Y(i+1,j+1) Y(i+1,j)], [0 0 0 0], [0.5 0.5 0.5]);
            else
                fill3([X(i,j) X(i,j+1) X(i+1,j+1) X(i+1,j)], [Y(i,j) Y(i,j+1) Y(i+1,j+1) Y(i+1,j)], [0 0 0 0], Z(i, j));
            end
        end
    end
end
colormap("parula")
c = colorbar;
clim([0 7]);
c.Label.String = 'Work [J]'; %フォントサイズも次の行で変更
c.Label.FontSize = 25;
xlabel('Initial hip angle [deg]', 'FontSize', 25)
ylabel('Initial knee angle [deg]', 'FontSize', 25)
xlim([min(theta2_vals) max(theta2_vals)])
ylim([min(theta3_vals) max(theta3_vals)])
h_axes = gca;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;
view(2)
exportgraphics(gca, [save_path 'stand_condition_plot_colormap.pdf'], 'ContentType', 'vector');




% clear

% load('results/20240726_init_condition_test_MinWork_per2mm/exp20240723_HighWalk_init_condition.mat');
% load('results/20240726_init_condition_test_MinWork_per2mm/exp20240726_init_condition_test_MinWork_per2mm_length_and_gain_combination.mat');
% save_path = 'results/20240726_init_condition_test_MinWork_per2mm/';

% L_CFL = length_and_gain_combination(1,1);
% L_Ci = length_and_gain_combination(1,2);
% L_CFLT = length_and_gain_combination(1,3);
% L_GEo = length_and_gain_combination(1,4);
% L_GE = length_and_gain_combination(1,5);

% init_theta2_ = 10:3:70;
% init_theta2_ = -init_theta2_/180*pi;
% init_theta3_ = 40:3:100;
% init_theta3_ = -init_theta3_/180*pi;
% [init_theta2_list, init_theta3_list] = ndgrid(init_theta2_, init_theta3_);
% init_condition_list_for_lim = [init_theta2_list(:), init_theta3_list(:)];
% r = 0.015;
% L_met = 90.739/1000;
% limit_ankle_ex = 120/180*pi;

% evaluate = []; % 1列目：theta2, 2列目：theta3, 3列目：work (立てなかった場合は-1, 除外された場合は-2)

% % Evaluate initial conditions
% for i = size(init_condition_list_for_lim,1):-1:1
%     theta2_ = init_condition_list_for_lim(i,1);
%     theta3_ = init_condition_list_for_lim(i,2);
%     theta4_ = 90/180*pi-(90/180*pi+theta2_+theta3_)-(asin(r/L_met)-0.03);
%     if theta4_ > limit_ankle_ex %初期の足関節角度がlimit_ankle_exよりも大きくなる時は除外
%         evaluate(end+1,1) = rad2deg(theta2_);
%         evaluate(end,2) = rad2deg(theta3_);
%         evaluate(end,3) = -2;
%     end
% end

% for i = 1:size(init_condition_list,1)
%     theta2 = init_condition_list(i,1); %hip
%     theta3 = init_condition_list(i,2); %knee
    
%     clear t q muscle_tension torque_muscle L_wire k_wire c_wire general_q general_dq l_muscle_list l_link_list data_Q power ankle_lim_index
%     filename = sprintf('results/20240726_init_condition_test_MinWork_per2mm/exp20240726_init_condition_test_MinWork_per2mm_%d_Hip%d_Knee%d_CFL%d_Ci%d_CFLT%d_GEo%d_GE%d.mat', i, -rad2deg(theta2), -rad2deg(theta3), L_CFL*1000, L_Ci*1000, L_CFLT*1000, L_GEo*1000, L_GE*1000);
%     load(filename);

%     % t = 10s のときのy座標
%     y_hip = q(1000,2)-l_link_list(6)*cos(q(1000,5));

%     if y_hip > 0
%         power = q(50:end,20).*data_Q(50:end,11);
%         work= trapz(t(50:end,1),power);
%         evaluate(end+1,1) = rad2deg(theta2);
%         evaluate(end,2) = rad2deg(theta3);
%         evaluate(end,3) = work;
%     elseif y_hip <= 0
%         evaluate(end+1,1) = rad2deg(theta2);
%         evaluate(end,2) = rad2deg(theta3);
%         evaluate(end,3) = -1;
%     end
% end

% % グリッドデータの準備
% theta2_vals = unique(evaluate(:,1));
% theta3_vals = unique(evaluate(:,2));
% work_grid = nan(length(theta2_vals), length(theta3_vals));

% for i = 1:size(evaluate,1)
%     row = find(theta2_vals == evaluate(i,1));
%     col = find(theta3_vals == evaluate(i,2));
%     work_grid(row, col) = evaluate(i,3);
% end

% % メッシュデータの作成
% [X, Y] = meshgrid(theta2_vals, theta3_vals);
% Z = work_grid';

% % プロット
% figure
% hold on
% unable_to_stand_patch = [];
% exclude_patch = [];
% for i = 1:size(Z, 1)-1
%     for j = 1:size(Z, 2)-1
%         if ~isnan(Z(i, j))
%             if Z(i, j) == -1
%                 h = fill3([X(i,j) X(i,j+1) X(i+1,j+1) X(i+1,j)], [Y(i,j) Y(i,j+1) Y(i+1,j+1) Y(i+1,j)], [0 0 0 0], 'r');
%                 unable_to_stand_patch = [unable_to_stand_patch, h];
%             elseif Z(i, j) == -2
%                 h = fill3([X(i,j) X(i,j+1) X(i+1,j+1) X(i+1,j)], [Y(i,j) Y(i,j+1) Y(i+1,j+1) Y(i+1,j)], [0 0 0 0], [0.5 0.5 0.5]);
%                 exclude_patch = [exclude_patch, h];
%             else
%                 fill3([X(i,j) X(i,j+1) X(i+1,j+1) X(i+1,j)], [Y(i,j) Y(i,j+1) Y(i+1,j+1) Y(i+1,j)], [0 0 0 0], Z(i, j));
%             end
%         end
%     end
% end
% colormap("parula")
% c = colorbar;
% clim([0 7]);
% c.Label.String = 'Work [J]';

% % 凡例を追加
% legend([unable_to_stand_patch(1), exclude_patch(1)], {'Unable to stand', 'Exclude from initial value'}, 'Location', 'northoutside');
% lgd.Position = [0.58 0.93 0.3 0.06];

% xlabel('Hip angle [deg]')
% ylabel('Knee angle [deg]')
% axis square
% xlim([min(theta2_vals) max(theta2_vals)])
% ylim([min(theta3_vals) max(theta3_vals)])
% view(2)
% exportgraphics(gca, [save_path 'stand_condition_plot_colormap.pdf'], 'ContentType', 'vector');