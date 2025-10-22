clear all
% load('results/20240718_HighWalk_init/exp20240718_HighWalk_init_PID_length_and_gain_combination.mat');
% save_path = 'results/20240718_HighWalk_init/';
load('results/20240712_MuscleLengthTest_PID/exp20240712_MuscleLengthTest_PID_length_and_gain_combination.mat');
save_path = 'results/20240712_MuscleLengthTest_PID/';

graph_view = false;
graph_save = false;

exp_num = size(length_and_gain_combination,1);
% exp_num = 300;

max_tension_CFL = zeros(size(exp_num,1),3); %index, tension, stand or not
max_torque = zeros(size(exp_num,1),5); %index, hip, knee, ankle, stand or not
torque_hip_knee_ratio = zeros(size(exp_num,1),2); %index, ratio
muscle_lengths = zeros(size(exp_num,1),6); %index, L_CFL, L_Ci, L_CFLT, L_GEo, L_GE
y_hip_coordination = zeros(size(exp_num,1),2); %index, y_hip
work = zeros(size(exp_num,1),3); %index, evaluate_value, stand or not

ankle_lim_time = zeros(exp_num,1);
ankle_lim_flag = false;
torque_hip_ankle_lim = zeros(exp_num,2);
torque_knee_ankle_lim = zeros(exp_num,2);

knee_angle_end = zeros(size(exp_num,1),3); %index, angle, stand or not

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

        if j == 1000
            y_hip_coordination(i,2) = general_q(2)-l_link_list(6)*cos(general_q(5));
        end
    end

    power = q(50:end,20).*data_Q(50:end,11);
    work(i,1) = i;
    work(i,2) = trapz(t(50:end,1),power);

    max_tension_CFL(i,1) = i;
    max_tension_CFL(i,2) = max(data_Q(50:end,11));
    max_torque(i,1) = i;
    max_torque(i,2) = min(torque_muscle(100:1000,6));
    %負の値で出てくるので，最小値を取る
    max_torque(i,3) = min(torque_muscle(100:1000,7));
    max_torque(i,4) = min(torque_muscle(100:1000,8));
    torque_hip_knee_ratio(i,1) = i;
    torque_hip_knee_ratio(i,2) = max_torque(i,3)/max_torque(i,2);

    knee_angle_end(i,1) = i;
    knee_angle_end(i,2) = q(end,7);

    if y_hip_coordination(i,2) > 0.0
        max_tension_CFL(i,3) = 1;
        max_torque(i,5) = 1;
        work(i,3) = 1;
        knee_angle_end(i,3) = 1;
    else
        max_tension_CFL(i,3) = 0;
        max_torque(i,5) = 0;
        work(i,3) = 0;
        if knee_angle_end(i,2) < 0
            knee_angle_end(i,3) = 2; %膝が曲がっている
        else
            knee_angle_end(i,3) = 3; %膝過伸展
        end
    end

    for j = 1:size(q,1)
        if ankle_lim_flag == false && q(j,8) > 120/180*pi
            ankle_lim_time(i) = t(j,1);
            % torque_hip_ankle_lim(i,1) = i;
            % torque_hip_ankle_lim(i,2) = -min(torque_muscle(j:end,6));
            torque_knee_ankle_lim(i,1) = i;
            torque_knee_ankle_lim(i,2) = max(-torque_muscle(j:end,7));
            torque_hip_ankle_lim(i,1) = i;
            if torque_knee_ankle_lim(i,2) == 0
                torque_hip_ankle_lim(i,2) = -torque_muscle(j,6);
                ankle_lim_flag = true;
            else
                % torque_knee_ankle_lim(i,2)となるindexを取得
                ankle_lim_index = find(-torque_muscle(j:end,7) == torque_knee_ankle_lim(i,2));
                torque_hip_ankle_lim(i,2) = -torque_muscle(ankle_lim_index+j-1,6);
                ankle_lim_flag = true;
            end
        end
    end
    ankle_lim_flag = false;
end

standing_index = find(work(:,3) == 1);
not_standing_index = find(work(:,3) == 0);

% standing_indexの中から仕事が最も小さい10個を取得
[~, top10_tension_index] = sort(work(standing_index, 2));
top10_standing_tension_index = standing_index(top10_tension_index(1:min(10, length(top10_tension_index))));

% 仕事が最も大きい10個を取得
[~, top10_work_index] = sort(work(standing_index, 2), 'descend');
top10_standing_work_index = standing_index(top10_work_index(1:min(10, length(top10_work_index))));



%切り抜く上位の数を指定
top_num = 10;
[~, top_tension_index] = sort(work(standing_index, 2));
top_standing_tension_index = standing_index(top_tension_index(1:min(top_num, length(top_tension_index))));

% standing_indexの中からtop_standing_tension_indexに含まれるものを除外
standing_index= standing_index(~ismember(standing_index, top_standing_tension_index));


exp_num = size(muscle_lengths, 1);
coorddata = [1 2 3 4 5];
% plot_tableの行を定義
plot_table = zeros(exp_num + 2, 6);  % +2は最大値と最小値用
plot_table(3:end, 1) = muscle_lengths(:,3);
plot_table(3:end, 2) = muscle_lengths(:,4);
plot_table(3:end, 3) = muscle_lengths(:,5);
plot_table(3:end, 4) = muscle_lengths(:,6);
plot_table(3:end, 5) = work(:,2);
% グループ分けのラベル付け (6列目)
plot_table(3:end, 6) = 0;  % 他のデータを除外対象に設定
plot_table(top_standing_tension_index + 2, 6) = 1; 
% 最大値と最小値を最初の2行に追加
plot_table(1, 1:5) = [max(muscle_lengths(:,3)), max(muscle_lengths(:,4)), max(muscle_lengths(:,5)), max(muscle_lengths(:,6)), 5.2];
plot_table(2, 1:5) = [min(muscle_lengths(:,3)), min(muscle_lengths(:,4)), min(muscle_lengths(:,5)), min(muscle_lengths(:,6)), 0];
plot_table(1:2, 6) = 2;  % 新しいグループID（最大・最小用）
% 'top30'と最大最小のみをフィルタリング
filtered_table = plot_table(plot_table(:,6) > 0, :);
% parallelplotの作成

% figure();
% p = parallelplot(filtered_table, 'CoordinateData', coorddata, 'GroupData', filtered_table(:,6), 'LineWidth', 1.5);
% p.CoordinateTickLabels = {'Ci [m]', 'CFLT [m]', 'GEo [m]', 'GE [m]', 'Work [J]'};
% % 各グループの色と透明度設定
% p.Color = [1 1 1; 0.4660 0.6740 0.1880];  % 'top30'と'最大最小' の色設定
% % 図の保存
% if graph_save == true
%     filename = ['top' num2str(top_num) '_min_max_parallelcoords.pdf'];
%     exportgraphics(gcf, [save_path, filename], 'ContentType', 'vector');
% end


% 条件分けされたtable作成
exp_num = size(muscle_lengths, 1);
coorddata = [1 2 3 4 5];
% plot_tableの行を定義
plot_table = zeros(exp_num + 2, 6);  % +2は最大値と最小値用
plot_table(3:end, 1) = muscle_lengths(:,3);
plot_table(3:end, 2) = muscle_lengths(:,4);
plot_table(3:end, 3) = muscle_lengths(:,5);
plot_table(3:end, 4) = muscle_lengths(:,6);
plot_table(3:end, 5) = work(:,2);
% グループ分けのラベル付け (6列目)
plot_table(:,6) = 0;  % デフォルトで除外対象
plot_table(not_standing_index + 2, 6) = 1;  % 'not standing'
plot_table(standing_index + 2, 6) = 2;      % 'standing'
% 最大値と最小値を最初の2行に追加
plot_table(1, 1:5) = [max(muscle_lengths(:,3)), max(muscle_lengths(:,4)), max(muscle_lengths(:,5)), max(muscle_lengths(:,6)), 5.2];
plot_table(2, 1:5) = [min(muscle_lengths(:,3)), min(muscle_lengths(:,4)), min(muscle_lengths(:,5)), min(muscle_lengths(:,6)), 0];
plot_table(1:2, 6) = 3;  % 新しいグループID（最大・最小用）
% 'not standing'のみのデータをフィルタリング
not_standing_table = plot_table(plot_table(:,6) == 1 | plot_table(:,6) == 3, :);
% 'standing'のみのデータをフィルタリング
% standing_table = plot_table(plot_table(:,6) == 2 | plot_table(:,6) == 3, :);

standing_table = plot_table(plot_table(:,6) == 2, :);

% 'not standing'のみのparallelplot
% figure();
% p1 = parallelplot(not_standing_table(:,1:5), 'CoordinateData', coorddata, 'GroupData', not_standing_table(:,6));
% p1.CoordinateTickLabels = {'Ci [m]', 'CFLT [m]', 'GEo [m]', 'GE [m]', 'Work [J]'};
% % 色の設定： 'not standing'は赤、最大・最小は白
% p1.Color = [1 1 1; 0.6350 0.0780 0.1840];  % 赤色と白色
% % 'not standing'グラフの保存
% if graph_save == true
%     exportgraphics(gcf, [save_path, 'not_standing_with_limits_parallelcoords.pdf'], 'ContentType', 'vector');
% end


% % 'standing'のみのparallelplot
% figure();
% p2 = parallelplot(standing_table(:,1:5), 'CoordinateData', coorddata, 'GroupData', standing_table(:,6));
% p2.CoordinateTickLabels = {'Ci [m]', 'CFLT [m]', 'GEo [m]', 'GE [m]', 'Work [J]'};
% % 色の設定： 'standing'は柔らかい緑、最大・最小は白
% p2.Color = [1 1 1; 0 0.4470 0.7410];  % 緑色と白色
% % 'standing'グラフの保存
% if graph_save == true
%     exportgraphics(gcf, [save_path, 'standing_with_limits_parallelcoords.pdf'], 'ContentType', 'vector');
% end




% カラーマップの準備
% プロット
figure();
cmap = colormap('turbo');
hold on; % 複数のプロットを重ねる

x = [1, 2, 3, 4];

Ci_list = [0.036, 0.038, 0.04, 0.042, 0.044];
Ci_list_norm = (Ci_list - min(Ci_list)) / (max(Ci_list) - min(Ci_list));
CFLT_list = [0.095, 0.098, 0.101, 0.104, 0.107];
CFLT_list_norm = (CFLT_list - min(CFLT_list)) / (max(CFLT_list) - min(CFLT_list));
GEo_list = [0.035, 0.037, 0.039, 0.041, 0.043];
GEo_list_norm = (GEo_list - min(GEo_list)) / (max(GEo_list) - min(GEo_list));
GE_list = [0.176, 0.179, 0.182, 0.185, 0.188, 0.191, 0.194, 0.197, 0.200];
GE_list_norm = (GE_list - min(GE_list)) / (max(GE_list) - min(GE_list));

new_standing_table = zeros(size(standing_table, 1), 5);
% standing_tableのmuscle_lengthsの値を0~1に正規化
for i = 1:size(standing_table, 1)
    new_standing_table(i, 1) = (standing_table(i, 1) - min(Ci_list)) / (max(Ci_list) - min(Ci_list));
    new_standing_table(i, 2) = (standing_table(i, 2) - min(CFLT_list)) / (max(CFLT_list) - min(CFLT_list));
    new_standing_table(i, 3) = (standing_table(i, 3) - min(GEo_list)) / (max(GEo_list) - min(GEo_list));
    new_standing_table(i, 4) = (standing_table(i, 4) - min(GE_list)) / (max(GE_list) - min(GE_list));
    new_standing_table(i, 5) = standing_table(i, 5);
end
% new_standing_tableを5列目の値で昇順にソート
new_standing_table = sortrows(new_standing_table, 5);

y_offset_scale = 0.0002; % y軸方向のオフセットの幅
y_offset_values = zeros(size(standing_table, 1), 4);

Ci_count = [1, 1, 1, 1, 1];
CFLT_count = [1, 1, 1, 1, 1];
GEo_count = [1, 1, 1, 1, 1];
GE_count = [1, 1, 1, 1, 1, 1, 1, 1, 1];
% offsetの設定
for i = 1:size(standing_table, 1)
    for j = 1:5
        if new_standing_table(i, 1) == Ci_list_norm(j)
            y_offset_values(i, 1) = y_offset_scale * Ci_count(j);
            Ci_count(j) = Ci_count(j) + 1;
        end
        if new_standing_table(i, 2) == CFLT_list_norm(j)
            y_offset_values(i, 2) = y_offset_scale * CFLT_count(j);
            CFLT_count(j) = CFLT_count(j) + 1;
        end
        if new_standing_table(i, 3) == GEo_list_norm(j)
            y_offset_values(i, 3) = y_offset_scale * GEo_count(j);
            GEo_count(j) = GEo_count(j) + 1;
        end
    end
    for j = 1:9
        if new_standing_table(i, 4) == GE_list_norm(j)
            y_offset_values(i, 4) = y_offset_scale * GE_count(j);
            GE_count(j) = GE_count(j) + 1;
        end
    end   
end

% standing_tableについて各データをプロット
for i = 1:size(standing_table, 1)

    y_values = new_standing_table(i, 1:4) + y_offset_values(i, 1:4);
    
    % 仕事に応じた色を取得
    val_scaled = (new_standing_table(i, 5) - min(new_standing_table(1:end, 5))) / (max(new_standing_table(1:end, 5)) - min(new_standing_table(1:end, 5)));
    color_index = round(val_scaled * (size(cmap, 1) - 1)) + 1;
    color = cmap(color_index, :);
    % 線をプロット
    plot(x, y_values, '-', 'Color', color);
    
    % マーカーをプロット
    scatter(x, y_values, 10, color, 'filled');
end

caxis([min(new_standing_table(1:end, 5)), max(new_standing_table(1:end, 5))]);
c = colorbar;
c.Label.String = 'Work [J]';
c.Label.FontSize = 25;
colormap(cmap); % カラーバーにカラーマップを適用
ylim([-0.1, 1.1]);
hold off;
xticks([1, 2, 3, 4]);
xticklabels({'Ci [m]', 'CFLT[m]', 'GEo[m]', 'GE[m]'});
yticks([]);  % y軸のメモリ（目盛り）を削除
yticklabels([]); % y軸のラベルを削除
exportgraphics(gcf, [save_path, 'standing_parallelcoords_color.pdf'], 'ContentType', 'vector');


dev_3 =floor(size(standing_table, 1)/3);
cmap = colormap('turbo');
figure();
hold on
Ci_count = [1, 1, 1, 1, 1];
CFLT_count = [1, 1, 1, 1, 1];
GEo_count = [1, 1, 1, 1, 1];
GE_count = [1, 1, 1, 1, 1, 1, 1, 1, 1];
% offsetの設定
for i = 1:dev_3
    for j = 1:5
        if new_standing_table(i, 1) == Ci_list_norm(j)
            y_offset_values(i, 1) = y_offset_scale * Ci_count(j);
            Ci_count(j) = Ci_count(j) + 1;
        end
        if new_standing_table(i, 2) == CFLT_list_norm(j)
            y_offset_values(i, 2) = y_offset_scale * CFLT_count(j);
            CFLT_count(j) = CFLT_count(j) + 1;
        end
        if new_standing_table(i, 3) == GEo_list_norm(j)
            y_offset_values(i, 3) = y_offset_scale * GEo_count(j);
            GEo_count(j) = GEo_count(j) + 1;
        end
    end
    for j = 1:9
        if new_standing_table(i, 4) == GE_list_norm(j)
            y_offset_values(i, 4) = y_offset_scale * GE_count(j);
            GE_count(j) = GE_count(j) + 1;
        end
    end   
end
for i = 1:dev_3

    y_values = new_standing_table(i, 1:4) + y_offset_values(i, 1:4)*3;
    
    % 仕事に応じた色を取得
    % val_scaled = (new_standing_table(i, 5) - min(new_standing_table(1:end, 5))) / (max(new_standing_table(1:end, 5)) - min(new_standing_table(1:end, 5)));
    val_scaled = (new_standing_table(i, 5) - min(new_standing_table(1:end, 5))) / (max(new_standing_table(1:end, 5)) - min(new_standing_table(1:end, 5)));
    color_index = round(val_scaled * (size(cmap, 1) - 1)) + 1;
    color = cmap(color_index, :);
    % 線をプロット
    % plot(x, y_values, '-', 'Color', color, 'LineWidth', 1.5);
    plot(x, y_values, '-', 'Color', color);
    
    % マーカーをプロット
    scatter(x, y_values, 10, color, 'filled');
end

caxis([min(new_standing_table(3:end, 5)), max(new_standing_table(3:end, 5))]);
c = colorbar;
c.Label.String = 'Work [J]';
c.Label.FontSize = 25;
colormap(cmap); % カラーバーにカラーマップを適用
ylim([-0.1, 1.1]);
hold off;
xticks([1, 2, 3, 4]);
xticklabels({'Ci [m]', 'CFLT[m]', 'GEo[m]', 'GE[m]'});
yticks([]);  % y軸のメモリ（目盛り）を削除
yticklabels([]); % y軸のラベルを削除
exportgraphics(gcf, [save_path, 'standing_parallelcoords_color_dev3_1.pdf'], 'ContentType', 'vector');



cmap = colormap('turbo');
figure();
hold on
Ci_count = [1, 1, 1, 1, 1];
CFLT_count = [1, 1, 1, 1, 1];
GEo_count = [1, 1, 1, 1, 1];
GE_count = [1, 1, 1, 1, 1, 1, 1, 1, 1];
% offsetの設定
for i = dev_3+1:dev_3*2
    for j = 1:5
        if new_standing_table(i, 1) == Ci_list_norm(j)
            y_offset_values(i, 1) = y_offset_scale * Ci_count(j);
            Ci_count(j) = Ci_count(j) + 1;
        end
        if new_standing_table(i, 2) == CFLT_list_norm(j)
            y_offset_values(i, 2) = y_offset_scale * CFLT_count(j);
            CFLT_count(j) = CFLT_count(j) + 1;
        end
        if new_standing_table(i, 3) == GEo_list_norm(j)
            y_offset_values(i, 3) = y_offset_scale * GEo_count(j);
            GEo_count(j) = GEo_count(j) + 1;
        end
    end
    for j = 1:9
        if new_standing_table(i, 4) == GE_list_norm(j)
            y_offset_values(i, 4) = y_offset_scale * GE_count(j);
            GE_count(j) = GE_count(j) + 1;
        end
    end   
end
for i = dev_3+1:dev_3*2

    y_values = new_standing_table(i, 1:4) + y_offset_values(i, 1:4)*3;
    
    % 仕事に応じた色を取得
    val_scaled = (new_standing_table(i, 5) - min(new_standing_table(1:end, 5))) / (max(new_standing_table(1:end, 5)) - min(new_standing_table(1:end, 5)));
    color_index = round(val_scaled * (size(cmap, 1) - 1)) + 1;
    color = cmap(color_index, :);
    % 線をプロット
    % plot(x, y_values, '-', 'Color', color, 'LineWidth', 1.5);
    plot(x, y_values, '-', 'Color', color);
    
    % マーカーをプロット
    scatter(x, y_values, 10, color, 'filled');
end

caxis([min(new_standing_table(1:end, 5)), max(new_standing_table(1:end, 5))]);
c = colorbar;
c.Label.String = 'Work [J]';
c.Label.FontSize = 25;
colormap(cmap); % カラーバーにカラーマップを適用
ylim([-0.1, 1.1]);
hold off;
xticks([1, 2, 3, 4]);
xticklabels({'Ci [m]', 'CFLT[m]', 'GEo[m]', 'GE[m]'});
yticks([]);  % y軸のメモリ（目盛り）を削除
yticklabels([]); % y軸のラベルを削除
exportgraphics(gcf, [save_path, 'standing_parallelcoords_color_dev3_2.pdf'], 'ContentType', 'vector');


cmap = colormap('turbo');
figure();
hold on
Ci_count = [1, 1, 1, 1, 1];
CFLT_count = [1, 1, 1, 1, 1];
GEo_count = [1, 1, 1, 1, 1];
GE_count = [1, 1, 1, 1, 1, 1, 1, 1, 1];
% offsetの設定
for i = dev_3*2+1:size(standing_table, 1)
    for j = 1:5
        if new_standing_table(i, 1) == Ci_list_norm(j)
            y_offset_values(i, 1) = y_offset_scale * Ci_count(j);
            Ci_count(j) = Ci_count(j) + 1;
        end
        if new_standing_table(i, 2) == CFLT_list_norm(j)
            y_offset_values(i, 2) = y_offset_scale * CFLT_count(j);
            CFLT_count(j) = CFLT_count(j) + 1;
        end
        if new_standing_table(i, 3) == GEo_list_norm(j)
            y_offset_values(i, 3) = y_offset_scale * GEo_count(j);
            GEo_count(j) = GEo_count(j) + 1;
        end
    end
    for j = 1:9
        if new_standing_table(i, 4) == GE_list_norm(j)
            y_offset_values(i, 4) = y_offset_scale * GE_count(j);
            GE_count(j) = GE_count(j) + 1;
        end
    end   
end
for i = dev_3*2+1:size(standing_table, 1)

    y_values = new_standing_table(i, 1:4) + y_offset_values(i, 1:4)*3;
    
    % 仕事に応じた色を取得
    val_scaled = (new_standing_table(i, 5) - min(new_standing_table(1:end, 5))) / (max(new_standing_table(1:end, 5)) - min(new_standing_table(1:end, 5)));
    color_index = round(val_scaled * (size(cmap, 1) - 1)) + 1;
    color = cmap(color_index, :);
    % 線をプロット
    % plot(x, y_values, '-', 'Color', color, 'LineWidth', 1.5);
    plot(x, y_values, '-', 'Color', color);
    
    % マーカーをプロット
    scatter(x, y_values, 10, color, 'filled');
end

caxis([min(new_standing_table(1:end, 5)), max(new_standing_table(1:end, 5))]);
c = colorbar;
c.Label.String = 'Work [J]';
c.Label.FontSize = 25;
colormap(cmap); % カラーバーにカラーマップを適用
ylim([-0.1, 1.1]);
hold off;
xticks([1, 2, 3, 4]);
xticklabels({'Ci [m]', 'CFLT[m]', 'GEo[m]', 'GE[m]'});
yticks([]);  % y軸のメモリ（目盛り）を削除
yticklabels([]); % y軸のラベルを削除
exportgraphics(gcf, [save_path, 'standing_parallelcoords_color_dev3_3.pdf'], 'ContentType', 'vector');


figure();
cmap = colormap('winter');
new_top10_table = zeros(size(filtered_table(3:end, :), 1), 5);
% top10_standing_tableのmuscle_lengthsの値を0~1に正規化
for i = 1:size(new_top10_table, 1)
    new_top10_table(i, 1) = (filtered_table(i+2, 1) - min(Ci_list)) / (max(Ci_list) - min(Ci_list));
    new_top10_table(i, 2) = (filtered_table(i+2, 2) - min(CFLT_list)) / (max(CFLT_list) - min(CFLT_list));
    new_top10_table(i, 3) = (filtered_table(i+2, 3) - min(GEo_list)) / (max(GEo_list) - min(GEo_list));
    new_top10_table(i, 4) = (filtered_table(i+2, 4) - min(GE_list)) / (max(GE_list) - min(GE_list));
    new_top10_table(i, 5) = filtered_table(i+2, 5);
end
% new_top10_tableを5列目の値で昇順にソート
new_top10_table = sortrows(new_top10_table, 5);

y_offset_scale = 0.0002; % y軸方向のオフセットの幅
y_offset_values = zeros(size(new_top10_table, 1), 4);

Ci_count = [1, 1, 1, 1, 1];
CFLT_count = [1, 1, 1, 1, 1];
GEo_count = [1, 1, 1, 1, 1];
GE_count = [1, 1, 1, 1, 1, 1, 1, 1, 1];
% offsetの設定
for i = 1:size(new_top10_table, 1)
    for j = 1:5
        if new_top10_table(i, 1) == Ci_list_norm(j)
            y_offset_values(i, 1) = y_offset_scale * Ci_count(j);
            Ci_count(j) = Ci_count(j) + 1;
        end
        if new_top10_table(i, 2) == CFLT_list_norm(j)
            y_offset_values(i, 2) = y_offset_scale * CFLT_count(j);
            CFLT_count(j) = CFLT_count(j) + 1;
        end
        if new_top10_table(i, 3) == GEo_list_norm(j)
            y_offset_values(i, 3) = y_offset_scale * GEo_count(j);
            GEo_count(j) = GEo_count(j) + 1;
        end
    end
    for j = 1:9
        if new_top10_table(i, 4) == GE_list_norm(j)
            y_offset_values(i, 4) = y_offset_scale * GE_count(j);
            GE_count(j) = GE_count(j) + 1;
        end
    end   
end

% top10_standing_tableについて各データをプロット
hold on
for i = 1:size(new_top10_table, 1)

    y_values = new_top10_table(i, 1:4) + y_offset_values(i, 1:4)*10;
    
    % 仕事に応じた色を取得
    val_scaled = (new_top10_table(i, 5) - min(new_top10_table(1:end, 5))) / (max(new_top10_table(1:end, 5)) - min(new_top10_table(1:end, 5)));
    color_index = round(val_scaled * (size(cmap, 1) - 1)) + 1;
    color = cmap(color_index, :);
    % 線をプロット
    plot(x, y_values, '-', 'Color', color);
    
    % マーカーをプロット
    scatter(x, y_values, 10, color, 'filled');
end

caxis([min(filtered_table(3:end, 5)), max(filtered_table(3:end, 5))]);
c = colorbar;
c.Label.String = 'Work [J]';
c.Label.FontSize = 25;
colormap(cmap); % カラーバーにカラーマップを適用
ylim([-0.1, 1.1]);
hold off;
xticks([1, 2, 3, 4]);
xticklabels({'Ci [m]', 'CFLT[m]', 'GEo[m]', 'GE[m]'});
yticks([]);  % y軸のメモリ（目盛り）を削除
yticklabels([]); % y軸のラベルを削除
exportgraphics(gcf, [save_path, 'top10_standing_parallelcoords_color.pdf'], 'ContentType', 'vector');


figure();
cmap = colormap('autumn');
cmap = flipud(cmap);
new_not_standing_table = zeros(size(not_standing_table(3:end, :), 1), 5);

% not_standing_tableのmuscle_lengthsの値を0~1に正規化
for i = 1:size(new_not_standing_table, 1)
    new_not_standing_table(i, 1) = (not_standing_table(i+2, 1) - min(Ci_list)) / (max(Ci_list) - min(Ci_list));
    new_not_standing_table(i, 2) = (not_standing_table(i+2, 2) - min(CFLT_list)) / (max(CFLT_list) - min(CFLT_list));
    new_not_standing_table(i, 3) = (not_standing_table(i+2, 3) - min(GEo_list)) / (max(GEo_list) - min(GEo_list));
    new_not_standing_table(i, 4) = (not_standing_table(i+2, 4) - min(GE_list)) / (max(GE_list) - min(GE_list));
    new_not_standing_table(i, 5) = not_standing_table(i+2, 5);
end

% new_not_standing_tableを5列目の値で昇順にソート
new_not_standing_table = sortrows(new_not_standing_table, 5);

y_offset_scale = 0.0002; % y軸方向のオフセットの幅
y_offset_values = zeros(size(new_not_standing_table, 1), 4);

Ci_count = [1, 1, 1, 1, 1];
CFLT_count = [1, 1, 1, 1, 1];
GEo_count = [1, 1, 1, 1, 1];
GE_count = [1, 1, 1, 1, 1, 1, 1, 1, 1];

% offsetの設定
for i = 1:size(new_not_standing_table, 1)
    for j = 1:5
        if new_not_standing_table(i, 1) == Ci_list_norm(j)
            y_offset_values(i, 1) = y_offset_scale * Ci_count(j);
            Ci_count(j) = Ci_count(j) + 1;
        end
        if new_not_standing_table(i, 2) == CFLT_list_norm(j)
            y_offset_values(i, 2) = y_offset_scale * CFLT_count(j);
            CFLT_count(j) = CFLT_count(j) + 1;
        end
        if new_not_standing_table(i, 3) == GEo_list_norm(j)
            y_offset_values(i, 3) = y_offset_scale * GEo_count(j);
            GEo_count(j) = GEo_count(j) + 1;
        end
    end
    for j = 1:9
        if new_not_standing_table(i, 4) == GE_list_norm(j)
            y_offset_values(i, 4) = y_offset_scale * GE_count(j);
            GE_count(j) = GE_count(j) + 1;
        end
    end   
end

% not_standing_tableについて各データをプロット
hold on
for i = 1:size(new_not_standing_table, 1)

    y_values = new_not_standing_table(i, 1:4) + y_offset_values(i, 1:4)*2;
    
    % 仕事に応じた色を取得
    val_scaled = (new_not_standing_table(i, 5) - min(new_not_standing_table(1:end, 5))) / (max(new_not_standing_table(1:end, 5)) - min(new_not_standing_table(1:end, 5)));
    color_index = round(val_scaled * (size(cmap, 1) - 1)) + 1;
    color = cmap(color_index, :);
    % 線をプロット
    plot(x, y_values, '-', 'Color', color);
    
    % マーカーをプロット
    scatter(x, y_values, 10, color, 'filled');
end

caxis([min(not_standing_table(3:end, 5)), max(not_standing_table(3:end, 5))]);
c = colorbar;
c.Label.String = 'Work [J]';
c.Label.FontSize = 25;
colormap(cmap); % カラーバーにカラーマップを適用
ylim([-0.1, 1.1]);
hold off;
xticks([1, 2, 3, 4]);
xticklabels({'Ci [m]', 'CFLT[m]', 'GEo[m]', 'GE[m]'});
yticks([]);  % y軸のメモリ（目盛り）を削除
yticklabels([]); % y軸のラベルを削除
exportgraphics(gcf, [save_path, 'not_standing_parallelcoords_color.pdf'], 'ContentType', 'vector');

% not_standing_tableの中から，膝が過伸展したものを抽出
plot_table = zeros(exp_num, 6);  % 
plot_table(1:end, 1) = muscle_lengths(:,3);
plot_table(1:end, 2) = muscle_lengths(:,4);
plot_table(1:end, 3) = muscle_lengths(:,5);
plot_table(1:end, 4) = muscle_lengths(:,6);
plot_table(1:end, 5) = work(:,2);
% グループ分けのラベル付け (6列目)
plot_table(:,6) = 0;  % デフォルトで除外対象
for i = 1:exp_num
    if knee_angle_end(i,3) == 2
        plot_table(i, 6) = 2;  % 立てず，かつ膝が曲がってる
    elseif knee_angle_end(i,3) == 3
        plot_table(i, 6) = 3;  % 立てず，かつ膝が過伸展
    end
end
not_standing_table_knee_flex = plot_table(plot_table(:,6) == 2, :);
not_standing_table_hyperextension = plot_table(plot_table(:,6) == 3, :);


% not_standing_table_knee_flexを正規化
new_not_standing_table_knee_flex = zeros(size(not_standing_table_knee_flex(1:end, :), 1), 5);
for i = 1:size(new_not_standing_table_knee_flex, 1)
    new_not_standing_table_knee_flex(i, 1) = (not_standing_table_knee_flex(i, 1) - min(Ci_list)) / (max(Ci_list) - min(Ci_list));
    new_not_standing_table_knee_flex(i, 2) = (not_standing_table_knee_flex(i, 2) - min(CFLT_list)) / (max(CFLT_list) - min(CFLT_list));
    new_not_standing_table_knee_flex(i, 3) = (not_standing_table_knee_flex(i, 3) - min(GEo_list)) / (max(GEo_list) - min(GEo_list));
    new_not_standing_table_knee_flex(i, 4) = (not_standing_table_knee_flex(i, 4) - min(GE_list)) / (max(GE_list) - min(GE_list));
    new_not_standing_table_knee_flex(i, 5) = not_standing_table_knee_flex(i, 5);
end
% new_not_standing_table_knee_flexを5列目の値で昇順にソート
new_not_standing_table_knee_flex = sortrows(new_not_standing_table_knee_flex, 5);

y_offset_scale = 0.0002; % y軸方向のオフセットの幅
y_offset_values = zeros(size(new_not_standing_table_knee_flex, 1), 4);

Ci_count = [1, 1, 1, 1, 1];
CFLT_count = [1, 1, 1, 1, 1];
GEo_count = [1, 1, 1, 1, 1];
GE_count = [1, 1, 1, 1, 1, 1, 1, 1, 1];
% offsetの設定
for i = 1:size(new_not_standing_table_knee_flex, 1)
    for j = 1:5
        if new_not_standing_table_knee_flex(i, 1) == Ci_list_norm(j)
            y_offset_values(i, 1) = y_offset_scale * Ci_count(j);
            Ci_count(j) = Ci_count(j) + 1;
        end
        if new_not_standing_table_knee_flex(i, 2) == CFLT_list_norm(j)
            y_offset_values(i, 2) = y_offset_scale * CFLT_count(j);
            CFLT_count(j) = CFLT_count(j) + 1;
        end
        if new_not_standing_table_knee_flex(i, 3) == GEo_list_norm(j)
            y_offset_values(i, 3) = y_offset_scale * GEo_count(j);
            GEo_count(j) = GEo_count(j) + 1;
        end
    end
    for j = 1:9
        if new_not_standing_table_knee_flex(i, 4) == GE_list_norm(j)
            y_offset_values(i, 4) = y_offset_scale * GE_count(j);
            GE_count(j) = GE_count(j) + 1;
        end
    end   
end



% not_standing_table_knee_flexについて各データをプロット
figure();
cmap = colormap('autumn');
cmap = flipud(cmap);
hold on
for i = 1:size(new_not_standing_table_knee_flex, 1)

    y_values = new_not_standing_table_knee_flex(i, 1:4) + y_offset_values(i, 1:4)*5;
    
    % 仕事に応じた色を取得
    val_scaled = (new_not_standing_table_knee_flex(i, 5) - min(new_not_standing_table_knee_flex(1:end, 5))) / (max(new_not_standing_table_knee_flex(1:end, 5)) - min(new_not_standing_table_knee_flex(1:end, 5)));
    color_index = round(val_scaled * (size(cmap, 1) - 1)) + 1;
    color = cmap(color_index, :);
    % 線をプロット
    plot(x, y_values, '-', 'Color', color);
    
    % マーカーをプロット
    scatter(x, y_values, 10, color, 'filled');
end

caxis([min(not_standing_table_knee_flex(1:end, 5)), max(not_standing_table_knee_flex(1:end, 5))]);
c = colorbar;
c.Label.String = 'Work [J]';
c.Label.FontSize = 25;
colormap(cmap); % カラーバーにカラーマップを適用
ylim([-0.1, 1.1]);
hold off;
xticks([1, 2, 3, 4]);
xticklabels({'Ci [m]', 'CFLT[m]', 'GEo[m]', 'GE[m]'});
yticks([]);  % y軸のメモリ（目盛り）を削除
yticklabels([]); % y軸のラベルを削除
exportgraphics(gcf, [save_path, 'not_standing_knee_flex_parallelcoords_color.pdf'], 'ContentType', 'vector');





% not_standing_table_hyperextensionを正規化
new_not_standing_table_hyperextension = zeros(size(not_standing_table_hyperextension(1:end, :), 1), 5);
for i = 1:size(new_not_standing_table_hyperextension, 1)
    new_not_standing_table_hyperextension(i, 1) = (not_standing_table_hyperextension(i, 1) - min(Ci_list)) / (max(Ci_list) - min(Ci_list));
    new_not_standing_table_hyperextension(i, 2) = (not_standing_table_hyperextension(i, 2) - min(CFLT_list)) / (max(CFLT_list) - min(CFLT_list));
    new_not_standing_table_hyperextension(i, 3) = (not_standing_table_hyperextension(i, 3) - min(GEo_list)) / (max(GEo_list) - min(GEo_list));
    new_not_standing_table_hyperextension(i, 4) = (not_standing_table_hyperextension(i, 4) - min(GE_list)) / (max(GE_list) - min(GE_list));
    new_not_standing_table_hyperextension(i, 5) = not_standing_table_hyperextension(i, 5);
end
% new_not_standing_table_hyperextensionを5列目の値で昇順にソート
new_not_standing_table_hyperextension = sortrows(new_not_standing_table_hyperextension, 5);

y_offset_scale = 0.0002; % y軸方向のオフセットの幅
y_offset_values = zeros(size(new_not_standing_table_hyperextension, 1), 4);

Ci_count = [1, 1, 1, 1, 1];
CFLT_count = [1, 1, 1, 1, 1];
GEo_count = [1, 1, 1, 1, 1];
GE_count = [1, 1, 1, 1, 1, 1, 1, 1, 1];
% offsetの設定
for i = 1:size(new_not_standing_table_hyperextension, 1)
    for j = 1:5
        if new_not_standing_table_hyperextension(i, 1) == Ci_list_norm(j)
            y_offset_values(i, 1) = y_offset_scale * Ci_count(j);
            Ci_count(j) = Ci_count(j) + 1;
        end
        if new_not_standing_table_hyperextension(i, 2) == CFLT_list_norm(j)
            y_offset_values(i, 2) = y_offset_scale * CFLT_count(j);
            CFLT_count(j) = CFLT_count(j) + 1;
        end
        if new_not_standing_table_hyperextension(i, 3) == GEo_list_norm(j)
            y_offset_values(i, 3) = y_offset_scale * GEo_count(j);
            GEo_count(j) = GEo_count(j) + 1;
        end
    end
    for j = 1:9
        if new_not_standing_table_hyperextension(i, 4) == GE_list_norm(j)
            y_offset_values(i, 4) = y_offset_scale * GE_count(j);
            GE_count(j) = GE_count(j) + 1;
        end
    end   
end

% not_standing_table_hyperextensionについて各データをプロット
figure();
cmap = colormap('autumn');
% 色の上下を反転
cmap = flipud(cmap);
hold on
for i = 1:size(new_not_standing_table_hyperextension, 1)

    y_values = new_not_standing_table_hyperextension(i, 1:4) + y_offset_values(i, 1:4)*2;
    
    % 仕事に応じた色を取得
    val_scaled = (new_not_standing_table_hyperextension(i, 5) - min(new_not_standing_table_hyperextension(1:end, 5))) / (max(new_not_standing_table_hyperextension(1:end, 5)) - min(new_not_standing_table_hyperextension(1:end, 5)));
    color_index = round(val_scaled * (size(cmap, 1) - 1)) + 1;
    color = cmap(color_index, :);
    % 線をプロット
    plot(x, y_values, '-', 'Color', color);
    
    % マーカーをプロット
    scatter(x, y_values, 10, color, 'filled');
end

caxis([min(not_standing_table_hyperextension(1:end, 5)), max(not_standing_table_hyperextension(1:end, 5))]);
c = colorbar;
c.Label.String = 'Work [J]';
c.Label.FontSize = 25;
colormap(cmap); % カラーバーにカラーマップを適用
ylim([-0.1, 1.1]);
hold off;
xticks([1, 2, 3, 4]);
xticklabels({'Ci [m]', 'CFLT[m]', 'GEo[m]', 'GE[m]'});
yticks([]);  % y軸のメモリ（目盛り）を削除
yticklabels([]); % y軸のラベルを削除
exportgraphics(gcf, [save_path, 'not_standing_hyperextension_parallelcoords_color.pdf'], 'ContentType', 'vector');


% % not_standing_indexから膝が過伸展したものを抽出し，別の配列に保存．not_standing_indexからは除外
% hyperextension_index = [];
% for i = 1:exp_num
%     if knee_angle_end(i,3) == 3
%         % hyperextension_indexの末尾に追加
%         hyperextension_index = [hyperextension_index, i];
%         % not_standing_indexから除外
%         not_standing_index = not_standing_index(not_standing_index ~= i);
%     end
% end
% hyperextension_index = hyperextension_index';
% save('top10_standing_not_standing_MuscleLengthTest_PID.mat', 'top10_standing_tension_index', 'not_standing_index', 'standing_index', 'hyperextension_index');