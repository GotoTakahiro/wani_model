clear
load('results/20240718_HighWalk_init/exp20240718_HighWalk_init_PID_length_and_gain_combination.mat');
save_path = 'results/20240718_HighWalk_init/';
% load('results/20240712_MuscleLengthTest_PID/exp20240712_MuscleLengthTest_PID_length_and_gain_combination.mat');
% save_path = 'results/20240712_MuscleLengthTest_PID/';

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
    filename = sprintf('results/20240718_HighWalk_init/exp20240718_HighWalk_init_PID_%d_P%d_I%d_D%d_CFL%d_Ci%d_CFLT%d_GEo%d_GE%d.mat',i,Pgain,Igain,Dgain,L_CFL*1000,L_Ci*1000,L_CFLT*1000,L_GEo*1000,L_GE*1000);
    % filename = sprintf('results/20240712_MuscleLengthTest_PID/exp20240712_MuscleLengthTest_PID_%d_P%d_I%d_D%d_CFL%d_Ci%d_CFLT%d_GEo%d_GE%d.mat',i,Pgain,Igain,Dgain,L_CFL*1000,L_Ci*1000,L_CFLT*1000,L_GEo*1000,L_GE*1000);
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


% 全ての条件をプロット
% plot_table = zeros(exp_num,6);
% % % top_standing_tension_indexに含まれるものは'top30', それ以外の立てたものは'standing', 立てないものは'not standing'として，plot_tableの6列目に格納
% plot_table(:,1) = muscle_lengths(:,3);
% plot_table(:,2) = muscle_lengths(:,4);
% plot_table(:,3) = muscle_lengths(:,5);
% plot_table(:,4) = muscle_lengths(:,6);
% plot_table(:,5) = work(:,2);
% plot_table(:,6) = 1;
% plot_table(top_standing_tension_index,6) = 2;
% plot_table(not_standing_index,6) = 0;

% coorddata = [1 2 3 4 5];
% % figure(12)
% figure();
% % p = parallelplot(plot_table, 'CoordinateData', coorddata, 'GroupData', plot_table(:,6));
% p = parallelplot(plot_table(top_standing_tension_index,:), 'CoordinateData', coorddata);
% p.CoordinateTickLabels = {'Ci [m]', 'CFLT [m]', 'GEo [m]', 'GE [m]', 'Work [J]'};
% % p.Color = {'#A2142F','#0072BD','#EDB120'};
% p.Color = {'#0072BD'};
% if graph_save == true
%     exportgraphics(gca, [save_path, 'min_work_parallelcoords.pdf'], 'ContentType', 'vector');
% end


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
% plot_table(1, 1:5) = [max(muscle_lengths(:,3)), max(muscle_lengths(:,4)), max(muscle_lengths(:,5)), max(muscle_lengths(:,6)), max(work(:,2))];
% plot_table(2, 1:5) = [min(muscle_lengths(:,3)), min(muscle_lengths(:,4)), min(muscle_lengths(:,5)), min(muscle_lengths(:,6)), min(work(:,2))];
plot_table(1, 1:5) = [max(muscle_lengths(:,3)), max(muscle_lengths(:,4)), max(muscle_lengths(:,5)), max(muscle_lengths(:,6)), 5.2];
plot_table(2, 1:5) = [min(muscle_lengths(:,3)), min(muscle_lengths(:,4)), min(muscle_lengths(:,5)), min(muscle_lengths(:,6)), 0];
plot_table(1:2, 6) = 2;  % 新しいグループID（最大・最小用）
% 'top30'と最大最小のみをフィルタリング
filtered_table = plot_table(plot_table(:,6) > 0, :);
% parallelplotの作成
figure();
p = parallelplot(filtered_table, 'CoordinateData', coorddata, 'GroupData', filtered_table(:,6), 'LineWidth', 1.5);
p.CoordinateTickLabels = {'Ci [m]', 'CFLT [m]', 'GEo [m]', 'GE [m]', 'Work [J]'};
% 各グループの色と透明度設定
p.Color = [1 1 1; 0.4660 0.6740 0.1880];  % 'top30'と'最大最小' の色設定
% 図の保存
if graph_save == true
    filename = ['top' num2str(top_num) '_min_max_parallelcoords.pdf'];
    exportgraphics(gcf, [save_path, filename], 'ContentType', 'vector');
end


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
% plot_table(1, 1:5) = [max(muscle_lengths(:,3)), max(muscle_lengths(:,4)), max(muscle_lengths(:,5)), max(muscle_lengths(:,6)), max(work(:,2))];
% plot_table(2, 1:5) = [min(muscle_lengths(:,3)), min(muscle_lengths(:,4)), min(muscle_lengths(:,5)), min(muscle_lengths(:,6)), min(work(:,2))];
plot_table(1, 1:5) = [max(muscle_lengths(:,3)), max(muscle_lengths(:,4)), max(muscle_lengths(:,5)), max(muscle_lengths(:,6)), 5.2];
plot_table(2, 1:5) = [min(muscle_lengths(:,3)), min(muscle_lengths(:,4)), min(muscle_lengths(:,5)), min(muscle_lengths(:,6)), 0];
plot_table(1:2, 6) = 3;  % 新しいグループID（最大・最小用）
% 'not standing'のみのデータをフィルタリング
not_standing_table = plot_table(plot_table(:,6) == 1 | plot_table(:,6) == 3, :);

% 'standing'のみのデータをフィルタリング
standing_table = plot_table(plot_table(:,6) == 2 | plot_table(:,6) == 3, :);
% 'not standing'のみのparallelplot
figure();
p1 = parallelplot(not_standing_table(:,1:5), 'CoordinateData', coorddata, 'GroupData', not_standing_table(:,6));
p1.CoordinateTickLabels = {'Ci [m]', 'CFLT [m]', 'GEo [m]', 'GE [m]', 'Work [J]'};
% 色の設定： 'not standing'は赤、最大・最小は白
p1.Color = [1 1 1; 0.6350 0.0780 0.1840];  % 赤色と白色
% 'not standing'グラフの保存
if graph_save == true
    exportgraphics(gcf, [save_path, 'not_standing_with_limits_parallelcoords.pdf'], 'ContentType', 'vector');
end


% 'standing'のみのparallelplot
figure();
p2 = parallelplot(standing_table(:,1:5), 'CoordinateData', coorddata, 'GroupData', standing_table(:,6));
p2.CoordinateTickLabels = {'Ci [m]', 'CFLT [m]', 'GEo [m]', 'GE [m]', 'Work [J]'};
% 色の設定： 'standing'は柔らかい緑、最大・最小は白
p2.Color = [1 1 1; 0 0.4470 0.7410];  % 緑色と白色
% 'standing'グラフの保存
if graph_save == true
    exportgraphics(gcf, [save_path, 'standing_with_limits_parallelcoords.pdf'], 'ContentType', 'vector');
end
% not_standing_indexから膝が過伸展したものを抽出し，別の配列に保存．not_standing_indexからは除外
hyperextension_index = [];
for i = 1:exp_num
    if knee_angle_end(i,3) == 3
        % hyperextension_indexの末尾に追加
        hyperextension_index = [hyperextension_index, i];
        % not_standing_indexから除外
        not_standing_index = not_standing_index(not_standing_index ~= i);
    end
end
hyperextension_index = hyperextension_index';
save('top10_standing_not_standing_MuscleLengthTest_PID.mat', 'top10_standing_tension_index', 'not_standing_index', 'standing_index', 'hyperextension_index');

% knee_angle_end, work, max_torqueをparallelplotでプロット．top10は緑，それ以外の立てたものは青，立てないものは赤
plot_table = zeros(exp_num+2, 6);  % +2は最大値と最小値用
plot_table(3:end,1) = abs(max_torque(:,2));
plot_table(3:end,2) = abs(max_torque(:,3));
plot_table(3:end,3) = abs(max_torque(:,4));
plot_table(3:end,4) = work(:,2);
plot_table(3:end,5) = knee_angle_end(:,2)*180/pi;  % rad -> deg
% グループ分けのラベル付け (6列目)
% knee_angle_end(i,3)を用いて分ける．top10を1，立ててかつknee_angle_end(i,3)が1のものを2，立てないかつknee_angle_end(i,3)が2のものを3，それ以外を0とする
plot_table(:,6) = 0;  % デフォルトで除外対象
plot_table(top10_standing_tension_index+2, 6) = 1;  % 'top10'
for i = 1:exp_num
    if knee_angle_end(i,3) == 2
        plot_table(i+2, 6) = 2;  % 立てず，かつ膝が曲がってる
    elseif knee_angle_end(i,3) == 3
        plot_table(i+2, 6) = 3;  % 立てず，かつ膝が過伸展
    end
end
max_table = [40, 16, 1.5, 5.5, 60];
min_table = [0, 0, 0, 0, -100];
plot_table(1, 1:5) = max_table;
plot_table(2, 1:5) = min_table;
plot_table(1:2, 6) = 4;  %　最大・最小用のグループID

% 最大最小は白で
coorddata = [1 2 3 4 5];
figure();
p = parallelplot(plot_table, 'CoordinateData', coorddata, 'GroupData', plot_table(:,6));
p.CoordinateTickLabels = {'Max hip torque [Nm]', 'Max knee torque [Nm]', 'Max ankle torque [Nm]', 'Work [J]', 'Knee angle at end [deg]'};
% p.Color = {'#0072BD','#77AC30','#D95319','#7E2F8E'};
% p.Color = {"#FFFFFF",'#D95319','#0072BD','#7E2F8E','#77AC30'};
p.Color = {"#FFFFFF",'#D95319','#0072BD','#77AC30','#77AC30'};
if graph_save == true
    exportgraphics(gca, [save_path, 'min_work_knee_angle_end_work_max_torque_parallelcoords.pdf'], 'ContentType', 'vector');
end


%仕事が最小であるインデックスとワイヤ長を横に並べて表示
disp('Min work index and wire length');
format long
disp([work(top10_standing_tension_index, 1), work(top10_standing_tension_index, 2), muscle_lengths(top10_standing_tension_index, 2:6)]);

%仕事が最大であるインデックスとワイヤ長を横に並べて表示
disp('Max work index and wire length');
disp([work(top10_standing_work_index, 1), work(top10_standing_work_index, 2), muscle_lengths(top10_standing_work_index, 2:6)]);
format short