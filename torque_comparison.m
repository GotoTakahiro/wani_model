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
line_width = 1;
figure()
hold on
for i = 1:exp_num
    if ismember(i,top10_standing_tension_index)
        p1 = plot(hip_angle_all(i,start_num:end),hip_torque_all(i,start_num:end),'LineWidth', line_width, "Color","#77AC30");
    elseif ismember(i,standing_index)
        p2 = plot(hip_angle_all(i,start_num:end),hip_torque_all(i,start_num:end),'LineWidth', line_width,"Color","#0072BD");
        p2.Color(4) = 0.3;
    elseif ismember(i,not_standing_index)
        p3 = plot(hip_angle_all(i,start_num:end),hip_torque_all(i,start_num:end),'LineWidth', line_width,"Color","#D95319");
        p3.Color(4) = 0.3;
    elseif ismember(i,hyperextension_index)
        p4 = plot(hip_angle_all(i,start_num:end),hip_torque_all(i,start_num:end),'LineWidth', line_width,"Color","#7E2F8E");
        p4.Color(4) = 0.3;
    end
end
set(gca, 'XDir', 'reverse');
hold off
xlim([-100, -10]);
ylim([-40, 0]);
h_axes = gca;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;
xlabel('Hip angle [deg]', 'FontSize', 25);
ylabel('Hip torque [Nm]','FontSize', 25);
% legend([p1 p2 p3 p4], '(i)', '(ii)', '(iii)', '(iv)', 'FontSize', 25, 'Location', 'best');
legend([p1 p2 p3], '(i)', '(ii)', '(iii)', 'FontSize', 25, 'Location', 'best');
if graph_save
    exportgraphics(gca, [save_path 'hip_torque_comparison.pdf'], 'ContentType', 'vector');
    saveas(gcf, [save_path 'hip_torque_comparison.png']);
end



figure()
hold on
for i = 1:exp_num
    if ismember(i,top10_standing_tension_index)
        p1 = plot(hip_angle_all(i,start_num:end),knee_torque_all(i,start_num:end),'LineWidth', line_width, "Color","#77AC30");
    elseif ismember(i,standing_index)
        p2 = plot(hip_angle_all(i,start_num:end),knee_torque_all(i,start_num:end),'LineWidth', line_width,"Color","#0072BD");
        p2.Color(4) = 0.3;
    elseif ismember(i,not_standing_index)
        p3 = plot(hip_angle_all(i,start_num:end),knee_torque_all(i,start_num:end),'LineWidth', line_width,"Color","#D95319");
        p3.Color(4) = 0.3;
    elseif ismember(i,hyperextension_index)
        p4 = plot(hip_angle_all(i,start_num:end),knee_torque_all(i,start_num:end),'LineWidth', line_width,"Color","#7E2F8E");
        p4.Color(4) = 0.3;
    end
end
set(gca, 'XDir', 'reverse');
hold off
xlim([-100, -10]);
ylim([-18, 0]);
h_axes = gca;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;
xlabel('Hip angle [deg]','FontSize', 25);
ylabel('Knee torque [Nm]','FontSize', 25);
% legend([p1 p2 p3 p4], '(i)', '(ii)', '(iii)', '(iv)', 'FontSize', 25, 'Location', 'best');
legend([p1 p2 p3], '(i)', '(ii)', '(iii)', 'FontSize', 25, 'Location', 'best');
if graph_save
    exportgraphics(gca, [save_path 'knee_torque_comparison.pdf'], 'ContentType', 'vector');
    saveas(gcf, [save_path 'knee_torque_comparison.png']);
end


hip_knee_torque_ratio = knee_torque_all./hip_torque_all;
figure()
hold on
for i = 1:exp_num
    if ismember(i,top10_standing_tension_index)
        p1 = plot(hip_angle_all(i,start_num:end),hip_knee_torque_ratio(i,start_num:end),'LineWidth', line_width,"Color","#77AC30");
    elseif ismember(i,standing_index)
        p2 = plot(hip_angle_all(i,start_num:end),hip_knee_torque_ratio(i,start_num:end),'LineWidth', line_width,"Color","#0072BD");
        p2.Color(4) = 0.3;
    elseif ismember(i,not_standing_index)
        p3 = plot(hip_angle_all(i,start_num:end),hip_knee_torque_ratio(i,start_num:end),'LineWidth', line_width,"Color","#D95319");
        p3.Color(4) = 0.3;
    elseif ismember(i,hyperextension_index)
        p4 = plot(hip_angle_all(i,start_num:end),hip_knee_torque_ratio(i,start_num:end),'LineWidth', line_width,"Color","#7E2F8E");
        p4.Color(4) = 0.3;
    end
end
set(gca, 'XDir', 'reverse');
hold off
xlim([-100, -10]);
ylim([0, 0.7]);
h_axes = gca;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;
xlabel('Hip angle [deg]', 'FontSize', 25);
ylabel('Knee/Hip torque ratio [-]', 'FontSize', 25);
% legend([p1 p2 p3 p4], '(i)', '(ii)', '(iii)', '(iv)', 'FontSize', 25, 'Location', 'best');
legend([p1 p2 p3], '(i)', '(ii)', '(iii)', 'FontSize', 25, 'Location', 'best');
if graph_save
    exportgraphics(gca, [save_path 'knee_hip_torque_ratio_comparison.pdf'], 'ContentType', 'vector');
    saveas(gcf, [save_path 'knee_hip_torque_ratio_comparison.png']);
end