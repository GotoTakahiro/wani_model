clear

load('results/20240724_init_condition_test_MinWork/exp20240724_HighWalk_init_condition.mat');
load('results/20240724_init_condition_test_MinWork/exp20240724_init_condition_test_Minwork_length_and_gain_combination.mat');
save_path = 'results/20240724_init_condition_test_MinWork/';

L_CFL = length_and_gain_combination(1,1);
L_Ci = length_and_gain_combination(1,2);
L_CFLT = length_and_gain_combination(1,3);
L_GEo = length_and_gain_combination(1,4);
L_GE = length_and_gain_combination(1,5);

standing_theta2 = [];
standing_theta3 = [];
not_standing_theta2 = [];
not_standing_theta3 = [];

init_theta2_ = 10:3:70;
init_theta2_ = -init_theta2_/180*pi;
init_theta3_ = 40:3:100;
init_theta3_ = -init_theta3_/180*pi;
[init_theta2_list, init_theta3_list] = ndgrid(init_theta2_, init_theta3_);
init_condition_list_for_lim = [init_theta2_list(:), init_theta3_list(:)];
ankle_lim_theta2 = [];
ankle_lim_theta3 = [];
r = 0.015;
L_met = 90.739/1000;
limit_ankle_ex = 120/180*pi;

for i = size(init_condition_list_for_lim,1):-1:1
    theta2_ = init_condition_list_for_lim(i,1);
    theta3_ = init_condition_list_for_lim(i,2);
    theta4_ = 90/180*pi-(90/180*pi+theta2_+theta3_)-(asin(r/L_met)-0.03);
    if theta4_ > limit_ankle_ex %初期の足関節角度がlimit_ankle_exよりも大きくなる時は除外
        ankle_lim_theta2(end+1) = rad2deg(theta2_);
        ankle_lim_theta3(end+1) = rad2deg(theta3_);
    end
end

for i = 1:size(init_condition_list,1)
    theta2 = init_condition_list(i,1); %hip
    theta3 = init_condition_list(i,2); %knee
    
    clear t q muscle_tension torque_muscle L_wire k_wire c_wire general_q general_dq l_muscle_list l_link_list data_Q power ankle_lim_index
    filename = sprintf('results/20240724_init_condition_test_MinWork/exp20240724_init_condition_test_Minwork_%d_Hip%d_Knee%d_CFL%d_Ci%d_CFLT%d_GEo%d_GE%d.mat', i, -rad2deg(theta2), -rad2deg(theta3), L_CFL*1000, L_Ci*1000, L_CFLT*1000, L_GEo*1000, L_GE*1000);
    load(filename);

    % t = 10s のときのy座標
    y_hip = q(1000,2)-l_link_list(6)*cos(q(1000,5));

    if y_hip > 0
        standing_theta2(end+1) = rad2deg(theta2);
        standing_theta3(end+1) = rad2deg(theta3);
    elseif y_hip <= 0
        not_standing_theta2(end+1) = rad2deg(theta2);
        not_standing_theta3(end+1) = rad2deg(theta3);
    end
end

min_theta2 = rad2deg(min(init_condition_list(:,1)));
max_theta2 = rad2deg(max(init_condition_list(:,1)));
min_theta3 = rad2deg(min(init_condition_list(:,2)));
max_theta3 = rad2deg(max(init_condition_list(:,2)));

figure
scatter(standing_theta2,standing_theta3,'filled')
hold on
scatter(not_standing_theta2,not_standing_theta3,'filled')
%黒でプロット
scatter(ankle_lim_theta2,ankle_lim_theta3,'filled','k')
xlim([min_theta2-5 max_theta2+5])
ylim([min_theta3-5 max_theta3+5])
xlabel('Hip angle [deg]')
ylabel('Knee angle [deg]')
legend('Standing','Not standing','Ankle limit')
exportgraphics(gca, [save_path 'stand_condition_plot.pdf'], 'ContentType', 'vector');