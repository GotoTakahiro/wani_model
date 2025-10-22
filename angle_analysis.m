clear
load('results/20240712_MuscleLengthTest_PID/exp20240712_MuscleLengthTest_PID_length_and_gain_combination.mat');
save_path = 'results/20240712_MuscleLengthTest_PID/';

graph_view = true;
graph_save = true;

exp_num = size(length_and_gain_combination,1);

tension_flag = false;
tension_start_index = zeros(exp_num,1);

t_tension = cell(exp_num,1);
alpha_CFLT_sum = cell(exp_num,1);
alpha_GEo_sum = cell(exp_num,1);
alpha_GE_sum = cell(exp_num,1);
alpha_Ci_sum = cell(exp_num,1);
angle_femur_GEo = cell(exp_num,1);
ratio_CFLT_GEo = cell(exp_num,1);
y_hip_end = zeros(exp_num,1);
standing = zeros(exp_num,2);
plot_flag_CFLT = false;
plot_flag_GEo = false;
plot_flag_GE = false;
plot_flag_Ci = false;
plot_flag_ratio = false;
plot_flag_angle_femur_GEo = false;

work = zeros(exp_num,2);

for i = 1:exp_num
    L_CFL = length_and_gain_combination(i,1);
    L_Ci = length_and_gain_combination(i,2);
    L_CFLT = length_and_gain_combination(i,3);
    L_GEo = length_and_gain_combination(i,4);
    L_GE = length_and_gain_combination(i,5);
    Pgain = length_and_gain_combination(i,6);
    Igain = length_and_gain_combination(i,7);
    Dgain = length_and_gain_combination(i,8);
    clear t q muscle_tension torque_muscle L_wire k_wire c_wire general_q general_dq l_muscle_list l_link_list data_Q power coordinates_x coordinates_y angle_wire_pulley 
    filename = sprintf('results/20240712_MuscleLengthTest_PID/exp20240712_MuscleLengthTest_PID_%d_P%d_I%d_D%d_CFL%d_Ci%d_CFLT%d_GEo%d_GE%d.mat',i,Pgain,Igain,Dgain,L_CFL*1000,L_Ci*1000,L_CFLT*1000,L_GEo*1000,L_GE*1000);
    load(filename);

    [coordinates_x, coordinates_y, angle_wire_pulley] = calc_coordinate_for_plot(q,r,l_link_list);
    y_hip_end(i) = coordinates_y(1000,2);
    standing(i,1) = i;
    if y_hip_end(i) > 0
        standing(i,2) = 1;
    end

    % -data_Q(:,11)が5以上になったところから最後までを取得
    for j = 300:size(data_Q,1)
        if tension_flag == false && -data_Q(j,11) > 5
            tension_start_index(i) = j;
            tension_flag = true;
        end
    end
    

    tension_flag = false;
    coordinates_x = coordinates_x(tension_start_index(i):end,:);
    coordinates_y = coordinates_y(tension_start_index(i):end,:);
    angle_wire_pulley = angle_wire_pulley(tension_start_index(i):end,:);

    t_tension{i} = t(tension_start_index(i):end,1) - t(tension_start_index(i),1);
    alpha_CFLT_sum{i} = calc_alpha_CFLT(coordinates_x, coordinates_y);
    alpha_GEo_sum{i} = calc_alpha_GEo(coordinates_x, coordinates_y);
    alpha_Ci_sum{i} = calc_alpha_Ci(coordinates_x, coordinates_y);
    alpha_GE_sum{i} = calc_alpha_GE(coordinates_x, coordinates_y);
    ratio_CFLT_GEo{i} = alpha_CFLT_sum{i}./alpha_GEo_sum{i};
    angle_femur_GEo{i} = calc_angle_femur_GEo(coordinates_x, coordinates_y);

    % flag = false;
    % for j = 300:400
    %     if alpha_GE_sum{i}(j) < 10 && flag == false && standing(i,2) == 0
    %         disp(i);
    %         flag = true;
    %     end
    % end
    % flag = false;

    power = q(50:end,20).*data_Q(50:end,11);
    work(i,1) = i;
    work(i,2) = trapz(t(50:end,1),power);
end

%仕事が小さく，かつ立ち上がることのできた上位30個のstanding(i,1)を2にする
% 仕事が小さく、かつ立ち上がることのできた上位30個のstanding(i,1)を2にする
standing_index = find(standing(:,2) == 1);
[~,sorted_index] = sort(work(standing_index,2));
% 取り出したい上位30個のインデックス
top_30_index = standing_index(sorted_index(1:30));
% standing_indexを利用して直接standingを更新
standing(top_30_index, 2) = 2;

% for i = 1:exp_num
%     if standing(i,2) == 2
%         for j = 100:200
%             if alpha_CFLT_sum{i}(j) < 110
%                 disp(i);
%                 disp(y_hip_end(i));
%                 break;
%             end
%         end
%     end
% end



% figure(1)
% for i = 1:exp_num
%     if standing(i,2) == 1
%         plot(t_tension{i},alpha_CFLT_sum{i},'Color','#0072BD'); %立てたら青
%         if plot_flag_CFLT == false
%             hold on
%             plot_flag_CFLT = true;
%         end
%     end
% end
% for i = 1:exp_num
%     if standing(i,2) == 0
%         plot(t_tension{i},alpha_CFLT_sum{i},'Color','#D95319'); %立てなかったら赤
%     end
% end
% xlabel('Time [s]');
% ylabel('alpha CFLT [deg]');
% if graph_save
%     exportgraphics(gca, [save_path 'alpha_CFLT.pdf'],'ContentType','vector');
% end
% plot_flag_CFLT = false;

% figure(2)
% for i = 1:exp_num
%     if standing(i,2) == 1
%         plot(t_tension{i},alpha_GEo_sum{i},'Color','#0072BD'); %立てたら青
%         if plot_flag_GEo == false
%             hold on
%             plot_flag_GEo = true;
%         end
%     end
% end
% for i = 1:exp_num
%     if standing(i,2) == 0
%         plot(t_tension{i},alpha_GEo_sum{i},'Color','#D95319'); %立てなかったら赤
%     end
% end
% xlabel('Time [s]');
% ylabel('alpha GEo [deg]');
% if graph_save
%     exportgraphics(gca, [save_path 'alpha_GEo.pdf'],'ContentType','vector');
% end
% plot_flag_GEo = false;

% figure(3)
% for i = 1:exp_num
%     if standing(i,2) == 1
%         plot(t_tension{i},alpha_GE_sum{i},'Color','#0072BD'); %立てたら青
%         if plot_flag_GE == false
%             hold on
%             plot_flag_GE = true;
%         end
%     end
% end
% for i = 1:exp_num
%     if standing(i,2) == 0
%         plot(t_tension{i},alpha_GE_sum{i},'Color','#D95319'); %立てなかったら赤
%     end
% end
% xlabel('Time [s]');
% ylabel('alpha GE [deg]');
% if graph_save
%     exportgraphics(gca, [save_path 'alpha_GE.pdf'],'ContentType','vector');
% end
% plot_flag_GE = false;

% figure(4)
% for i = 1:exp_num
%     if standing(i,2) == 1
%         plot(t_tension{i},alpha_Ci_sum{i},'Color','#0072BD'); %立てたら青
%         if plot_flag_Ci == false
%             hold on
%             plot_flag_Ci = true;
%         end
%     end
% end
% for i = 1:exp_num
%     if standing(i,2) == 0
%         plot(t_tension{i},alpha_Ci_sum{i},'Color','#D95319'); %立てなかったら赤
%     end
% end
% xlabel('Time [s]');
% ylabel('alpha Ci [deg]');
% if graph_save
%     exportgraphics(gca, [save_path 'alpha_Ci.pdf'],'ContentType','vector');
% end
% plot_flag_Ci = false;

% figure(5)
% for i = 1:exp_num
%     if standing(i,2) == 2
%         plot(t_tension{i},alpha_CFLT_sum{i},'Color','#0072BD');
%         if plot_flag_CFLT == false
%             hold on
%             plot_flag_CFLT = true;
%         end
%     end
% end
% xlabel('Time [s]');
% ylabel('alpha CFLT [deg]');
% if graph_save
%     exportgraphics(gca, [save_path 'alpha_CFLT_top30.pdf'],'ContentType','vector');
% end
% plot_flag_CFLT = false;

% figure(6)
% for i = 1:exp_num
%     if standing(i,2) == 2
%         plot(t_tension{i},alpha_GEo_sum{i},'Color','#0072BD');
%         if plot_flag_GEo == false
%             hold on
%             plot_flag_GEo = true;
%         end
%     end
% end
% xlabel('Time [s]');
% ylabel('alpha GEo [deg]');
% if graph_save
%     exportgraphics(gca, [save_path 'alpha_GEo_top30.pdf'],'ContentType','vector');
% end


% figure(3)
% for i = 1:exp_num
%     if standing(i) == 1
%         plot(t_tension{i},ratio_CFLT_GEo{i},'Color','#0072BD'); %立てたら青
%         if plot_flag_ratio == false
%             hold on
%             plot_flag_ratio = true;
%         end
%     end
% end
% for i = 1:exp_num
%     if standing(i) == 0
%         plot(t_tension{i},ratio_CFLT_GEo{i},'Color','#D95319'); %立てなかったら赤
%     end
% end
% xlabel('Time [s]');
% ylabel('alpha CFLT / alpha GEo');
% if graph_save
%     exportgraphics(gca, [save_path 'ratio_CFLT_GEo.pdf'],'ContentType','vector');
% end


figure
for i = 1:exp_num
    if standing(i,2) == 1
        plot(t_tension{i},angle_femur_GEo{i},'Color','#0072BD'); %立てたら青
        if plot_flag_angle_femur_GEo == false
            hold on
            plot_flag_angle_femur_GEo = true;
        end
    end
end
for i = 1:exp_num
    if standing(i,2) == 0
        plot(t_tension{i},angle_femur_GEo{i},'Color','#D95319'); %立てなかったら赤
    end
end
for i = 1:exp_num
    if standing(i,2) == 2
        plot(t_tension{i},angle_femur_GEo{i},'Color','#77AC30','LineWidth',2); %立てたら緑
    end
end
xlabel('Time [s]');
ylabel('angle femur GEo [deg]');
if graph_save == true
    exportgraphics(gca, [save_path 'angle_femur_GEo.pdf'],'ContentType','vector');
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