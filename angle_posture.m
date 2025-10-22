clear all
L_fem = 160.563/1000;
L_tib = 166.190/1000;
L_met = 90.739/1000;
L_frame = 0.30;
L_4th_troch = 49.597/1000;
L_GE_origin = 152.265/1000;
% L_fem = 0.18;
% L_tib = 0.17;
% L_met = 0.1;
% L_frame = 0.30;
% L_4th_troch = 0.18*1/3;
% L_GE_origin = 0.18*0.95;
%プーリ半径
r = 0.015;
l_link_list = zeros(7,1);
l_link_list(1) = L_fem;
l_link_list(2) = L_tib;
l_link_list(3) = L_met;
l_link_list(4) = L_4th_troch;
l_link_list(5) = L_GE_origin;
l_link_list(6) = L_frame;
l_link_list(7) = r;

phi = linspace(0,2*pi,100);
init_theta2_ = 10:2:70;
init_theta2_ = -init_theta2_/180*pi;
init_theta3_ = 30:2:110;
init_theta3_ = -init_theta3_/180*pi;
[init_theta2, init_theta3] = ndgrid(init_theta2_, init_theta3_);
init_condition_list = [init_theta2(:), init_theta3(:)];
disp(size(init_condition_list))

theta4_lim = 120/180*pi;

for i = size(init_condition_list,1):-1:1
    theta2_ = init_condition_list(i,1);
    theta3_ = init_condition_list(i,2);
    theta4_ = 90/180*pi-(90/180*pi+theta2_+theta3_)-(asin(r/L_met)-0.03);
    if theta4_ > theta4_lim
        %条件を満たす要素を削除
        init_condition_list(i,:) = [];
        % init_condition_list(i,2) = [];
    end
end
disp(size(init_condition_list))
% init_condition_list = [init_theta2(:) init_theta3(:)];
% [init_theta2, init_theta3] = ndgrid(init_theta2_, init_theta3_);
% disp(size(init_theta3_));
% init_condition_list = [init_theta2(:) init_theta3(:)];

%ワイヤ長さ
L_Ci = 0.044;
L_CFLT = 0.107;
L_GEo = 0.037;
L_GE = 0.188;

x_lim_neg = -0.5;
x_lim_pos = 0.5;
y_lim_neg = -0.5;
y_lim_pos = 0.5;

disp(size(init_condition_list,1))
v = VideoWriter('angle_posture.mp4','MPEG-4');
v.Quality = 100;
v.FrameRate = 100;
open(v)
ax = figure;

for i = 1:size(init_condition_list,1)
    theta1 = 90/180*pi;
    theta2 = init_condition_list(i,1);
    theta3 = init_condition_list(i,2);
    theta4 = 90/180*pi-(theta1+theta2+theta3)-(asin(r/L_met)-0.03);
    % theta_CFL = 85/180*pi;
    theta_CFL = 85/180*pi+theta2/5;
    L_CFL = 0.35;
    x1 = -0.30;
    y1 = 0;

    x2 = x1 + L_CFL*sin(theta_CFL);
    y2 = y1 - L_CFL*cos(theta_CFL);
    x3 = x2 + L_CFLT*sin(theta1+theta2);
    y3 = y2 - L_CFLT*cos(theta1+theta2);
    
    x_hip = x1 + L_frame*sin(theta1);
    y_hip = y1 - L_frame*cos(theta1);
    x_knee = x_hip + L_fem*sin(theta1+theta2);
    y_knee = y_hip - L_fem*cos(theta1+theta2);
    x_ankle = x_knee + L_tib*sin(theta1+theta2+theta3);
    y_ankle = y_knee - L_tib*cos(theta1+theta2+theta3);
    x_toe = x_ankle + L_met*sin(theta1+theta2+theta3+theta4);
    y_toe = y_ankle - L_met*cos(theta1+theta2+theta3+theta4);

    x_Ci = x_hip + L_4th_troch*sin(theta1+theta2);
    y_Ci = y_hip - L_4th_troch*cos(theta1+theta2);
    x_GE_origin = x_hip + L_GE_origin*sin(theta1+theta2);
    y_GE_origin = y_hip - L_GE_origin*cos(theta1+theta2);

    plot([x1 x_hip],[y1 y_hip],'-','MarkerSize',2,'LineWidth',2); %CFL_origin to hip
    hold on
    plot([x_hip x_knee],[y_hip y_knee],'-','MarkerSize',2,'LineWidth',2); %hip to knee
    plot([x_knee x_ankle],[y_knee y_ankle],'-','MarkerSize',2,'LineWidth',2); %knee to ankle
    plot([x_ankle x_toe],[y_ankle y_toe],'-','MarkerSize',2,'LineWidth',2); %ankle to toe
    plot((r*cos(phi)+x_ankle) , (r*sin(phi)+y_ankle),'-','MarkerSize',2,'LineWidth',2); %pulley
    plot([x1 x2],[y1 y2],'-r','MarkerSize',2,'LineWidth',2); %CFL
    plot([x2 x_Ci],[y2 y_Ci],'-r','MarkerSize',2,'LineWidth',2); %Ci
    plot([x2 x3],[y2 y3],'-r','MarkerSize',2,'LineWidth',2); %CFLT
    plot([x3 x_GE_origin],[y3 y_GE_origin],'-r','MarkerSize',2,'LineWidth',2); %GEo
    drawnow limitrate % drawnowを入れるとアニメーションになる
    hold off
    grid on
    xlim([x_lim_neg x_lim_pos]);
    ylim([y_lim_neg y_lim_pos]);
    axis square
    frame= getframe(ax); % アニメーションのフレームをゲットする 
    text_str = ['theta2 = ' num2str(rad2deg(theta2)) ' theta3 = ' num2str(rad2deg(theta3)) ' theta4 = ' num2str(rad2deg(theta4))];
    RGB = insertText(frame.cdata,[1 50],text_str,FontSize=18,TextColor="Black");
    writeVideo(v,RGB);
end
close(v)