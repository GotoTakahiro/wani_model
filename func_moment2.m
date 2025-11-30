%尻尾周りの地面反力と壁のモーメントの釣り合いを計算
%
function [moment_CoM] = func_moment2(q,m_list,l_link_list,GRF,r,F)

    % 座標を計算
    %1:hip, 2:4th trochanter, 3:GE origin, 4:knee, 5:ankle, 6:toem 7:CFTLT branch, 8:Y shaped branch
    [coordinates_x, coordinates_y, angle_wire_pulley] = calc_coordinate_for_plot(q,r,l_link_list); 
    COM_link=zeros(size(q(:,1),1),2);
    moment_CoM=zeros(size(q(:,1),1),16);
    for j = 1:size(q,1)
        %qの1~10までを抜き出し
        general_q = q(j,1:10).';
        %重心座標を求める
        COM_link(j,:) = (calc_COM_link(m_list,l_link_list,general_q))';
    end

    
    L_femur=l_link_list(1);
    L_tibia=l_link_list(2);
    L_met=l_link_list(3);
    %L_4th_troch=l_link_list(4);
    %L_GE_origin=l_link_list(5);
    L_frame=l_link_list(6);

    M_hip=m_list(4);
    M_frame=m_list(5);
    M_fem=m_list(6);
    M_tib=m_list(7);
    M_met_pulley=m_list(8);
    M_met_rod=m_list(9);
    %総質量
    M_all=m_list(1)+m_list(2)+m_list(3)+m_list(4)+m_list(5)+m_list(6)+m_list(7)+m_list(8)+m_list(9);
    g=-9.81;

    %metatarsalの端部から重心までの長さ
    L_met_CoM = (L_met*M_met_rod/2)/(M_met_pulley + M_met_rod);
    %hipの端部から重心までの長さ
    L_frame_CoM = L_frame-(L_frame*M_frame/2 + M_hip*L_frame)/(M_frame + M_hip);
    %tibiaの重心までの長さ
    L_tibia_CoM = L_tibia/2;
    %femurの重心までの長さ
    L_femur_CoM = L_femur/2;
    %各張力のxy方向成分を計算するための式
    F_dist = func_F_dist(q,r,l_link_list,F);

    %data_accel_GRF(12,13):hip (14,15):heel (16,17):toe
    %重心周りのモーメント計算
    for i=1:size(q(:,1),1)
        %coordinates(i,4):ankleの座標
        %coordinates(i,1):CFL_originの座標
        %CFL起始点から踵までのベクトル
        vector1=[(coordinates_x(i,4)-coordinates_x(i,1)) (coordinates_y(i,4)-r-coordinates_y(i,1))];
        %踵の床反力のx方向成分の力
        moment_CoM(i,1)=GRF(i,14)*(-vector1(1,2));
        %踵の床反力のy方向成分の力
        moment_CoM(i,2)=GRF(i,15)*vector1(1,1);
        
        %coordinates(i,5):toeの座標
        %CFL起始点からつま先までのベクトル
        vector2=[(coordinates_x(i,5)-coordinates_x(i,1)) (coordinates_y(i,5)-coordinates_y(i,1))];
        %つま先の床反力のx方向成分の力
        moment_CoM(i,3)=GRF(i,16)*(-vector2(1,2));
        %つま先の床反力のy方向成分の力
        moment_CoM(i,4)=GRF(i,17)*vector2(1,1);

        %coordinates(i,2):hipの座標
        %CFL起始点から股関節までのベクトル
        vector3=[(coordinates_x(i,2)-coordinates_x(i,1)) (coordinates_y(i,2)-coordinates_y(i,1))];
        %股関節の床反力のx方向成分の力
        moment_CoM(i,5)=GRF(i,12)*(-vector3(1,2));
        %股関節の床反力のy方向成分の力
        moment_CoM(i,6)=GRF(i,13)*vector3(1,1);

        
        %CFL起始点から重心までのベクトル
        vector4=[(COM_link(i,1)-coordinates_x(i,1)) (COM_link(i,2)-coordinates_y(i,1))];
        %CFLのx方向成分の力
        moment_CoM(i,7)=0*(-vector4(1,2));
        %CFLのy方向成分の力
        moment_CoM(i,8)=M_all*g*vector4(1,1);

        %coordinates(i,10):4th trochanter（Ciの付着位置）の座標
        %CFL起始点からCi起始点までのベクトル
        vector5=[(coordinates_x(i,10)-coordinates_x(i,1)) (coordinates_y(i,10)-coordinates_y(i,1))];
        %Ciのx方向成分の力
        moment_CoM(i,9)=F_dist(i,1)*(-vector5(1,2));
        %Ciのy方向成分の力
        moment_CoM(i,10)=F_dist(i,2)*vector5(1,1);

        %coordinates(i,11):GE_origin（GEoの付着位置）の座標
        %CFL起始点からGEo起始点までのベクトル
        vector6=[(coordinates_x(i,11)-coordinates_x(i,1)) (coordinates_y(i,11)-coordinates_y(i,1))];
        %GEoのx方向成分の力
        moment_CoM(i,11)=F_dist(i,3)*(-vector6(1,2));
        %GEoのy方向成分の力
        moment_CoM(i,12)=F_dist(i,4)*vector6(1,1);

        %coordinates(i,12):toeからpulleyへ接線を引いたときの接点の座標
        %CFL起始点からGE停止点までのベクトル
        vector7=[(coordinates_x(i,12)-coordinates_x(i,1)) (coordinates_y(i,12)-coordinates_y(i,1))];
        %GEのx方向成分の力
        moment_CoM(i,13)=F_dist(i,5)*(-vector7(1,2));
        %GEのy方向成分の力
        moment_CoM(i,14)=F_dist(i,6)*vector7(1,1);

        %coordinates(i,8):M3からpulleyへ接線を引いたときの接点の座標
        %CFL起始点からGE接点までのベクトル
        vector8=[(coordinates_x(i,8)-coordinates_x(i,1)) (coordinates_y(i,8)-coordinates_y(i,1))];
        %GEのx方向成分の力（接点で作用すると仮定）
        moment_CoM(i,15)=F_dist(i,9)*(-vector8(1,2));
        %GEのy方向成分の力（接点で作用すると仮定）
        moment_CoM(i,16)=F_dist(i,10)*vector8(1,1);

        
        
        
       
    end
    
end