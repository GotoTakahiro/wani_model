%各リンクの並進方向の加速度を計算するための式
function [accel_CoM] = func_accel(phi,dphi,ddphi,m_list,l_link_list)

    
    L_femur=l_link_list(1);
    L_tibia=l_link_list(2);
    L_met=l_link_list(3);
    %L_4th_troch=l_link_list(4);
    %L_GE_origin=l_link_list(5);
    L_frame=l_link_list(6);

    M_hip=m_list(4);
    M_frame=m_list(5);
    %m_list(6) = M_fem;
    %m_list(7) = M_tib;
    M_met_pulley=m_list(8);
    M_met_rod=m_list(9);


    %metatarsalの端部から重心までの長さ
    L_met_CoM = (L_met*M_met_rod/2)/(M_met_pulley + M_met_rod);
    %hipの端部から重心までの長さ
    L_frame_CoM = L_frame-(L_frame*M_frame/2 + M_hip*L_frame)/(M_frame + M_hip);
    %tibiaの重心までの長さ
    L_tibia_CoM = L_tibia/2;
    %femurの重心までの長さ
    L_femur_CoM = L_femur/2;
    accel_CoM=zeros(size(phi(:,1),1),14);

    for in=1:size(phi(:,1),1)
        %ankleからmetatarsalの重心までのベクトル
        r_met_CoM=[L_met_CoM*sin(phi(in,4)) -L_met_CoM*cos(phi(in,4))];
        %ankleの重心の並進加速度
        accel_met_CoM=[-ddphi(in,4)*r_met_CoM(1,2) ddphi(in,4)*r_met_CoM(1,1)]+[-dphi(in,4)^2*r_met_CoM(1,1) -dphi(in,4)^2*r_met_CoM(1,2)];
        
        %ankleから足先までのベクトル
        r_foot_CoM=[L_met*sin(phi(in,4)) -L_met*cos(phi(in,4))];
        %ankleの重心の並進加速度
        accel_foot_CoM=[-ddphi(in,4)*r_foot_CoM(1,2) ddphi(in,4)*r_foot_CoM(1,1)]+[-dphi(in,4)^2*r_foot_CoM(1,1) -dphi(in,4)^2*r_foot_CoM(1,2)];
    
        %ankleからtibiaの重心までのベクトル
        r_tibia_CoM=[L_tibia_CoM*sin(-phi(in,3)) L_tibia_CoM*cos(phi(in,3))];
        %tibiaの重心の並進加速度
        accel_tibia_CoM=[-ddphi(in,3)*r_tibia_CoM(1,2) ddphi(in,3)*r_tibia_CoM(1,1)]+[-dphi(in,3)^2*r_tibia_CoM(1,1) -dphi(in,3)^2*r_tibia_CoM(1,2)];
        
        %kneeの位置までのベクトル
        r_knee_CoM=[L_tibia*sin(-phi(in,3)) L_tibia*cos(phi(in,3))];
        %kneeの重心の並進加速度
        accel_knee_CoM=[-ddphi(in,3)*r_knee_CoM(1,2) ddphi(in,3)*r_knee_CoM(1,1)]+[-dphi(in,3)^2*r_knee_CoM(1,1) -dphi(in,3)^2*r_knee_CoM(1,2)];
        
        %femurの位置までのベクトル
        r_femur_CoM=[L_femur_CoM*sin(-phi(in,2)) L_femur_CoM*cos(phi(in,2))];
        %femurの重心の並進加速度
        accel_femur_CoM=accel_knee_CoM+[-ddphi(in,2)*r_femur_CoM(1,2) ddphi(in,2)*r_femur_CoM(1,1)]+[-dphi(in,2)^2*r_femur_CoM(1,1) -dphi(in,2)^2*r_femur_CoM(1,2)];
        
        %hipの位置までのベクトル
        r_hip_CoM=[L_femur*sin(-phi(in,2)) L_femur*cos(phi(in,2))];
        %hipの重心の並進加速度
        accel_hip_CoM=accel_knee_CoM+[-ddphi(in,2)*r_hip_CoM(1,2) ddphi(in,2)*r_hip_CoM(1,1)]+[-dphi(in,2)^2*r_hip_CoM(1,1) -dphi(in,2)^2*r_hip_CoM(1,2)];
    
        %frameの重心位置までのベクトル
        r_frame_CoM=[L_frame_CoM*sin(-phi(in,1)) L_frame_CoM*cos(phi(in,1))];
        %frameの重心の並進加速度
        accel_frame_CoM=accel_hip_CoM+[-ddphi(in,1)*r_frame_CoM(1,2) ddphi(in,1)*r_frame_CoM(1,1)]+[-dphi(in,1)^2*r_frame_CoM(1,1) -dphi(in,1)^2*r_frame_CoM(1,2)];
    
        %frameの先端位置までのベクトル
        r_tail_CoM=[L_frame*sin(-phi(in,1)) L_frame*cos(phi(in,1))];
        %frameの重心の並進加速度
        accel_tail_CoM=accel_hip_CoM+[-ddphi(in,1)*r_tail_CoM(1,2) ddphi(in,1)*r_tail_CoM(1,1)]+[-dphi(in,1)^2*r_tail_CoM(1,1) -dphi(in,1)^2*r_tail_CoM(1,2)];
    
    
        %計算結果の格納
        accel_CoM(in,1)=accel_met_CoM(1,1);
        accel_CoM(in,2)=accel_met_CoM(1,2);
        accel_CoM(in,3)=accel_tibia_CoM(1,1);
        accel_CoM(in,4)=accel_tibia_CoM(1,2);
        accel_CoM(in,5)=accel_knee_CoM(1,1);
        accel_CoM(in,6)=accel_knee_CoM(1,2);
        accel_CoM(in,7)=accel_femur_CoM(1,1);
        accel_CoM(in,8)=accel_femur_CoM(1,2);
        accel_CoM(in,9)=accel_hip_CoM(1,1);
        accel_CoM(in,10)=accel_hip_CoM(1,2);
        accel_CoM(in,11)=accel_frame_CoM(1,1);
        accel_CoM(in,12)=accel_frame_CoM(1,2);
        accel_CoM(in,11)=accel_tail_CoM(1,1);
        accel_CoM(in,12)=accel_tail_CoM(1,2);
        accel_CoM(in,13)=accel_foot_CoM(1,1);
        accel_CoM(in,14)=accel_foot_CoM(1,2);
    end
    
end