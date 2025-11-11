%一般化力からCOMの質点に働く力を逆算するための関数
function [F] = func_F_COM(q,m_list,l_link_list,Q)

    Fp = zeros(2,1);
    F = zeros(size(q,1),2);

    L_fem=l_link_list(1);
    L_tib=l_link_list(2);
    L_met=l_link_list(3);
    %L_4th_troch=l_link_list(4);
    %L_GE_origin=l_link_list(5);
    L_frame=l_link_list(6);
    %r=l_link_list(7);
    M1=m_list(1);
    M2=m_list(2);
    M3=m_list(3);
    M_hip=m_list(4);
    M_frame=m_list(5);
    M_fem=m_list(6);
    M_tib=m_list(7);
    M_met_pulley=0.3;%m_list(8);
    M_met_rod=0.1;%m_list(9);
    M_met=0.4;%m_list(10);
    
    J_p=zeros(2,10);
    for i=1:size(q,1)
        theta1=q(i,5);
        theta2=q(i,6);
        theta3=q(i,7);
        theta4=q(i,8);
        theta_CFL=q(i,9);
        l_CFL=q(i,10);
        J_p(1,1)=(M1 + M2 + M_fem + M_frame + M_met + M_tib)/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(1,3)=M3/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(1,5)=(M_met*(L_fem*cos(theta1 + theta2) + L_frame*cos(theta1) + L_tib*cos(theta1 + theta2 + theta3) + (L_met*M_met_rod*cos(theta1 + theta2 + theta3 + theta4))/(2*(M_met_rod + M_met_pulley))) + M_tib*(L_fem*cos(theta1 + theta2) + L_frame*cos(theta1) + (L_tib*cos(theta1 + theta2 + theta3))/2) + M_fem*((L_fem*cos(theta1 + theta2))/2 + L_frame*cos(theta1)) + (M_frame*cos(theta1)*((L_frame*M_frame)/2 + L_frame*M_hip))/(M_frame + M_hip))/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(1,6)=(M_met*(L_fem*cos(theta1 + theta2) + L_tib*cos(theta1 + theta2 + theta3) + (L_met*M_met_rod*cos(theta1 + theta2 + theta3 + theta4))/(2*(M_met_rod + M_met_pulley))) + M_tib*(L_fem*cos(theta1 + theta2) + (L_tib*cos(theta1 + theta2 + theta3))/2) + (L_fem*M_fem*cos(theta1 + theta2))/2)/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(1,7)=(M_met*(L_tib*cos(theta1 + theta2 + theta3) + (L_met*M_met_rod*cos(theta1 + theta2 + theta3 + theta4))/(2*(M_met_rod + M_met_pulley))) + (L_tib*M_tib*cos(theta1 + theta2 + theta3))/2)/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(1,8)=(L_met*M_met*M_met_rod*cos(theta1 + theta2 + theta3 + theta4))/(2*(M_met_rod + M_met_pulley)*(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib));
        J_p(1,9)=(M2*l_CFL*cos(theta_CFL))/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(1,10)=(M2*sin(theta_CFL))/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(2,2)=(M1 + M2 + M_fem + M_frame + M_met + M_tib)/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(2,4)=M3/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(2,5)=(M_fem*((L_fem*sin(theta1 + theta2))/2 + L_frame*sin(theta1)) + M_met*(L_fem*sin(theta1 + theta2) + L_frame*sin(theta1) + L_tib*sin(theta1 + theta2 + theta3) + (L_met*M_met_rod*sin(theta1 + theta2 + theta3 + theta4))/(2*(M_met_rod + M_met_pulley))) + M_tib*(L_fem*sin(theta1 + theta2) + L_frame*sin(theta1) + (L_tib*sin(theta1 + theta2 + theta3))/2) + (M_frame*sin(theta1)*((L_frame*M_frame)/2 + L_frame*M_hip))/(M_frame + M_hip))/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(2,6)=(M_met*(L_fem*sin(theta1 + theta2) + L_tib*sin(theta1 + theta2 + theta3) + (L_met*M_met_rod*sin(theta1 + theta2 + theta3 + theta4))/(2*(M_met_rod + M_met_pulley))) + M_tib*(L_fem*sin(theta1 + theta2) + (L_tib*sin(theta1 + theta2 + theta3))/2) + (L_fem*M_fem*sin(theta1 + theta2))/2)/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(2,7)=(M_met*(L_tib*sin(theta1 + theta2 + theta3) + (L_met*M_met_rod*sin(theta1 + theta2 + theta3 + theta4))/(2*(M_met_rod + M_met_pulley))) + (L_tib*M_tib*sin(theta1 + theta2 + theta3))/2)/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(2,8)=(L_met*M_met*M_met_rod*sin(theta1 + theta2 + theta3 + theta4))/(2*(M_met_rod + M_met_pulley)*(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib));
        J_p(2,9)=(M2*l_CFL*sin(theta_CFL))/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);
        J_p(2,10)=-(M2*cos(theta_CFL))/(M1 + M2 + M3 + M_fem + M_frame + M_met + M_tib);



        M = J_p * J_p.';
        if rcond(M) < 1e-12
            % 正則化または擬似逆
            Fp = pinv(J_p.') * transpose(Q(i,:));   % ここは (J J^T)^(-1) J tau の代替：最小ノルム解
            % あるいは
            % Fp = (M + 1e-8*eye(2)) \ (J * tau);
        else
            Fp = M \ (J_p * transpose(Q(i,:)));     % (J J^T)^{-1} * (J * tau)
        end
        F(i,1)=Fp(1,1);
        F(i,2)=Fp(2,1);
    end

    
end