%一般化力からankleの質点に働く力を逆算するための関数
function [F] = func_F_ankle_dyn(q,dq,ddq,m_list,l_link_list)

    F = zeros(size(q,1),2);

    L_fem=l_link_list(1);
    L_tib=l_link_list(2);
    %L_4th_troch=l_link_list(4);
    %L_GE_origin=l_link_list(5);
    L_frame=l_link_list(6);
    %r=l_link_list(7);

    
    J_p=zeros(2,10);
    for i=1:size(q,1)
        theta1=q(i,5);
        theta2=q(i,6);
        theta3=q(i,7);

        J_p(1,1)=1;
        J_p(2,2)=1;
        J_p(1,5)=L_fem*cos(theta1 + theta2) + L_frame*cos(theta1) + L_tib*cos(theta1 + theta2 + theta3);
        J_p(1,6)=L_fem*cos(theta1 + theta2) + L_tib*cos(theta1 + theta2 + theta3);
        J_p(1,7)=L_tib*cos(theta1 + theta2 + theta3);
        J_p(2,5)=L_fem*sin(theta1 + theta2) + L_frame*sin(theta1) + L_tib*sin(theta1 + theta2 + theta3);
        J_p(2,6)=L_fem*sin(theta1 + theta2) + L_tib*sin(theta1 + theta2 + theta3);
        J_p(2,7)=L_tib*sin(theta1 + theta2 + theta3);
        
        M = Inertial_matrix(m_list,l_link_list,transpose(q(i,1:10)));
        C = Coriolis_matrix(m_list,l_link_list,transpose(q(i,1:10)),transpose(dq(i,1:10))); %コリオリ力・遠心力の計算
        T_dyn=M*transpose(ddq(i,1:10))+C;


        Mp = J_p * J_p.';
        if rcond(Mp) < 1e-12
            % 正則化または擬似逆
            Fp = pinv(J_p.') * T_dyn;   % ここは (J J^T)^(-1) J tau の代替：最小ノルム解
            % あるいは
            % Fp = (M + 1e-8*eye(2)) \ (J * tau);
        else
            Fp = Mp \ (J_p * T_dyn);     % (J J^T)^{-1} * (J * tau)
        end
        F(i,1)=Fp(1,1);
        F(i,2)=Fp(2,1);
    end

    
end