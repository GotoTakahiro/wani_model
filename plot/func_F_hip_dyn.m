%一般化力からhipの質点に働く動的な力を逆算するための関数
function [F] = func_F_hip_dyn(q,dq,ddq,m_list,l_link_list)

    
    F = zeros(size(q,1),2);
    L_frame=l_link_list(6);

    for i=1:size(q,1)

        M = Inertial_matrix(m_list,l_link_list,transpose(q(i,1:10)));
        C = Coriolis_matrix(m_list,l_link_list,transpose(q(i,1:10)),transpose(dq(i,1:10))); %コリオリ力・遠心力の計算
        T_dyn=M*transpose(ddq(i,1:10))+C;
        J_p =[1, 0, 0, 0, L_frame*cos(q(i,5)), 0, 0, 0, 0, 0;0, 1, 0, 0, L_frame*sin(q(i,5)), 0, 0, 0, 0, 0];
        Mp = J_p * J_p.'; %正方行列の計算
        if rcond(Mp) < 1e-12
            % 正則化または擬似逆
            Fp = pinv(J_p.') * T_dyn;   % ここは (J J^T)^(-1) J tau の代替：最小ノルム解
            % あるいは
            % Fp = (M + 1e-8*eye(2)) \ (J * tau);
              warning('1')
        else
            Fp = Mp \ (J_p * T_dyn);     % (J J^T)^{-1} * (J * tau)
        end
        F(i,1)=Fp(1,1);
        F(i,2)=Fp(2,1);
    end

    
end