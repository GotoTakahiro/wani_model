%一般化力からhipの質点に働く力を逆算するための関数
function [F] = func_F_hip(q,l_link_list,Q)

    Fp = zeros(2,1);
    F = zeros(size(q,1),2);
    L_frame=l_link_list(6);
    for i=1:size(q,1)
        J_p =[1, 0, 0, 0, L_frame*cos(q(i,5)), 0, 0, 0, 0, 0;0, 1, 0, 0, L_frame*sin(q(i,5)), 0, 0, 0, 0, 0];
        Mp = J_p * J_p.';
        if rcond(Mp) < 1e-12
            % 正則化または擬似逆
            Fp = pinv(J_p.') * transpose(Q(i,:));   % ここは (J J^T)^(-1) J tau の代替：最小ノルム解
            % あるいは
            % Fp = (M + 1e-8*eye(2)) \ (J * tau);
              warning('1')
        else
            Fp = Mp \ (J_p * transpose(Q(i,:)));     % (J J^T)^{-1} * (J * tau)
        end
        F(i,1)=Fp(1,1);
        F(i,2)=Fp(2,1);
    end

    
end