function [F] = func_GRF2(q,k_ground,c_ground,mu,y_fixed,y_hip,y_heel,y_toe,dx_heel_vec,dx_toe_vec,dx_hip_vec)

    F = zeros(8,1);

    %股関節が一定以上沈み込まないように支持する力
    if y_hip < y_fixed(1)
        F(2) = -k_ground*(y_hip-y_fixed(1)) - c_ground*dx_hip_vec(2); %y方向
        %支持力によって生じるレール水平方向の摩擦力
        %ロボットでは進行方向に滑るようにベアリングがついていたので，それを再現して股関節に関しては摩擦係数は0.1倍
        F(1) = mu*0.1*F(2)*(-dx_hip_vec(1)/abs(dx_hip_vec(1))); %x方向
    end

    %踵で受ける地面からの反力
    if y_heel < y_fixed(2)
        F(4) = -k_ground*(y_heel-y_fixed(2)) - c_ground*dx_heel_vec(2);
        F(3) = mu*F(4)*(-dx_heel_vec(1)/abs(dx_heel_vec(1)));
    end
    
    %つま先で受ける地面からの反力
    if y_toe < y_fixed(2)
        F(6) = -k_ground*(y_toe-y_fixed(2)) - c_ground*dx_toe_vec(2);
        F(5) = mu*F(6)*(-dx_toe_vec(1)/abs(dx_toe_vec(1)));
    end
end