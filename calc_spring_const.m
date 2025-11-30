function [k,c] = calc_spring_const(general_q,l_link_list,l_muscle_list,limit_list,default_wire_k,default_wire_c,k_frame,c_frame)
    global L_GE_present
    r = l_link_list(7);
    k = zeros(7,1);
    c = zeros(8,1);
    
    Coordinate_4th_troch = calc_coordinate_4th_troch(l_link_list,general_q);
    Coordinate_GE_origin = calc_coordinate_GE_origin(l_link_list,general_q);
    Coordinate_M3_to_pulley = calc_coordinate_M3_to_pulley(l_link_list,general_q);
    Coordinate_toe_to_pulley = calc_coordinate_toe_to_pulley(l_link_list,general_q);
    angle_wire_pulley = calc_angle_wire_pulley(l_link_list,general_q);
    cross_wire_pulley = calc_cross_wire_pulley(l_link_list,general_q);

    % プーリからY字分岐点に向かうベクトルと，プーリからつま先に向かうベクトルの外積を求め，正ならベクトルが交差する→回転方向が逆
    if cross_wire_pulley(3) > 0
        angle_wire_pulley = 2*pi-angle_wire_pulley;
    end
    %BP1の座標計算
    x2 = general_q(1) + general_q(10)*sin(general_q(9));
    y2 = general_q(2) - general_q(10)*cos(general_q(9));

    L_Ci_present = sqrt((Coordinate_4th_troch(1)-x2)^2 + (Coordinate_4th_troch(2)-y2)^2);
    L_CFLT_present = sqrt((x2-general_q(3))^2 + (y2-general_q(4))^2);
    L_GEo_present = sqrt((Coordinate_GE_origin(1)-general_q(3))^2 + (Coordinate_GE_origin(2)-general_q(4))^2);

    if cross_wire_pulley(3) < 0
        L_GE_present = sqrt((general_q(3)-Coordinate_M3_to_pulley(1))^2 + (general_q(4)-Coordinate_M3_to_pulley(2))^2) + r*angle_wire_pulley;
    else
        L_GE_present = sqrt((general_q(3)-Coordinate_toe_to_pulley(1))^2 + (general_q(4)-Coordinate_toe_to_pulley(2))^2);
    end

    %L_Ciがたるんでいたら係数0
    if L_Ci_present < l_muscle_list(1)
        k(1) = 0;
        c(1) = 0;
    else
        k(1) = default_wire_k;
        c(1) = default_wire_c;
    end
    %L_CFLTがたるんでいたら係数0
    if L_CFLT_present < l_muscle_list(2)
        k(2) = 0;
        c(2) = 0;
    else
        k(2) = default_wire_k;
        c(2) = default_wire_c;
    end
    %L_GEoがたるんでいたら係数0
    if L_GEo_present < l_muscle_list(3)
        k(3) = 0;
        c(3) = 0;
    else
        k(3) = default_wire_k;
        c(3) = default_wire_c;
    end
    %L_GEがたるんでいたら係数0
    if L_GE_present < l_muscle_list(4)
        k(4) = 0;
        c(4) = 0;
    else
        k(4) = default_wire_k;
        c(4) = default_wire_c;
    end

    % 足関節の背屈の制限
    if general_q(8) > limit_list(1)
        k(5) = k_frame;    %default_wire_k;
        %c(5) = c_frame/50;
    else
        k(5) = 0;
        %c(5) = 0;
    end
    %足関節の底屈の制限
    if general_q(8) < limit_list(2)
        k(6) = k_frame;    %default_wire_k;
        %c(6) = c_frame/50;
    else
        k(6) = 0;
        %c(6) = 0;
    end

    % 足関節の粘性は常に入れる
    %底屈の粘性は，底屈の条件を満たすときだけ生じるように変更
    c(5) = c_frame/50;      %ankle
    c(6) = c_frame/50;      %ankle

    % フレームの角度が常に一定（水平）に保つための係数，角度theta1に依存
    k(7) = k_frame;    %frame
    c(7) = c_frame/20; %frame
    %CFLの粘性係数c_CFLは0 CFLの粘弾性項は粘弾性エネルギーからではなく，l_CFLを一定に保つための力の方に含まれている
    c(8) = 0; %c_CFL
end