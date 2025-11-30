function [coordinates_x, coordinates_y, angle_wire_pulley] = calc_coordinate_for_plot(q,r,L_link)
    coordinates_x=zeros(size(q,1),12);
    coordinates_y=zeros(size(q,1),12);
    angle_wire_pulley=zeros(size(q,1),1);

    for i = 1:size(q,1)
        general_q = q(i,1:10).';
        %CFL_originの座標
        coordinates_x(i,1) = q(i,1); 
        coordinates_y(i,1) = q(i,2);
        %hipの座標
        coordinates_x(i,2) = coordinates_x(i,1) + L_link(6)*sin(q(i,5)); 
        coordinates_y(i,2) = coordinates_y(i,1) - L_link(6)*cos(q(i,5));
        %kneeの座標
        coordinates_x(i,3) = coordinates_x(i,2) + L_link(1)*sin(q(i,5)+q(i,6)); 
        coordinates_y(i,3) = coordinates_y(i,2) - L_link(1)*cos(q(i,5)+q(i,6));
        %ankleの座標
        coordinates_x(i,4) = coordinates_x(i,3) + L_link(2)*sin(q(i,5)+q(i,6)+q(i,7)); 
        coordinates_y(i,4) = coordinates_y(i,3) - L_link(2)*cos(q(i,5)+q(i,6)+q(i,7));
        %toeの座標
        coordinates_x(i,5) = coordinates_x(i,4) + L_link(3)*sin(q(i,5)+q(i,6)+q(i,7)+q(i,8));
        coordinates_y(i,5) = coordinates_y(i,4) - L_link(3)*cos(q(i,5)+q(i,6)+q(i,7)+q(i,8));
        %CFLT branching point（BP1）の座標
        coordinates_x(i,6) = q(i,1) + q(i,10)*sin(q(i,9)); 
        coordinates_y(i,6) = q(i,2) - q(i,10)*cos(q(i,9));
        %Y shaped（BP2）の座標
        coordinates_x(i,7) = q(i,3); 
        coordinates_y(i,7) = q(i,4);
        %M3からpulleyへ接線を引いたときの接点
        temp = calc_coordinate_M3_to_pulley(L_link,general_q);
        coordinates_x(i,8) = temp(1,1); 
        coordinates_y(i,8) = temp(2,1);
        %描画のときに恐らく使っていない
        coordinates_x(i,9) = q(i,7); 
        coordinates_y(i,9) = q(i,8);
        %4th trochanter（Ciの付着位置）の座標
        coordinates_x(i,10) = coordinates_x(i,2) + L_link(4)*sin(q(i,5)+q(i,6)); 
        coordinates_y(i,10) = coordinates_y(i,2) - L_link(4)*cos(q(i,5)+q(i,6));
        %GE_origin（GEoの付着位置）の座標  
        coordinates_x(i,11) = coordinates_x(i,2) + L_link(5)*sin(q(i,5)+q(i,6)); 
        coordinates_y(i,11) = coordinates_y(i,2) - L_link(5)*cos(q(i,5)+q(i,6));
        %toeからpulleyへ接線を引いたときの接点
        temp2 = calc_coordinate_toe_to_pulley(L_link,general_q);
        coordinates_x(i,12) = temp2(1,1); 
        coordinates_y(i,12) = temp2(2,1);

        %プーリーにかかっているワイヤの角度
        angle_wire_pulley(i) = calc_angle_wire_pulley(L_link,general_q);
    end