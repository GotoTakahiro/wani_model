%各張力のxy方向成分を計算するための式
function [F_dist] = func_F_dist(q,r,l_link_list,F)

    F_dist=zeros(size(q(:,1),1),12);
    [coordinates_x, coordinates_y, angle_wire_pulley] = calc_coordinate_for_plot(q,r,l_link_list);
    %F_all=[F_CFL_PID;F_Ci_all;F_CFLT_all;F_GEo_all;F_GE_all];
    for i=1:size(q(:,1),1)
        %coordinates(i,6):CFLT branching point（BP1）の座標
        %coordinates(i,10):4th trochanter（Ciの付着位置）の座標
        %Ciの停止位置を基準とした単位方向ベクトルを求める
        u_Ci=[(coordinates_x(i,6)-coordinates_x(i,10)) (coordinates_y(i,6)-coordinates_y(i,10))]/sqrt((coordinates_x(i,6)-coordinates_x(i,10))^2+(coordinates_y(i,6)-coordinates_y(i,10))^2);
        F_dist(i,1)=F(i,2)*u_Ci(1,1);
        F_dist(i,2)=F(i,2)*u_Ci(1,2);
        %coordinates(i,7):Y shaped（BP2）の座標
        %coordinates(i,11):GE_origin（GEoの付着位置）の座標
        %GEoの停止位置を基準とした単位方向ベクトルを求める
        u_GEo=[(coordinates_x(i,7)-coordinates_x(i,11)) (coordinates_y(i,7)-coordinates_y(i,11))]/sqrt((coordinates_x(i,7)-coordinates_x(i,11))^2+(coordinates_y(i,7)-coordinates_y(i,11))^2);
        F_dist(i,3)=F(i,4)*u_GEo(1,1);
        F_dist(i,4)=F(i,4)*u_GEo(1,2);
        %coordinates(i,12):toeからpulleyへ接線を引いたときの接点の座標
        %coordinates(i,5):toeの座標
        %GEの停止位置を基準とした単位方向ベクトルを求める
        u_GE=[(coordinates_x(i,12)-coordinates_x(i,5)) (coordinates_y(i,12)-coordinates_y(i,5))]/sqrt((coordinates_x(i,12)-coordinates_x(i,5))^2+(coordinates_y(i,12)-coordinates_y(i,5))^2);
        F_dist(i,5)=F(i,5)*u_GE(1,1);
        F_dist(i,6)=F(i,5)*u_GE(1,2);
        %coordinates(i,6):CFLT branching point（BP1）の座標
        %coordinates(i,1):CFL_originの座標
        %CFLの停止位置を基準とした単位方向ベクトルを求める
        u_CFL=[(coordinates_x(i,6)-coordinates_x(i,1)) (coordinates_y(i,6)-coordinates_y(i,1))]/sqrt((coordinates_x(i,6)-coordinates_x(i,1))^2+(coordinates_y(i,6)-coordinates_y(i,1))^2);
        F_dist(i,7)=-F(i,1)*u_CFL(1,1);
        F_dist(i,8)=-F(i,1)*u_CFL(1,2);
        %coordinates(i,7):Y shaped（BP2）の座標
        %coordinates(i,8):M3からpulleyへ接線を引いたときの接点の座標
        %GEの接点位置を基準とした単位方向ベクトルを求める
        u_GE_test=[(coordinates_x(i,7)-coordinates_x(i,8)) (coordinates_y(i,7)-coordinates_y(i,8))]/sqrt((coordinates_x(i,7)-coordinates_x(i,8))^2+(coordinates_y(i,7)-coordinates_y(i,8))^2);
        F_dist(i,9)=F(i,5)*u_GE_test(1,1);
        F_dist(i,10)=F(i,5)*u_GE_test(1,2);
        %GEの張力によって生じるプーリー中心からの力
        F_dist(i,11)=F_dist(i,5)+F_dist(i,9);
        F_dist(i,12)=F_dist(i,6)+F_dist(i,10);


    end
    
end