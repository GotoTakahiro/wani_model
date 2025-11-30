%各張力のxy方向成分を計算するための式
function [musclelength] = func_musclelength(q,r,l_link_list)
    
    musclelength=zeros(size(q(:,1),1),5);
    [coordinates_x, coordinates_y, angle_wire_pulley] = calc_coordinate_for_plot(q,r,l_link_list);
    %F_all=[F_CFL_PID;F_Ci_all;F_CFLT_all;F_GEo_all;F_GE_all];
    for i=1:size(q(:,1),1)

        %CFLの長さを計算
        musclelength(i,1)=sqrt((coordinates_x(i,6)-coordinates_x(i,1))^2+(coordinates_y(i,6)-coordinates_y(i,1))^2);
        
        %coordinates(i,6):CFLT branching point（BP1）の座標
        %coordinates(i,10):4th trochanter（Ciの付着位置）の座標
        %Ciの長さを計算
        musclelength(i,2)=sqrt((coordinates_x(i,6)-coordinates_x(i,10))^2+(coordinates_y(i,6)-coordinates_y(i,10))^2);
        
        %coordinates(i,7):Y shaped（BP2）の座標
        %coordinates(i,11):GE_origin（GEoの付着位置）の座標
        %GEoの長さを計算
        musclelength(i,3)=sqrt((coordinates_x(i,7)-coordinates_x(i,11))^2+(coordinates_y(i,7)-coordinates_y(i,11))^2);
        %CFLTの長さを計算
        musclelength(i,4)=sqrt((coordinates_x(i,7)-coordinates_x(i,6))^2+(coordinates_y(i,7)-coordinates_y(i,6))^2);

        %coordinates(i,7):Y shaped（BP2）の座標
        %coordinates(i,8):M3からpulleyへ接線を引いたときの接点の座標
        %GEの接点位置を基準とした単位方向ベクトルを求める
        musclelength(i,5)=sqrt((coordinates_x(i,7)-coordinates_x(i,8))^2+(coordinates_y(i,7)-coordinates_y(i,8))^2);

        %coordinates(i,12):toeからpulleyへ接線を引いたときの接点の座標
        %coordinates(i,5):toeの座標
        %GEの停止位置を基準とした単位方向ベクトルを求める
        % u_GE=[(coordinates_x(i,12)-coordinates_x(i,5)) (coordinates_y(i,12)-coordinates_y(i,5))]/sqrt((coordinates_x(i,12)-coordinates_x(i,5))^2+(coordinates_y(i,12)-coordinates_y(i,5))^2);
        % musclelength(i,5)=F(i,5)*u_GE(1,1);
        % musclelength(i,6)=F(i,5)*u_GE(1,2);

    end
    
end