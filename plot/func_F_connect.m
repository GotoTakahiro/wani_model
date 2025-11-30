%各関節を分離したときののxy方向の接続力を計算するための関数
function [F_connect] = func_F_connect(q,r,l_link_list,m_list,accel_CoM,F_dist,GRF)

    F_connect=zeros(size(q(:,1),1),8);
    %[coordinates_x, coordinates_y, angle_wire_pulley] = calc_coordinate_for_plot(q,r,l_link_list);
    g=9.81;
    %重力の計算
    F_frame= (m_list(4)+m_list(5))*g;
    F_fem= m_list(6)*g;
    F_tib= m_list(7)*g;
    F_met= m_list(10)*g;
    

    %F_all=[F_CFL_PID;F_Ci_all;F_CFLT_all;F_GEo_all;F_GE_all];
    for i=1:size(q(:,1),1)
        %足関節でのx,y方向の接続力
        F_connect(i,1)=-F_dist(i,9)-GRF(i,3)-GRF(i,5)+m_list(10)*accel_CoM(i,1);
        F_connect(i,2)=-F_dist(i,10)-GRF(i,4)-GRF(i,6)+m_list(10)*accel_CoM(i,2)+F_met;
        

        %膝関節でのx,y方向の接続力
        F_connect(i,3)=F_connect(i,1)+m_list(7)*accel_CoM(i,3);
        F_connect(i,4)=F_connect(i,2)+m_list(7)*accel_CoM(i,4)+F_tib;
        %股関節でのx,y方向の接続力
        F_connect(i,5)=F_connect(i,3)-F_dist(i,1)-F_dist(i,3)+m_list(6)*accel_CoM(i,7)-GRF(i,1);
        F_connect(i,6)=F_connect(i,4)-F_dist(i,2)-F_dist(i,4)+m_list(6)*accel_CoM(i,8)+F_fem-GRF(i,2);
        %
        %frameでのx,y方向の接続力（0になるはず）
        F_connect(i,7)=F_connect(i,5)-F_dist(i,7)+(m_list(4)+m_list(5))*accel_CoM(i,11);
        F_connect(i,8)=F_connect(i,6)-F_dist(i,8)+(m_list(4)+m_list(5))*accel_CoM(i,12)+F_frame;

        %股関節でのx,y方向の接続力
        F_connect(i,9)=-F_dist(i,7)+(m_list(4)+m_list(5))*accel_CoM(i,11)-GRF(i,1);
        F_connect(i,10)=-F_dist(i,8)+(m_list(4)+m_list(5))*accel_CoM(i,12)+F_frame-GRF(i,2);

         %膝関節でのx,y方向の接続力
        F_connect(i,11)=F_connect(i,9)-F_dist(i,1)-F_dist(i,3)+m_list(6)*accel_CoM(i,7);
        F_connect(i,12)=F_connect(i,10)-F_dist(i,2)-F_dist(i,4)+m_list(6)*accel_CoM(i,8)+F_fem;

        %足関節でのx,y方向の接続力
        F_connect(i,13)=F_connect(i,11)+m_list(7)*accel_CoM(i,3);
        F_connect(i,14)=F_connect(i,12)+m_list(7)*accel_CoM(i,4)+F_tib;

        %つま先でのx,y方向の接続力
        F_connect(i,15)=F_connect(i,13)-F_dist(i,9)-GRF(i,3)-GRF(i,5)+m_list(10)*accel_CoM(i,1);
        F_connect(i,16)=F_connect(i,14)-F_dist(i,10)-GRF(i,4)-GRF(i,6)+m_list(10)*accel_CoM(i,2)+F_met;

        
        
    end
    
end