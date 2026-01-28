%各関節を分離したときののxy方向の接続力を計算するための関数　上体起点
function [M_connect] = func_M_connect(q,r,l_link_list,m_list,accel_CoM,F_dist,GRF,data_F_all,ddphi,data_T_all,F_connect,data_accel)

    M_connect=zeros(size(q(:,1),1),37);
    [coordinates_x, coordinates_y, angle_wire_pulley] = calc_coordinate_for_plot(q,r,l_link_list);
    g=9.81;
    %重力の計算
    F_frame= (m_list(4)+m_list(5))*g;
    F_fem= m_list(6)*g;
    F_tib= m_list(7)*g;
    F_met= m_list(10)*g;
    %慣性モーメントの計算
    J_list=calc_J_list(m_list, l_link_list); % [J_frame; J_fem; J_tib; J_met;]
    L_met_CoM = (l_link_list(3)*m_list(9)/2)/(m_list(8) + m_list(9));
    L_frame_CoM = (l_link_list(6)*m_list(5)/2 + l_link_list(6)*m_list(4))/(m_list(5) + m_list(4));
    J_list(1)=J_list(1)+(m_list(4)+m_list(5))*L_met_CoM^2;
    J_list(2)=J_list(2)+m_list(6)*l_link_list(1)^2/4;
    J_list(3)=J_list(3)+m_list(3)*l_link_list(2)^2/4;
    J_list(4)=J_list(4)+m_list(10)*L_frame_CoM^2;
    
    %重心の座標計算
    [COM_x,COM_y] = calc_coordinate_COM(q,l_link_list,m_list);
    

    
    %GRF=[hip_x,hip_y,heel_x,heel_y,toe_x,toe_y]
    %F_all=[F_CFL_PID;F_Ci_all;F_CFLT_all;F_GEo_all;F_GE_all];
    for i=1:size(q(:,1),1)
        COM=calc_COM(m_list, l_link_list, transpose(q(i,1:10)));
        %関節周りの慣性モーメントを計算
        M = Inertial_matrix(m_list,l_link_list,transpose(q(i,1:10)));
        Int=M*transpose(data_accel(i,2:11));
        %足関節でのx,y方向の接続力
        M_connect(i,1)=J_list(4)*ddphi(i,4)+m_list(10)*(accel_CoM(i,1)*(COM_y(i,4)-coordinates_y(i,4))+accel_CoM(i,2)*(COM_x(i,4)-coordinates_x(i,4)))...
                        -(GRF(i,3)+GRF(i,5))*r-GRF(i,6)*(coordinates_x(i,5)-coordinates_x(i,4))...
                        +data_F_all(i,5)*r+F_met*(COM_x(i,4)-coordinates_x(i,4))-data_T_all(i,93);

        M_connect(i,5)=(GRF(i,3)+GRF(i,5))*r; %地面反力の水平方向成分
        M_connect(i,6)=GRF(i,6)*(coordinates_x(i,5)-coordinates_x(i,4)); %地面反力の垂直方向成分
        M_connect(i,7)=-data_F_all(i,5)*r; %F_GE方向成分
        M_connect(i,8)=-F_met*(COM_x(i,4)-coordinates_x(i,4));%重力方向成分
        M_connect(i,9)=m_list(10)*accel_CoM(i,1)*(COM_y(i,4)-coordinates_y(i,4));%慣性モーメントx成分
        M_connect(i,10)=m_list(10)*accel_CoM(i,2)*(COM_x(i,4)-coordinates_x(i,4));%慣性モーメントy成分
        M_connect(i,32)=M_connect(i,5)+M_connect(i,6)+M_connect(i,7)+M_connect(i,8)+data_T_all(i,93);

        %膝関節周りのモーメントの計算
        M_connect(i,2)=M_connect(i,1)+J_list(3)*ddphi(i,3)+m_list(7)*(accel_CoM(i,3)*(COM_y(i,3)-coordinates_y(i,3))+accel_CoM(i,4)*(COM_x(i,3)-coordinates_x(i,3)))...
                        +F_tib*(COM_x(i,3)-coordinates_x(i,3))-F_connect(i,1)*(coordinates_y(i,4)-coordinates_y(i,3))+F_connect(i,2)*(coordinates_x(i,4)-coordinates_x(i,3));
        M_connect(i,11)=F_connect(i,1)*(coordinates_y(i,4)-coordinates_y(i,3)); %関節力の水平方向成分
        M_connect(i,12)=-F_connect(i,2)*(coordinates_x(i,4)-coordinates_x(i,3)); %関節力の垂直方向成分
        M_connect(i,13)=-F_tib*(COM_x(i,3)-coordinates_x(i,3));%重力方向成分
        M_connect(i,14)=m_list(7)*accel_CoM(i,3)*(COM_y(i,3)-coordinates_y(i,3));%慣性モーメントx成分
        M_connect(i,15)=m_list(7)*accel_CoM(i,4)*(COM_x(i,3)-coordinates_x(i,3));%慣性モーメントy成分
        M_connect(i,33)=M_connect(i,11)+M_connect(i,12)+M_connect(i,13);


        %股関節周りのモーメントの計算
        M_connect(i,3)=M_connect(i,2)+J_list(2)*ddphi(i,2)+m_list(6)*(accel_CoM(i,7)*(COM_y(i,2)-coordinates_y(i,2))+accel_CoM(i,8)*(COM_x(i,2)-coordinates_x(i,2)))...
                        +F_fem*(COM_x(i,2)-coordinates_x(i,2))-F_connect(i,3)*(coordinates_y(i,3)-coordinates_y(i,2))+F_connect(i,4)*(coordinates_x(i,3)-coordinates_x(i,2))...
                            +F_dist(i,1)*(coordinates_y(i,10)-coordinates_y(i,2))-F_dist(i,2)*(coordinates_x(i,10)-coordinates_x(i,2))... %F_Ciによるトルク
                            +F_dist(i,3)*(coordinates_y(i,11)-coordinates_y(i,2))-F_dist(i,4)*(coordinates_x(i,11)-coordinates_x(i,2)); %F_GEによるトルク
        M_connect(i,16)=F_connect(i,3)*(coordinates_y(i,3)-coordinates_y(i,2)); %関節力の水平方向成分
        M_connect(i,17)=-F_connect(i,4)*(coordinates_x(i,3)-coordinates_x(i,2)); %関節力の垂直方向成分
        M_connect(i,18)=-F_fem*(COM_x(i,2)-coordinates_x(i,2));%重力方向成分
        M_connect(i,19)=m_list(6)*accel_CoM(i,7)*(COM_y(i,2)-coordinates_y(i,2));%慣性モーメントx成分
        M_connect(i,20)=m_list(6)*accel_CoM(i,8)*(COM_x(i,2)-coordinates_x(i,2));%慣性モーメントy成分
        M_connect(i,21)=-F_dist(i,1)*(coordinates_y(i,10)-coordinates_y(i,2)); %Fciのx方向成分
        M_connect(i,22)=F_dist(i,2)*(coordinates_x(i,10)-coordinates_x(i,2));%Fciのy方向成分
        M_connect(i,23)=-F_dist(i,3)*(coordinates_y(i,11)-coordinates_y(i,2));%F_GEoのx成分
        M_connect(i,24)=F_dist(i,4)*(coordinates_x(i,11)-coordinates_x(i,2));%F_GEoのy成分
        M_connect(i,34)=M_connect(i,16)+M_connect(i,17)+M_connect(i,18)+M_connect(i,21)+M_connect(i,22)+M_connect(i,23)+M_connect(i,24);


        %原点周りのモーメントの計算
        M_connect(i,4)=M_connect(i,3)+J_list(1)*ddphi(i,1)+(m_list(4)+m_list(5))*(accel_CoM(i,11)*(COM_y(i,1)-coordinates_y(i,1))+accel_CoM(i,12)*(COM_x(i,1)-coordinates_x(i,1)))...
                        +F_frame*(COM_x(i,1)-coordinates_x(i,1))-F_connect(i,5)*(coordinates_y(i,2)-coordinates_y(i,1))+F_connect(i,6)*(coordinates_x(i,2)-coordinates_x(i,1))...
                            -data_T_all(i,91); %theta1を一定に保つ拘束トルク
        M_connect(i,25)=F_connect(i,5)*(coordinates_y(i,2)-coordinates_y(i,1)); %関節力の水平方向成分
        M_connect(i,26)=-F_connect(i,6)*(coordinates_x(i,2)-coordinates_x(i,1)); %関節力の垂直方向成分
        M_connect(i,27)=-F_frame*(COM_x(i,1)-coordinates_x(i,1));%重力方向成分
        M_connect(i,28)=(m_list(4)+m_list(5))*accel_CoM(i,11)*(COM_y(i,1)-coordinates_y(i,1));%慣性モーメントx成分
        M_connect(i,29)=(m_list(4)+m_list(5))*accel_CoM(i,12)*(COM_x(i,1)-coordinates_x(i,1));%慣性モーメントy成分
        M_connect(i,35)=M_connect(i,25)+M_connect(i,26)+M_connect(i,27)+data_T_all(i,91);

        %M_connect(i,1)=M_connect(i,1)-data_T_all(i,93);
        %つま先周りのモーメントの計算
        %時計回りのモーメント
        M_connect(i,36)=-(F_dist(i,9)*(coordinates_y(i,8)-coordinates_y(i,5))+F_dist(i,10)*(coordinates_x(i,8)-coordinates_x(i,5))+data_T_all(i,93));%+data_F_all(i,5)*r;%+M_connect(i,1));
        %反時計回りのモーメント
        M_connect(i,37)=F_connect(i,1)*(coordinates_y(i,4)-coordinates_y(i,5))+F_connect(i,2)*(coordinates_x(i,4)-coordinates_x(i,5))-F_met*(COM_x(i,4)-coordinates_x(i,5));
        %時計回りのモーメント
        M_connect(i,38)=-(F_dist(i,9)*(coordinates_y(i,8)-coordinates_y(i,5))+F_dist(i,10)*(coordinates_x(i,8)-coordinates_x(i,5)))+GRF(i,2)*(coordinates_x(i,2)-coordinates_x(i,5));%+GRF(i,4)*(coordinates_x(i,4)-coordinates_x(i,5)));
        %反時計回りのモーメント
        M_connect(i,39)=-(F_frame+F_tib+F_fem-GRF(i,2))*(COM(1)-coordinates_x(i,5))+GRF(i,1)*(coordinates_y(i,2)-coordinates_y(i,5));


    end
    
end