%F_GEの作用を検証するためのプログラム
clear all

%%%%%%%%%%%
% 3つの質点が論文ではM0, M1, M2になっているのが，こっちではM1, M2, M3になってます．
%%%%%%%%%%% 

syms x1 y1 x3 y3 theta1 theta2 theta3 theta4 theta_CFL l_CFL %CFL's origin, CFL's Branching point, Point of Y-shape junction, frame angle, hip angle, knee angle, ankle angle
syms dx1 dy1 dx3 dy3 dtheta1 dtheta2 dtheta3 dtheta4 dtheta_CFL dl_CFL %CFL's origin, CFL's Branching point, Point of Y-shape junction, frame angle, hip angle, knee angle, ankle angle
syms g
syms M1 M2 M3 M_hip M_fem M_tib M_met_pulley M_met_rod M_met M_frame %M1:CFL's origin, M2:CFL's Branching point, M3:Point of Y-shape junction, frame: caudal vertebral column（長くなるのでframe）
syms L_fem L_tib L_met L_frame %Length of femur, tibia, metatarsus, frame
syms L_4th_troch %Length from hip to 4th trochanter
syms L_GE_origin %Length from hip to origin of GE
syms L_Ci L_CFLT L_GEo L_GE L_CFL_target
% syms xh yh %Origin of CFL
syms r %Radius of pulley
syms k_Ci k_CFLT k_GEo k_GE k_ankle_flex k_ankle_ex k_frame %Spring constants
syms c_Ci c_CFLT c_GEo c_GE c_ankle_flex c_ankle_ex c_frame c_CFL %Damping constants
syms ankle_limit_flex ankle_limit_ex %Ankle joint limits
syms F_hip_x F_hip_y F_heel_x F_heel_y F_toe_x F_toe_y %Ground reaction forces
syms c_m1 c_m2 c_m3 c_theta1 c_theta2 c_theta3 c_theta4 %Damping constants
syms default_frame_angle
syms Pgain Igain Dgain error_CFL int_error_CFL derror_CFL

q = [x1; y1; x3; y3; theta1; theta2; theta3; theta4; theta_CFL; l_CFL];
dq = [dx1; dy1; dx3; dy3; dtheta1; dtheta2; dtheta3; dtheta4; dtheta_CFL; dl_CFL];
k_list = [k_Ci; k_CFLT; k_GEo; k_GE; k_ankle_ex; k_ankle_flex; k_frame];
c_list = [c_Ci; c_CFLT; c_GEo; c_GE; c_ankle_ex; c_ankle_flex; c_frame; c_CFL];
l_link_list = [L_fem; L_tib; L_met; L_4th_troch; L_GE_origin; L_frame; r];
m_list = [M1; M2; M3; M_hip; M_frame; M_fem; M_tib; M_met_pulley; M_met_rod; M_met];
l_muscle_list = [L_Ci; L_CFLT; L_GEo; L_GE];
limit_list = [ankle_limit_ex; ankle_limit_flex; default_frame_angle];
gain_list = [Pgain; Igain; Dgain];
% gain_list = [Pgain; Igain];

%慣性モーメント，剛体棒，円を仮定
J_fem = M_fem*L_fem^2/12;
J_tib = M_tib*L_tib^2/12;
J_met_pulley = M_met_pulley*r^2/2;
J_met_rod = M_met_rod*L_met^2/12;

%metatarsalの端部から重心までの長さ
L_met_CoM = (L_met*M_met_rod/2)/(M_met_pulley + M_met_rod);
%metatarsal全体の慣性モーメント
J_met = (J_met_rod + M_met_rod*(L_met/2 - L_met_CoM)^2) + (J_met_pulley + M_met_pulley*L_met_CoM^2);

%frameの端部から重心までの長さ
L_frame_CoM = (L_frame*M_frame/2 + M_hip*L_frame)/(M_frame + M_hip);
J_frame_only = M_frame*L_frame^2/12;
%frame全体の慣性モーメント
J_frame = (J_frame_only + M_frame*(L_frame_CoM - L_frame/2)^2) + (M_hip*(L_frame - L_frame_CoM)^2);

J_list = [J_frame; J_fem; J_tib; J_met;];

%m2の座標，速度をtheta_CFL, l_CFLで計算
x2 = x1 + l_CFL*sin(theta_CFL);
y2 = y1 - l_CFL*cos(theta_CFL);
dx2 = jacobian(x2, q)*dq;
dy2 = jacobian(y2, q)*dq;

%フレーム重心の座標，速度
x_frame_CoM = x1 + L_frame_CoM*sin(theta1);
y_frame_CoM = y1 - L_frame_CoM*cos(theta1);
dx_frame_CoM = jacobian(x_frame_CoM, q)*dq;
dy_frame_CoM = jacobian(y_frame_CoM, q)*dq;
dx_frame_CoM_vec = [dx_frame_CoM; dy_frame_CoM];

%股関節重心の座標，速度
x_hip = x1 + L_frame*sin(theta1);
y_hip = y1 - L_frame*cos(theta1);
dx_hip = jacobian(x_hip, q)*dq;
dy_hip = jacobian(y_hip, q)*dq;
dx_hip_vec = [dx_hip; dy_hip];

%大腿骨重心の座標，速度
x_fem_CoM = x_hip + L_fem/2*sin(theta1+theta2);
y_fem_CoM = y_hip - L_fem/2*cos(theta1+theta2);
dx_fem_CoM = jacobian(x_fem_CoM, q)*dq;
dy_fem_CoM = jacobian(y_fem_CoM, q)*dq;
dx_fem_CoM_vec = [dx_fem_CoM; dy_fem_CoM];
%膝重心の座標，速度
x_knee = x_hip + L_fem*sin(theta1+theta2);
y_knee = y_hip - L_fem*cos(theta1+theta2);
%脛骨重心の座標，速度
x_tib_CoM = x_knee + L_tib/2*sin(theta1+theta2+theta3);
y_tib_CoM = y_knee - L_tib/2*cos(theta1+theta2+theta3);
dx_tib_CoM = jacobian(x_tib_CoM, q)*dq;
dy_tib_CoM = jacobian(y_tib_CoM, q)*dq;
dx_tib_CoM_vec = [dx_tib_CoM; dy_tib_CoM];

%くるぶし重心の座標，速度
x_ankle = x_knee + L_tib*sin(theta1+theta2+theta3);
y_ankle = y_knee - L_tib*cos(theta1+theta2+theta3);
%踵重心の座標，速度
x_heel = x_ankle;
y_heel = y_ankle - r;
dx_heel = jacobian(x_heel, q)*dq;
dy_heel = jacobian(y_heel, q)*dq;
dx_heel_vec = [dx_heel; dy_heel];
%足先重心の座標，速度
x_met_CoM = x_ankle + L_met_CoM*sin(theta1+theta2+theta3+theta4);
y_met_CoM = y_ankle - L_met_CoM*cos(theta1+theta2+theta3+theta4);
dx_met_CoM = jacobian(x_met_CoM, q)*dq;
dy_met_CoM = jacobian(y_met_CoM, q)*dq;
dx_met_CoM_vec = [dx_met_CoM; dy_met_CoM];

%つま先重心の座標，速度
x_toe = x_ankle + L_met*sin(theta1+theta2+theta3+theta4);
y_toe = y_ankle - L_met*cos(theta1+theta2+theta3+theta4);
dx_toe = jacobian(x_toe, q)*dq;
dy_toe = jacobian(y_toe, q)*dq;
dx_toe_vec = [dx_toe; dy_toe];


%Coordinate of 4th trochanter
x_4th_troch = x_hip + L_4th_troch*sin(theta1+theta2);
y_4th_troch = y_hip - L_4th_troch*cos(theta1+theta2);
Coordinate_4th_troch = [x_4th_troch; y_4th_troch];
dx_4th_troch = jacobian(x_4th_troch, q)*dq;
dy_4th_troch = jacobian(y_4th_troch, q)*dq;
dx_4th_troch_vec = [dx_4th_troch; dy_4th_troch];

%Coordinate of GE origin
x_GE_origin = x_hip + L_GE_origin*sin(theta1+theta2);
y_GE_origin = y_hip - L_GE_origin*cos(theta1+theta2);
Coordinate_GE_origin = [x_GE_origin; y_GE_origin];
dx_GE_origin = jacobian(x_GE_origin, q)*dq;
dy_GE_origin = jacobian(y_GE_origin, q)*dq;
dx_GE_origin_vec = [dx_GE_origin; dy_GE_origin];

%Coordinates of tangent line drawn from M3 to the pulley (calculated in circle_point_of_tangency.m)
% M3（論文ではM2）を通り，足関節に配置したプーリに引かれる接線のうち，プーリの中心からM3に向かうベクトルを基準として紙面向かって反時計まわりに位置する方の接点の座標，速度
x_M3_to_pulley = (r^2*x3 - r^2*x_ankle + x_ankle*x3^2 - 2*x_ankle^2*x3 + x_ankle*y_ankle^2 + x_ankle*y3^2 + x_ankle^3 + r*y_ankle*(- r^2 + x_ankle^2 - 2*x_ankle*x3 + x3^2 + y_ankle^2 - 2*y_ankle*y3 + y3^2)^(1/2) - r*y3*(- r^2 + x_ankle^2 - 2*x_ankle*x3 + x3^2 + y_ankle^2 - 2*y_ankle*y3 + y3^2)^(1/2) - 2*x_ankle*y_ankle*y3)/(x_ankle^2 - 2*x_ankle*x3 + x3^2 + y_ankle^2 - 2*y_ankle*y3 + y3^2);
y_M3_to_pulley = (r^2*y3 - r^2*y_ankle + x_ankle^2*y_ankle + x3^2*y_ankle + y_ankle*y3^2 - 2*y_ankle^2*y3 + y_ankle^3 - r*x_ankle*(- r^2 + x_ankle^2 - 2*x_ankle*x3 + x3^2 + y_ankle^2 - 2*y_ankle*y3 + y3^2)^(1/2) + r*x3*(- r^2 + x_ankle^2 - 2*x_ankle*x3 + x3^2 + y_ankle^2 - 2*y_ankle*y3 + y3^2)^(1/2) - 2*x_ankle*x3*y_ankle)/(x_ankle^2 - 2*x_ankle*x3 + x3^2 + y_ankle^2 - 2*y_ankle*y3 + y3^2);
coodinate_M3_to_pulley = [x_M3_to_pulley; y_M3_to_pulley];
dx_M3_to_pulley = jacobian(x_M3_to_pulley, q)*dq;
dy_M3_to_pulley = jacobian(y_M3_to_pulley, q)*dq;
dx_M3_to_pulley_vec = [dx_M3_to_pulley; dy_M3_to_pulley];

%Coordinates of tangent line drawn from toe to the pulley (calculated in circle_point_of_tangency.m)
% つま先からプーリに引かれる接線のうち，プーリの中心からつま先に向かうベクトルを基準として紙面向かって時計まわりに位置する方の接点の座標
x_toe_to_pulley = (r^2*x_toe - r^2*x_ankle + x_ankle*x_toe^2 - 2*x_ankle^2*x_toe + x_ankle*y_ankle^2 + x_ankle*y_toe^2 + x_ankle^3 - r*y_ankle*(- r^2 + x_ankle^2 - 2*x_ankle*x_toe + x_toe^2 + y_ankle^2 - 2*y_ankle*y_toe + y_toe^2)^(1/2) + r*y_toe*(- r^2 + x_ankle^2 - 2*x_ankle*x_toe + x_toe^2 + y_ankle^2 - 2*y_ankle*y_toe + y_toe^2)^(1/2) - 2*x_ankle*y_ankle*y_toe)/(x_ankle^2 - 2*x_ankle*x_toe + x_toe^2 + y_ankle^2 - 2*y_ankle*y_toe + y_toe^2);
y_toe_to_pulley = (r^2*y_toe - r^2*y_ankle + x_ankle^2*y_ankle + x_toe^2*y_ankle + y_ankle*y_toe^2 - 2*y_ankle^2*y_toe + y_ankle^3 + r*x_ankle*(- r^2 + x_ankle^2 - 2*x_ankle*x_toe + x_toe^2 + y_ankle^2 - 2*y_ankle*y_toe + y_toe^2)^(1/2) - r*x_toe*(- r^2 + x_ankle^2 - 2*x_ankle*x_toe + x_toe^2 + y_ankle^2 - 2*y_ankle*y_toe + y_toe^2)^(1/2) - 2*x_ankle*x_toe*y_ankle)/(x_ankle^2 - 2*x_ankle*x_toe + x_toe^2 + y_ankle^2 - 2*y_ankle*y_toe + y_toe^2);
coodinate_toe_to_pulley = [x_toe_to_pulley; y_toe_to_pulley];

%calculate beta angle (Angle of wire through pulley)
%プーリーにかかっているワイヤの角度を求める
vec_pulley_to_M3 = [x3-x_M3_to_pulley; y3-y_M3_to_pulley];
vec_pulley_to_toe = [x_toe - x_toe_to_pulley; y_toe - y_toe_to_pulley];
cross_wire_pulley = cross([vec_pulley_to_M3; 0], [vec_pulley_to_toe; 0]);
beta = acos(dot(vec_pulley_to_M3, vec_pulley_to_toe)/(norm(vec_pulley_to_M3)*norm(vec_pulley_to_toe)));
%巻きつき角
angle_wire_pulley = pi - beta;

%単位ベクトル
uvec_m1_to_m2 = [(x2-x1); (y2-y1)]/sqrt((x2-x1)^2 + (y2-y1)^2);
uvec_4thtro_to_m2 = [(x2-x_4th_troch); (y2-y_4th_troch)]/sqrt((x2-x_4th_troch)^2 + (y2-y_4th_troch)^2);
uvec_m2_to_m3 = [(x3-x2); (y3-y2)]/sqrt((x3-x2)^2 + (y3-y2)^2);
uvec_GEo_to_m3 = [(x3-x_GE_origin); (y3-y_GE_origin)]/sqrt((x3-x_GE_origin)^2 + (y3-y_GE_origin)^2);
uvec_toe_to_pulley = -vec_pulley_to_toe/norm(vec_pulley_to_toe);
uvec_pulley_to_toe = vec_pulley_to_toe/norm(vec_pulley_to_toe);
uvec_m3_to_pulley = -vec_pulley_to_M3/norm(vec_pulley_to_M3);

%Kinetic energy，運動エネルギーの計算
T_M1 = 1/2*M1*(dx1^2 + dy1^2);
T_M2 = 1/2*M2*(dx2^2 + dy2^2);
T_M3 = 1/2*M3*(dx3^2 + dy3^2);
T_frame = 1/2*J_frame*dtheta1^2 + 1/2*(M_frame+M_hip)*(dx_frame_CoM^2 + dy_frame_CoM^2);
T_fem = 1/2*J_fem*(dtheta1+dtheta2)^2 + 1/2*M_fem*(dx_fem_CoM^2 + dy_fem_CoM^2);
T_tib = 1/2*J_tib*(dtheta1+dtheta2+dtheta3)^2 + 1/2*M_tib*(dx_tib_CoM^2 + dy_tib_CoM^2);
T_met = 1/2*J_met*(dtheta1+dtheta2+dtheta3+dtheta4)^2 + 1/2*M_met*(dx_met_CoM^2 + dy_met_CoM^2);
T = T_M1 + T_M2 + T_M3 + T_frame + T_fem + T_tib + T_met;

%Potential by gravity，位置エネルギーの計算
U_M1 = M1*g*y1;
U_M2 = M2*g*y2;
U_M3 = M3*g*y3;
U_frame = (M_frame+M_hip)*g*y_frame_CoM;
U_fem = M_fem*g*y_fem_CoM;
U_tib = M_tib*g*y_tib_CoM;
U_met = M_met*g*y_met_CoM;
U = U_M1 + U_M2 + U_M3 + U_frame + U_fem + U_tib + U_met;

%Spring displacement,各ワイヤが伸ばされた長さ
L_Ci_dis = sqrt((x_4th_troch - x2)^2 + (y_4th_troch - y2)^2) - L_Ci;
L_CFLT_dis = sqrt((x3 - x2)^2 + (y3- y2)^2) - L_CFLT;
L_GEo_dis = sqrt((x3 - x_GE_origin)^2 + (y3 - y_GE_origin)^2) - L_GEo;
% M3（論文ではM2）からプーリに引かれるワイヤの長さと，ワイヤとプーリの接点から足底（GE停止点）までの円弧で長さを評価．
L_GE_dis = (sqrt((x3 - x_M3_to_pulley)^2 + (y3 - y_M3_to_pulley)^2) + r*angle_wire_pulley) - L_GE;
L_GE_present = (sqrt((x3 - x_M3_to_pulley)^2 + (y3 - y_M3_to_pulley)^2) + r*angle_wire_pulley);
L_ankle_ex = theta4 - ankle_limit_ex;
L_ankle_flex = ankle_limit_flex - theta4;
L_frame_dis = theta1 - default_frame_angle;
%弾性力の計算
F_Ci_spring = k_Ci*L_Ci_dis;
F_CFLT_spring = k_CFLT*L_CFLT_dis;
F_GEo_spring = k_GEo*L_GEo_dis;
F_GE_spring = k_GE*L_GE_dis;

%CFLを目標長さに調整するようにPI制御でF_CFLを計算
%CFLの粘弾性についてもここで考慮
F_CFL_PID = Pgain*error_CFL + Igain*int_error_CFL + Dgain*derror_CFL;
% F_CFL_PI = Pgain*error_CFL + Igain*int_error_CFL;
%弾性力に単位ベクトルをかけ，(x,y)方向成分に分離
F_Ci_spring_vec = F_Ci_spring*uvec_4thtro_to_m2;
F_CFLT_spring_vec = F_CFLT_spring*uvec_m2_to_m3;
F_GEo_spring_vec = F_GEo_spring*uvec_GEo_to_m3;
F_GE_spring_vec = F_GE_spring*uvec_toe_to_pulley;
%GEがM3で作用するトルクの計算
F_GE_spring_vec_M3 = F_GE_spring*uvec_m3_to_pulley;
%T_GE=jacobian([x3; y3], q).'*F_GE_spring_vec_M3;
%厳密なGEによる発生トルクの計算
T_GE=-F_GE_spring*jacobian(L_GE_present, q).';

%Spring energy，弾性エネルギーの計算
U_Ci = 1/2*k_Ci*L_Ci_dis^2;
U_CFLT = 1/2*k_CFLT*L_CFLT_dis^2;
U_GEo = 1/2*k_GEo*L_GEo_dis^2;
U_GE =0;% 1/2*k_GE*L_GE_dis^2;
U_ankle_flex = 1/2*k_ankle_flex*L_ankle_flex^2;
U_ankle_ex = 1/2*k_ankle_ex*L_ankle_ex^2;
U_frame_spring = 1/2*k_frame*L_frame_dis^2;
U_spring = U_Ci + U_CFLT + U_GEo + U_GE + U_ankle_flex + U_ankle_ex + U_frame_spring;

%減衰関数を計算するために，ダンパに沿う方向の速度の(x,y)成分を求める
v_m2_to_4thtro = dot([dx2; dy2], -uvec_4thtro_to_m2);
v_4thtro_to_m2 = dot([dx_4th_troch; dy_4th_troch], uvec_4thtro_to_m2);
v_m2_to_m3 = dot([dx2; dy2], uvec_m2_to_m3);
v_m3_to_m2 = dot([dx3; dy3], -uvec_m2_to_m3);
v_m3_to_GEo = dot([dx3; dy3], -uvec_GEo_to_m3);
v_GEo_to_m3 = dot([dx_GE_origin; dy_GE_origin], uvec_GEo_to_m3);
v_m3_to_pulley = dot([dx3; dy3], uvec_m3_to_pulley);
v_pulley_to_m3 = dot([dx_M3_to_pulley; dy_M3_to_pulley], -uvec_m3_to_pulley);

% 質点や起始停止点に作用する減衰力を計算+(x,y)成分に分離
F_m2_damper_vec = - c_Ci*(v_m2_to_4thtro+v_4thtro_to_m2)*(-uvec_4thtro_to_m2) - c_CFLT*(v_m2_to_m3+v_m3_to_m2)*uvec_m2_to_m3;
F_m3_damper_vec = - c_CFLT*(v_m3_to_m2+v_m2_to_m3)*(-uvec_m2_to_m3) - c_GEo*(v_m3_to_GEo+v_GEo_to_m3)*(-uvec_GEo_to_m3) - c_GE*(v_m3_to_pulley+v_pulley_to_m3)*uvec_m3_to_pulley;
F_4thtro_damper_vec = - c_Ci*(v_4thtro_to_m2+v_m2_to_4thtro)*uvec_4thtro_to_m2;
F_GEo_damper_vec = - c_GEo*(v_GEo_to_m3+v_m3_to_GEo)*uvec_GEo_to_m3;
F_GE_pulley_damper_vec = - c_GE*(v_pulley_to_m3+v_m3_to_pulley)*(-uvec_m3_to_pulley);
%各ワイヤが発揮する減衰力
F_Ci_damper = c_Ci*(v_m2_to_4thtro+v_4thtro_to_m2);
F_CFLT_damper = c_CFLT*(v_m2_to_m3+v_m3_to_m2);
F_GEo_damper = c_GEo*(v_m3_to_GEo+v_GEo_to_m3);
F_GE_pulley_damper = c_GE*(v_m3_to_pulley+v_pulley_to_m3);

%減衰力が作用する点のヤコビアンを計算
jac_m2 = jacobian([x2; y2], q);
jac_m3 = jacobian([x3; y3], q);
jac_4th_troch = jacobian([x_4th_troch; y_4th_troch], q);
jac_GE_origin = jacobian([x_GE_origin; y_GE_origin], q);
jac_m3_to_pulley = jacobian([x_M3_to_pulley; y_M3_to_pulley], q);
jac_GE_insertion = jacobian([x_toe_to_pulley; y_toe_to_pulley], q);

% 仮想仕事の原理で各一般化座標に対する減衰力（一般化力）を計算，次元はトルク
D_m2 = jac_m2.'*F_m2_damper_vec;
D_m3 = jac_m3.'*F_m3_damper_vec;
D_4thtro = jac_4th_troch.'*F_4thtro_damper_vec;
D_GEo = jac_GE_origin.'*F_GEo_damper_vec;
D_GE_pulley = jac_m3_to_pulley.'*F_GE_pulley_damper_vec;

%粘弾性エネルギーの計算
dL_ankle_flex = dtheta4;
dL_ankle_ex = dtheta4;
dL_frame_dis = dtheta1;
dL_l_CFL = dl_CFL;
D_ankle_flex = 1/2*c_ankle_flex*dL_ankle_flex^2;
D_ankle_ex = 1/2*c_ankle_ex*dL_ankle_ex^2;
D_frame = 1/2*c_frame*dL_frame_dis^2;
%l_CFLの粘弾性による粘弾性エネルギーの計算
D_l_CFL = 1/2*c_CFL*dL_l_CFL^2; %c_CFL=0のため，0となる
%エネルギー散逸
D_lim = D_ankle_flex+D_ankle_ex + D_frame;
%一般化力への変換，単位はトルク
Jac_D_lim = jacobian(-D_lim, dq).';

%Lagrangian
%リンク系の運動E-リンク系の位置E-筋腱が元の長さを保とうとする張力による弾性E
L = T - U - U_spring;

%Coefficient matrix of the EoM
M = jacobian(jacobian(L, dq), dq);
C = jacobian(jacobian(L, dq), q)*dq;
G = jacobian(L, q).';
%話の単純化のためプーリーの粘性力を0に
D = D_m2 + D_m3 + D_4thtro + D_GEo  + Jac_D_lim; %+ D_GE_pulley

%地面反力の計算
Jac_hip = jacobian([x_hip; y_hip], q);
Jac_toe = jacobian([x_toe; y_toe], q);
Jac_heel = jacobian([x_heel; y_heel], q);

%地面反力
F_hip = [F_hip_x; F_hip_y];
F_toe = [F_toe_x; F_toe_y];
F_heel = [F_heel_x; F_heel_y];
F_list = [F_hip; F_heel; F_toe];

% 仮想仕事の原理で各一般化座標に対する地面反力による一般化力を計算，単位はトルク
Q_hip = Jac_hip.'*F_hip;
Q_toe = Jac_toe.'*F_toe;
Q_heel = Jac_heel.'*F_heel;
Q_CFL = [0; 0; 0; 0; 0; 0; 0; 0; 0; F_CFL_PID];
%GEがtheta4に対して発生させるトルクを計算
Q_GE = [0; 0; 0; 0; 0; 0; 0; -k_GE*L_GE_dis*r; 0; 0];
Q = Q_hip + Q_toe + Q_heel + Q_CFL+T_GE;

%張力のリストへの格納
F_Ci = F_Ci_spring;
F_CFLT = F_CFLT_spring;
F_GE_origin = F_GEo_spring;
F_GE_insertion = F_GE_spring;
muscle_tension = [F_Ci; F_CFLT; F_GE_origin; F_GE_insertion];

F_Ci_vec = F_Ci_spring_vec;
F_CFLT_vec = F_CFLT_spring_vec;
F_GE_origin_vec = F_GEo_spring_vec;
F_GE_insertion_vec = F_GE_spring_vec;

%筋腱の弾性力による一般化力を計算，単位はトルク
torque_4th_troch = jac_4th_troch.'*F_Ci_vec;
torque_GE_origin = jac_GE_origin.'*F_GE_origin_vec;
torque_GE_insertion = jac_GE_insertion.'*F_GE_insertion_vec;
torque_muscle = torque_4th_troch + torque_GE_origin + torque_GE_insertion;

%全体の重心を計算
COM_x = (M1*x1 + M2*x2 + M3*x3 + M_frame*x_frame_CoM + M_fem*x_fem_CoM + M_tib*x_tib_CoM + M_met*x_met_CoM)/(M1 + M2 + M3 + M_frame + M_fem + M_tib + M_met);
COM_y = (M1*y1 + M2*y2 + M3*y3 + M_frame*y_frame_CoM + M_fem*y_fem_CoM + M_tib*y_tib_CoM + M_met*y_met_CoM)/(M1 + M2 + M3 + M_frame + M_fem + M_tib + M_met);
COM = [COM_x; COM_y];


%②F_Ciの発揮力（伸展方向を正ととる，M2が基準点）
%相対速度を求める
rv_m2_to_4thtro = [dx_4th_troch-dx2;dy_4th_troch-dy2];
%M2 to 4th方向の速度成分
rv_Ci = dot(rv_m2_to_4thtro, -uvec_4thtro_to_m2);
F_Ci_all = -k_Ci*L_Ci_dis-c_Ci*rv_Ci;
%x,y成分に分解
F_Ci_all_vec = F_Ci_all*(-uvec_4thtro_to_m2);
%F_Ciが4th_trochで発生させるトルク
T_Ci_4th=jacobian([x_4th_troch; y_4th_troch], q).'*(-F_Ci_all_vec);
%F_CiがM2で発生させるトルク
T_Ci_M2=jacobian([x2; y2], q).'*F_Ci_all_vec;


%③F_CFLTの発揮力（伸展方向を正ととる，M2が基準点）
%相対速度を求める
rv_m2_to_m3 = [dx3-dx2;dy3-dy2];
%M2 to M3方向の速度成分
rv_CFLT = dot(rv_m2_to_m3, uvec_m2_to_m3);
F_CFLT_all = -k_CFLT*L_CFLT_dis-c_CFLT*rv_CFLT;
%x,y成分に分解
F_CFLT_all_vec = F_CFLT_all*uvec_m2_to_m3;
%F_CFLTがM2で発生させるトルク
T_CFLT_M2=jacobian([x2; y2], q).'*F_CFLT_all_vec;
%F_CFLTがM3で発生させるトルク
T_CFLT_M3=jacobian([x3; y3], q).'*(-F_CFLT_all_vec);

%④F_GEoの発揮力（伸展方向を正ととる，M3が基準点）
%相対速度を求める
rv_m3_to_GEo = [dx_GE_origin-dx3;dy_GE_origin-dy3];
%M3 to GEo方向の速度成分
rv_GEo = dot(rv_m3_to_GEo, -uvec_GEo_to_m3);
F_GEo_all = -k_GEo*L_GEo_dis-c_GEo*rv_GEo;
%x,y成分に分解
F_GEo_all_vec = F_GEo_all*(-uvec_GEo_to_m3);
%F_GEoがM3で発生させるトルク
T_GEo_M2=jacobian([x3; y3], q).'*F_GEo_all_vec;
%F_GEoがGEoriginで発生させるトルク
T_GEo_GEorigin=jacobian([x_GE_origin; y_GE_origin], q).'*(-F_GEo_all_vec);

%⑤F_GEの発揮力（伸展方向を正ととる，M3が基準点）
%相対速度を求める
rv_m3_to_GE = [dx_M3_to_pulley-dx3;dy_M3_to_pulley-dy3];
%M3 to GEo方向の速度成分
rv_GE = dot(rv_m3_to_GE, uvec_m3_to_pulley);
%一端粘性項は0にする
F_GE_all = -k_GE*L_GE_dis; %-c_GE*rv_GE;
%x,y成分に分解
F_GE_all_vec = F_GE_all*uvec_m3_to_pulley;
F_GE_all_vec2 = F_GE_all*uvec_pulley_to_toe;
%F_GEがM3で発生させるトルク
T_GE_M3=jacobian([x3; y3], q).'*F_GE_all_vec;
%F_GEがtoe_pulleyで発生させるトルク
%力の向きもpulley_to_toeを正方向にとったとき
T_GE_toe_pulley=jacobian([x_toe_to_pulley; y_toe_to_pulley], q).'*(-F_GE_all_vec2);
%F_GEがM3_pulleyで発生させるトルク
T_GE_M3_pulley=jacobian([x_M3_to_pulley; y_M3_to_pulley], q).'*(-F_GE_all_vec);
%F_GEがプーリーに対して生じさせるトルク
T_GE_pulley=F_GE_all*r;

%⑥frame(theta1)の角度を固定しようとする拘束力（M1に直接作用するトルク）
T_frame = -k_frame*L_frame_dis-c_frame*dL_frame_dis;
%⑦theta4が底屈しないようにするためのトルク
T_ankle_flex = -k_ankle_flex*L_ankle_flex-c_ankle_flex*dL_ankle_flex;
%⑧theta4が背屈しないようにするためのトルク
T_ankle_ex = -k_ankle_ex*L_ankle_ex-c_ankle_ex*dL_ankle_ex;
%ワイヤ張力とトルクの格納
F_all=[F_CFL_PID;F_Ci_all;F_CFLT_all;F_GEo_all;F_GE_all];
T_all=[T_Ci_4th;T_Ci_M2;T_CFLT_M2;T_CFLT_M3;T_GEo_M2;T_GEo_GEorigin;T_GE_M3;T_GE_toe_pulley;T_GE_M3_pulley;T_frame;T_ankle_flex;T_ankle_ex;T_GE_pulley];

%重力項のトルク計算
%フレーム＋股関節
Mg_frame=[0;-(M_frame+M_hip)*g];
Torque_frame=jacobian([x_frame_CoM;y_frame_CoM], q).'*Mg_frame;
%大腿骨
Mg_fem=[0;-M_fem*g];
Torque_fem=jacobian([x_fem_CoM;y_fem_CoM], q).'*Mg_fem;
%膝関節
Mg_tib=[0;-M_tib*g];
Torque_tib=jacobian([x_tib_CoM;y_tib_CoM], q).'*Mg_tib;
%足先
Mg_met=[0;-M_met*g];
Torque_met=jacobian([x_met_CoM;y_met_CoM], q).'*Mg_met;
%質点M2
Mg_M2=[0;-M2*g];
Torque_M2=jacobian([x2;y2], q).'*Mg_M2;
%質点M3
Mg_M3=[0;-M3*g];
Torque_M3=jacobian([x3;y3], q).'*Mg_M3;
%トルクの格納
T_gravity_all=[Torque_M2;Torque_M3;Torque_frame;Torque_fem;Torque_tib;Torque_met];
matlabFunction(T_gravity_all,'File','calc_torque_gravity_all', 'Vars', {m_list,l_link_list, l_muscle_list, q,g});


disp('writing...')
matlabFunction(M,'File','Inertial_matrix', 'Vars', {m_list, l_link_list, q});
matlabFunction(C,'File','Coriolis_matrix', 'Vars', {m_list, l_link_list, q, dq});
matlabFunction(G,'File','Stiffness_matrix', 'Vars', {m_list, l_link_list, l_muscle_list, k_list, limit_list, g, q, dq});
matlabFunction(D,'File','Damping_matrix', 'Vars', {l_link_list, c_list, q, dq});
matlabFunction(Q,'File','External_force_matrix', 'Vars', {F_list, l_link_list, q, k_list,l_muscle_list, gain_list, error_CFL, int_error_CFL, derror_CFL});
matlabFunction(dx_toe_vec,'File','calc_dx_toe_vec', 'Vars', {l_link_list, q, dq});
matlabFunction(dx_heel_vec,'File','calc_dx_heel_vec', 'Vars', {l_link_list, q, dq});
matlabFunction(dx_hip_vec,'File','calc_dx_hip_vec', 'Vars', {l_link_list, q, dq});
matlabFunction(Coordinate_4th_troch,'File','calc_coordinate_4th_troch', 'Vars', {l_link_list, q});
matlabFunction(Coordinate_GE_origin,'File','calc_coordinate_GE_origin', 'Vars', {l_link_list, q}); 
matlabFunction(coodinate_M3_to_pulley,'File','calc_coordinate_M3_to_pulley', 'Vars', {l_link_list, q});
matlabFunction(coodinate_toe_to_pulley,'File','calc_coordinate_toe_to_pulley', 'Vars', {l_link_list, q});
matlabFunction(angle_wire_pulley,'File','calc_angle_wire_pulley', 'Vars', {l_link_list, q});
matlabFunction(cross_wire_pulley,'File','calc_cross_wire_pulley', 'Vars', {l_link_list, q});
matlabFunction(y_heel,'File','calc_y_heel', 'Vars', {l_link_list, q});
matlabFunction(y_toe,'File','calc_y_toe', 'Vars', {l_link_list, q});
matlabFunction(torque_muscle,'File','calc_torque_muscle', 'Vars', {l_link_list, l_muscle_list, k_list, c_list, q, dq});
matlabFunction(muscle_tension,'File','calc_muscle_tension', 'Vars', {l_link_list, l_muscle_list, k_list, c_list, q, dq});
matlabFunction(COM,'File','calc_COM', 'Vars', {m_list, l_link_list, q});

matlabFunction(F_all,'File','calc_force_all', 'Vars', {l_link_list, l_muscle_list, k_list, c_list, q,dq, gain_list, error_CFL, int_error_CFL, derror_CFL});
matlabFunction(T_all,'File','calc_torque_force_all', 'Vars', { l_link_list, l_muscle_list, k_list, c_list,limit_list, q,dq, gain_list, error_CFL, int_error_CFL, derror_CFL});

matlabFunction(dx_frame_CoM_vec, 'File', 'calc_dx_frame_CoM_vec', 'Vars', {m_list,l_link_list, q, dq});
matlabFunction(dx_fem_CoM_vec, 'File', 'calc_dx_fem_CoM_vec', 'Vars', {l_link_list, q, dq});
matlabFunction(dx_tib_CoM_vec, 'File', 'calc_dx_tib_CoM_vec', 'Vars', {l_link_list, q, dq});
matlabFunction(dx_met_CoM_vec, 'File', 'calc_dx_met_CoM_vec', 'Vars', {m_list,l_link_list, q, dq});
matlabFunction(J_list, 'File', 'calc_J_list', 'Vars', {m_list, l_link_list});
%matlabFunction(momentum_list, 'File', 'calc_momentum_list', 'Vars', {m_list, l_link_list, q, dq});
%matlabFunction(angular_moment_list, 'File', 'calc_angular_moment_list', 'Vars', {m_list, l_link_list, q, dq});
%matlabFunction(momentum_COM_list, 'File', 'calc_momentum_COM_list', 'Vars', {m_list, l_link_list, q, dq});
%ht=matlabFunction(momentum_COM_list, 'File', 'calc_momentum_COM_list', 'Vars', {m_list, l_link_list, q, dq});