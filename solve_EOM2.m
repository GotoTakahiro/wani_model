%外力のデータを格納
function solve_EOM = solve_EOM2(tmax,tspace,tspan,initial_condition,x_fixed,y_fixed,k_ground,c_ground,mu,l_link_list,l_muscle_list,limit_list,m_list,default_wire_k,default_wire_c,g,t_CFL,k_frame,c_frame,default_frame_angle,filename,end_CFL,CFL_alpha,t_end_exp,default_CFL,gain_list)

    t_message = 0.5;

    % 外乱の設定 unit: N
    max_disturbance = 4;
    dis_x_vec = [0,0,-1,0.5,0.8];
    dis_y_vec = [1,-1,0,0,0];
    x_disturbance = max_disturbance*dis_x_vec;
    y_disturbance = max_disturbance*dis_y_vec;
    % x_disturbance = max_disturbance*(2*rand(5,1)-1);
    % y_disturbance = max_disturbance*(2*rand(5,1)-1);

    disp(['Ci: ', num2str(l_muscle_list(1)), ' CFLT: ', num2str(l_muscle_list(2)), ' GEo: ', num2str(l_muscle_list(3)), ' GE: ', num2str(l_muscle_list(4))]);

    % 保存用の配列を定義
    data_error_CFL = [];
    data_l_CFL= 0;
    data_k_c_wire = [];
    data_accel_GRF = [];
    data_accel = [];
    data_force= [];
    data_Q = [];
    dqdt_rec = [];
    k_list = [];
    c_list = [];
    GRF = [];
    GRF_external = [];
    Q = zeros(10,1);
    Q_hip_pull = zeros(10,1);
    Q_hip_up = zeros(10,1);
    Q_heel = zeros(10,1);
    Q_toe = zeros(10,1);
    error_CFL = 0;
    target_CFL=0;
    pullForce_x=0;
    pullForce_y=0;
    T_gravity_all=zeros(60,1);
    F_all=zeros(5,1);
    T_all=zeros(93,1);


    data_accel_GRF(1,:) = [0, zeros(1,10), zeros(1,6)];
    data_accel(1,:) = [0, zeros(1,10)];
    data_force(1,:) = [0, zeros(1,10)];
    data_Q(1,:) = [0, zeros(1,10)];
    data_Q_hip_pull(1,:)=zeros(1,10);
    data_Q_hip_up(1,:)=zeros(1,10);
    data_Q_toe(1,:)=zeros(1,10);
    data_Q_heel(1,:)=zeros(1,10);
    data_k_c_wire(1,:) = [0, zeros(1,7), zeros(1,8)];
    data_error_CFL(1,:) = 0;
    data_l_CFL(1,:) = 0;
    data_T_gravity_all(1,:)=zeros(1,60);
    data_F_all(1,:)=zeros(1,5);
    data_T_all(1,:)=zeros(1,93);

    function [ode, output] = getHandles()
        ode = @func_EOM;
        output = @get_data;

        function [dqdt] = func_EOM(t,q)

            if isempty(t_message)
                t_message = 0.5;
            end

            %目標CFL長さの設定（CFLを徐々に縮めていく）
            if t < t_CFL
                target_CFL = default_CFL;
            elseif t >= t_CFL && t < end_CFL
                target_CFL = default_CFL - CFL_alpha*(t-t_CFL);
            elseif t >= end_CFL
                target_CFL = default_CFL - CFL_alpha*(end_CFL-t_CFL) - CFL_alpha*(t_end_exp-end_CFL)*(1-exp(-(t-end_CFL)));
            end

            % CFLの収縮速度誤差の計算．
            % derror_CFLはCFLの収縮速度誤差．CFLの目標長さの式（論文参照）を時間微分したものから得る．
            error_CFL = target_CFL - q(10);
            int_error_CFL = q(21);
            if t < t_CFL
                derror_CFL = 0;
            elseif t >= t_CFL && t < end_CFL
                derror_CFL = (- CFL_alpha) - q(20);
            elseif t >= end_CFL
                derror_CFL = (- CFL_alpha*(t_end_exp-end_CFL)*exp(-(t-end_CFL))) - q(20);
            end
                
            % 外力(股関節を引く力，地面反力，支持力)を考えるために必要な値を計算
            y_heel = calc_y_heel(l_link_list, q(1:10));
            y_toe = calc_y_toe(l_link_list, q(1:10));
            dx_heel_vec = calc_dx_heel_vec(l_link_list, q(1:10), q(11:20));
            dx_toe_vec = calc_dx_toe_vec(l_link_list, q(1:10), q(11:20));
            y_hip = q(2)-l_link_list(6)*cos(q(5));
            x_hip = q(1)+l_link_list(6)*sin(q(5));
            dx_hip_vec = calc_dx_hip_vec(l_link_list, q(1:10), q(11:20));

            % 股関節を後ろに引く力を設定
            if x_hip > x_fixed(1)
                pullForce_x = -10*(1-exp(-40*(abs(x_hip)))); %x方向に引っ張る力．
                pullForce_y = 0;
            else
                pullForce_x = -k_ground*(x_hip-x_fixed(1)) - c_ground*dx_hip_vec(1); %ストッパー．股関節が初期値位置よりも後ろに動かないようにする．
                pullForce_y = 0;
            end

            % 筋腱の起始停止点からたるみを考慮．たるんでいればばね定数と粘性係数が０になる．
            [k_list, c_list] = calc_spring_const(q(1:10),l_link_list,l_muscle_list,limit_list,default_wire_k,default_wire_c,k_frame,c_frame);
            % k_list(2) = 0;
            % c_list(2) = 0;

            % 運動方程式の係数行列．calc_EOM_workでシンボリックで計算した後関数化したものを呼び出す．
            M = Inertial_matrix(m_list,l_link_list,q(1:10));
            C = Coriolis_matrix(m_list,l_link_list,q(1:10),q(11:20));
            G = Stiffness_matrix(m_list, l_link_list, l_muscle_list, k_list, limit_list, g, q(1:10), q(11:20));
            D = Damping_matrix(l_link_list, c_list, q(1:10), q(11:20));

            % 地面反力
            GRF_external = func_GRF(q,k_ground,c_ground,mu,y_fixed,y_hip,y_heel,y_toe,dx_heel_vec,dx_toe_vec,dx_hip_vec);

            % 地面反力を後ろに引っ張る力も加えて再定義
            GRF = [GRF_external(1)+pullForce_x;GRF_external(2)+pullForce_y;GRF_external(3);GRF_external(4);GRF_external(5);GRF_external(6)]; 
            GRF1 = [pullForce_x;pullForce_y;0;0;0;0]; 
            GRF2 = [GRF_external(1);GRF_external(2);0;0;0;0]; 
            GRF3 = [0;0;GRF_external(3);GRF_external(4);0;0]; 
            GRF4 = [0;0;0;0;GRF_external(5);GRF_external(6);]; 
            %重力によるトルクの計算
            T_gravity_all = calc_torque_gravity_all(m_list,l_link_list, l_muscle_list, q(1:10),g); %m_list,l_link_list, l_muscle_list, q,g
            %ワイヤー張力の計算
            F_all = calc_force_all(l_link_list, l_muscle_list, k_list, c_list, q(1:10),q(11:20), gain_list,error_CFL,int_error_CFL,derror_CFL);
            %ワイヤ張力によるトルクと，拘束力によるトルクの計算
            T_all = calc_torque_force_all(l_link_list, l_muscle_list, k_list, c_list,limit_list, q(1:10),q(11:20), gain_list, error_CFL, int_error_CFL, derror_CFL);
            % 外力ベクトル．地面反力なども考慮し，各一般化座標に作用するトルクを計算．
            Q = External_force_matrix(GRF,l_link_list,q(1:10),gain_list,error_CFL,int_error_CFL,derror_CFL);
            %各成分の寄与分を格納
            %⑨，10，股関節を引く力によるトルク
            Q_hip_pull = External_force_matrix(GRF1,l_link_list,q(1:10),gain_list,error_CFL,int_error_CFL,derror_CFL);
            %11,股関節を支えるトルク
            Q_hip_up = External_force_matrix(GRF2,l_link_list,q(1:10),gain_list,error_CFL,int_error_CFL,derror_CFL);
            %12,かかとで発生する地面反力によるトルク
            Q_heel = External_force_matrix(GRF3,l_link_list,q(1:10),gain_list,error_CFL,int_error_CFL,derror_CFL);
            %13,つま先で発生する地面反力によるトルク
            Q_toe = External_force_matrix(GRF4,l_link_list,q(1:10),gain_list,error_CFL,int_error_CFL,derror_CFL);
            % 一般化座標の微分．1-10が速度，11-20が加速度．
            dqdt = [q(11);q(12);q(13);q(14);q(15);q(16);q(17);q(18);q(19);q(20);M\(-C+G+D+Q);error_CFL];
            dqdt_rec = dqdt(11:20,1).'; %一般化加速度を記録

        end

        % 0.01秒（tspace）ごとに呼び出され，データを記録．
        function status = get_data(t,q,flag)
            status = 0;
            % if flag ~= strcmp(flag,'init')|strcmp(flag,'done')
            if isempty(flag) || strcmp(flag,'')
                % disp(size(t));
                % disp(t)
                % disp(size(dqdt_rec));
                % disp(size(GRF));
                if size(t) == [1 2] %配列の形がたまにバグるのでその補正
                    t = t.';
                    data_accel_GRF(end+1,:) = [t(1), dqdt_rec, GRF.'];
                    data_accel(end+1,:) = [t(1), dqdt_rec];
                    data_force(end+1,:) = [t(1), pullForce_x,pullForce_y,GRF_external.'];
                    data_k_c_wire(end+1,:) = [t(1), k_list.', c_list.'];
                    data_Q(end+1,:) = [t(1), Q.'];
                    data_accel_GRF(end+1,:) = [t(2), dqdt_rec, GRF.'];
                    data_accel(end+1,:) = [t(2), dqdt_rec];
                    data_force(end+1,:) = [t(2), pullForce_x,pullForce_y,GRF_external.'];
                    data_k_c_wire(end+1,:) = [t(2), k_list.', c_list.'];
                    data_Q(end+1,:) = [t(2), Q.'];
                    data_error_CFL(end+1,:) = error_CFL;
                    data_l_CFL(end+1,:) = target_CFL;
                    data_Q_hip_pull(end+1,:)=Q_hip_pull.';
                    data_Q_hip_up(end+1,:)=Q_hip_up.';
                    data_Q_heel(end+1,:)=Q_heel.';
                    data_Q_toe(end+1,:)=Q_toe.';
                    data_T_gravity_all(end+1,:)=T_gravity_all.';
                    data_F_all(end+1,:)=F_all.';
                    data_T_all(end+1,:)=T_all.';



                else
                    data_accel_GRF(end+1,:) = [t, dqdt_rec, GRF.'];
                    data_accel(end+1,:) = [t, dqdt_rec];
                    data_force(end+1,:) = [t, pullForce_x,pullForce_y,GRF_external.'];
                    data_k_c_wire(end+1,:) = [t, k_list.', c_list.'];
                    data_Q(end+1,:) = [t, Q.'];
                    data_error_CFL(end+1,:) = error_CFL;
                    data_l_CFL(end+1,:) = target_CFL;
                    data_Q_hip_pull(end+1,:)=Q_hip_pull.';
                    data_Q_hip_up(end+1,:)=Q_hip_up.';
                    data_Q_heel(end+1,:)=Q_heel.';
                    data_Q_toe(end+1,:)=Q_toe.';
                    data_T_gravity_all(end+1,:)=T_gravity_all.';
                    data_F_all(end+1,:)=F_all.';
                    data_T_all(end+1,:)=T_all.';
                end
                status = 0;
                disp(t);
            end
        end
    
    end

    [ode, output] = getHandles();

    % 初期値を定義
    q = initial_condition;
    options = odeset('RelTol',1e-3,'AbsTol',1e-5,'OutputFcn',output);  % 計算に時間がかかるので誤差を大きく設定．
    % options = [];

    % ode45で運動方程式を解く
    [t,q] = ode45(ode,tspan,q,options);
    r = l_link_list(7);

    % データを保存
    save(filename,'t','q','default_wire_k','default_wire_c','l_link_list','l_muscle_list','m_list','r','g','k_ground','c_ground','mu','k_frame','c_frame','default_frame_angle','limit_list','x_fixed','y_fixed','t_CFL','data_accel_GRF','data_k_c_wire','filename','x_disturbance','y_disturbance','data_Q','gain_list','data_error_CFL','initial_condition','data_accel','data_force','data_l_CFL',...
        "data_Q_hip_pull","data_Q_hip_up","data_Q_heel","data_Q_toe","data_T_all","data_F_all","data_T_gravity_all");
    solve_EOM = true;
end