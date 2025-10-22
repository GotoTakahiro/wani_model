# はじめに
運動方程式の導出方法や考え方は木村のNOLTAの論文か修論を見てください．大体そこに書いてます．
シンボリックで運動方程式の導出などを行って，それを使って数値計算してます．

使用しているシンボリックについて，3つの質点が論文ではM0, M1, M2となっているが，ここではM1, M2, M3としている．変えてもらっても大丈夫だが，使用するファイル全てで整合を取る必要がある．

基本はmain.mを実行するだけ．この中で一般化座標の初期値や筋腱長などを定義している．main.mでは結果（一般化座標の時系列データなど）がmatファイルで保存される．
解析用のmake_video_nolta.mやstand_condition_colormap.mでは，条件ごとの結果（matファイル）をフォルダに分けて保存し，そのフォルダをmake_video_nolta.mやstand_condition_colormap.mの中で指定する．
どこで何をしているかは大体コメントつけたつもりなので，コードの中も見てください．

ode45の使い方については[ここ](https://union-of-dsc.esa.io/posts/271)を参照．

# 使い方
1. calc_EOM_work.mで運動方程式をシンボリックで計算．最下部のmatlabFunctionで必要なものを関数化している．
     ```
    matlabFunction(M,'File','Inertial_matrix', 'Vars', {m_list, l_link_list, q});
    ```
    なら，シンボリック計算したM（慣性行列）を返り値にもつ関数がInertial_matrix.mというファイル名で保存される．'Vars', {m_list, l_link_list, q}で引数を指定．

2. main.m（またはmain_init_condition.m）で一般化座標と一般化速度の初期値を設定．筋腱の自然長を設定．  
    ```
    [CFL, Ci, CFLT, GEo, GE, Pgain, Igain, Dgain] = ndgrid(L_CFL, L_Ci, L_CFLT, L_GEo, L_GE, Pgain_list, Igain_list, Dgain_list);
    length_and_gain_combination = [CFL(:), Ci(:), CFLT(:), GEo(:), GE(:), Pgain(:), Igain(:), Dgain(:)];
    ```
    で筋腱の長さとゲイン全ての組み合わせを配列に格納．L_Ciなどを配列ではなく単一の数値にすれば1つのパラメータの絞ってシミュレーション可能．  
    ```
    save('exp20240712_MuscleLengthTest_PID_length_and_gain_combination.mat','length_and_gain_combination');
    ```
    で筋腱の長さとゲインの組み合わせを保存．このファイルは後で仕事の大小などを踏まえたグラフなどを使用する．上書きしないようファイル名は要変更．複数のPCで分散してシミュレーションを行うときは，この配列を分割して実行する．  
    最下部のfor文で全ての筋腱長とゲインの組み合わせに対してシミュレーションを行う．
    ```
    simulation = solve_EOM(tmax,tspace,tspan,initial_condition,x_fixed,y_fixed,k_ground,c_ground,mu,l_link_list,l_muscle_list,limit_list,m_list,default_wire_k,default_wire_c,g,t_CFL,k_frame,c_frame,default_frame_angle,filename,end_CFL,CFL_alpha,t_end_exp,default_CFL,gain_list);
    ```
    でode45を実行しているsolve_EOM.mを呼び出し，シミュレーションを行う．
    ```
    save(filename,...);
    ```
    でシミュレーション結果を保存．

3. make_video_nolta.mなどでシミュレーション結果を解析．

中身を確認する時は，main.m → solve_EOM.mの順で見て，適宜solve_EOM.m内で使用されている関数を見るとわかりやすいと思います．ただ，matlabFunctionを使用して関数化したものは自動生成なので中身を見てもよくわかりません（Stiffness_matrix.mなどは特に．）．その場合はcalc_EOM_work.mを見て，何を導出しているのかを確認してください．

# ファイル構成
よく使うファイルの説明．基本的には以下のものでシミュレーションから解析まで可能．

## main.m
このファイルを実行するとシミュレーションが始まる．
このファイルの中でsolve_EOM.mを呼び出しており，そこで運動方程式を解いている．

## main_init_condition.m
基本はmain.mと同じだが，指定した範囲で初期姿勢を網羅的に変えながらシミュレーションを行う．

## solve_EOM.m
運動方程式を解くファイル．ode45を使用．matファイルに保存するパラメータの設定もここで．
デフォルトでほぼ全てのパラメータを保存するようになっている．
複数の関数が入れ子になっている．詳細は[ここ](https://union-of-dsc.esa.io/posts/271)を参照．ざっくりした入れ子の構成は，
```
function solve_EOM = solve_EOM(....)

    ... %保存用の配列などの定義

    function [ode, output] = getHandles()
        ode = @func_EOM;
        output = @get_data;

        function [dqdt] = func_EOM(t,q)
        ...
        end

        function status = get_data(t,q,flag)
        ...
        end

    end
    
    [ode, output] = getHandles();
    [t,q] = ode45(ode,tspan,q,options);
end
```
という感じ．`[t,q] = ode45(ode,tspan,q,options);`でode45を実行している．  `function [ode, output] = getHandles()`でodeとoutputのハンドルを取得し，ode45の引数odeは`function [dqdt] = func_EOM(t,q)`に対応．

matlabでは，大枠の関数を親関数，内部に入っている関数を子関数といい，子関数は親関数内で定義された変数にアクセスできる．そのため，`solve_EOM = solve_EOM(....)`内で色々データ保存用の変数などを定義して，子関数で使う形をとっている．

ファイル内のコメントも参照してください．

## calc_EOM_work.m
運動方程式とその計算に必要なものをシンボリックで解いているファイル．結果は関数化し，solve_EOM.mで呼び出している．
運動方程式の係数行列，踵やつま先などの座標，後で解析するための筋腱の張力や関節トルクを計算する関数などを作成．
筋腱の粘性の導入方法について，修論で記述している式と少し違う書き方をしているが，等価な式になっている．修論では各筋腱に作用する減衰力を個別に求め，それぞれを仮想仕事の原理で作用させていたが，ここでは筋腱の起始停止点の座標への力を個別に求めている．

## calc_spring_const.m
筋腱のばね定数と粘性係数を計算．デフォルトの値をmain.mで与え，この関数で筋腱が弛んでいれば0にする．

## func_GRF.m
地面反力を計算する関数．

## make_video_nolta.m
筋腱長のグラフやシミュレーションの動画を作成．calc_EOM_work.mで作成されたcalc_muscle_tension.mやcalc_torque_muscle.mを使用．詳細は本ファイルのコメントを参照．

## make_video_nolta_alt.m
基本はmake_video_nolta.mと同じだが，時間を指定してフェーズごとにグラフに背景色を追加可能．

## stand_condition_colormap.m
初期姿勢を網羅的に変えた場合に，そこから立てたもの，立てなかったもの，初期姿勢から除外されたものをカラーマップで表示．

## work_plot.m
筋腱長を網羅的に変化させた場合に，その中で仕事の小さいものや立てなかった条件などを抽出してpallarel plotを作成．仕事が少ない上位10，それ以外の立てたやつ，膝が伸展せずに立てないやつ，膝が過伸展して立てないやつのインデックスをまとめたmatファイルを作成．

## work_plot_rev.m
修論3章Fig. 3.2を作成．筋腱長を網羅的に変化させた場合に，その中で仕事の小さいものや立てなかった条件などを抽出し，仕事に応じたカラーカップでプロットの色を変えてグラフ化．仕事が少ない上位10，それ以外の立てたやつ，膝が伸展せずに立てないやつ，膝が過伸展して立てないやつのインデックスをまとめたmatファイルを作成．

## torque_comparison.m
筋腱長を網羅的に変化させた場合に，work_plot_rev.mなどで出力した立てたやつなどをまとめたmatファイルを読み込み，股関節角度とトルクの関係を色分けしてプロット．

## torque_comparison_rev.m
修論3章Fig. 3.3とFig. 3.4を作成．筋腱長を網羅的に変化させた場合に，work_plot_rev.mなどで出力した立てたやつなどをまとめたmatファイルを読み込み，股関節角度とトルクの関係を色分けしてプロット．平均値を実線で，最大値と最小値で囲まれる部分を領域としてプロット．