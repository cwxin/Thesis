clear all;
clc;

format long

%nonlinear

%real equation
%第一個測試函數 u:sin(pi*x)sin(2*pi*y);
%第二個測試函數 u:x^2+y^2
%we have equation
%第一個對應到的 f=-5*pi*pi*sin(pi*x)*sin(2*pi*y);
%第二個對應到的 f=4
%因為這個程式是求delta_u-sinh(u)=f的解，也就是求出u的數值解 
%(x/a)^2+(y/b)^2=1, and a>b 橢圓方程式(a=長軸,b=短軸)
%橢圓座標系為：
%x=A*cosh(phi)*cos(theta)
%y=A*sinh(phi)*sin(theta)
%在x-y座標下，橢圓可以被表示成
%x=a*cos(theta)
%y=b*sin(theta)
%因此可以知道a=A*cosh(phi), b=A*sinh(phi)
%那麼可以推測出(a^2-b^2)=A*(cosh^2(phi)-sinh^2(phi))=A
%這裡使用Central-Difference Formulas
%求出我們所切割出來的點的數值解
%因此原本要求的delta_u+sinh(u)=f，可以轉換成：
%Bu=f+sinh(u)
%原本利用matlab可以直接求得u=B^-1*(f-sinh(u))
%而這個程式是利用迭代法來求得u的數值解
%過程中利用傅立葉轉換，去取代用inv函式，進而加速得到結果

display('The center of elliptic equation is origin point:');
display('(x/a)^2+(y/b)^2=1, and a>b');  %橢圓方程式
a=input('a:');     %輸入橢圓長軸
b=input('b:');     %輸入橢圓短軸
phi_bounded=(atanh(b/a));   %從橢圓座標轉換後成笛卡爾座標(phi,theta)，因為我們domain是一個橢圓，因此會有一個邊界，此邊界會決定轉換後，phi的最大值
A=sqrt(a^2-b^2);            %算出常數A
N=input('Please input the odd number N of block on theta-axis(N>2):');  %輸入座標theta方向被切成幾個點
M=input('Please input the odd number M of block on phi-axis(M>2):');    %輸入座標phi方向被切幾個點
Total_point=M*N;        %計算出總共切成幾個點
fprintf('Total block:%d\n', M*N);   %表現出總共幾個點
delta_N=2*pi/N;                   %theta方向一個block的長度
delta_M=phi_bounded/(M+1/2);      %phi方向一個block的長度(這個設定的方式是為了不碰觸到phi等於零的地方，原因為在座標轉換的時候會有問題)
%以下是為了要造出矩陣B所需要的係數
coefficient_phi=1/((delta_M)^2);  
coefficient_theta=1/((delta_N)^2);
coefficient=(-2)*((1/(delta_M)^2)+1/((delta_N)^2));

%設立選擇控制
case_no = 1;      %選擇要測試的函數f(1 → f=-5*pi*pi*sin(pi*x)*sin(2*pi*y), 0 → f=4)
case_iterate = 1; %選擇迭代方法(1 → 固定點迭代, 0 → 牛頓法迭代)
fig_output = 0;   %選擇要不要開啟2D圖片的值(1→開啟,0→關閉)
speed_output= 0;  %選擇要不要開起profiler(1→開啟,0→關閉) 
c=1;              %c*sinh前的係數
[f, u_real, divide_matrix, u_real_boundary] = nonlinear_initial_setup(case_no,N,M,delta_M,delta_N,A,phi_bounded,c,a,b); 
%這個初始設定函式可以把真實解u、初始值f、Eigenvalue matrix算出來，為了畫出圖形，需要真實解的邊界

[u_k_1] = eliminate_intial(case_no,N,M,delta_M,delta_N,A,a,b);
u_k_1=Fastpoisson(u_k_1,divide_matrix,N,M); 
u_k_1=reshape(u_k_1,N*M,1);                      
%找出迭代的第一個初始值，改成u_k_1=(A^-1*f)

epsilon=1;          %預設迭代誤差使以下函式可以運作
iterate_number=0;   %預設迭代次數為0

%選擇要使用哪一種迭代方式
switch case_iterate
    %Fixed Point Iterate Method 
    %u_k=A^-1*(f-sinh(u_k_1))
    case 1
        tic;
        while epsilon>=10^-8 && isnan(epsilon)==0 && isinf(epsilon)==0
            [u_k]=fixed_point(f,u_k_1,divide_matrix,N,M,c);
            u_k=reshape(u_k,N*M,1);
            max(abs(u_k));
            error=u_k-u_k_1;
            epsilon=max(abs(error));
            u_k_1=u_k;
            iterate_number=iterate_number+1;
        end
        time2=toc;
    otherwise
        %Newton's Method
        %u_k = u_k_1-((A+cosh(u_k_1))^-1)*F = u_k_1-u
        %F=A*u_k_1+sinh(u_k_1)-f
        tic;
        [Array]=Make_matrix(N,M,coefficient,coefficient_theta,coefficient_phi);   %製造A矩陣
        time1=toc;
        tic;
        while epsilon>=10^-8 && isnan(epsilon)==0 && isinf(epsilon)==0
            [u]=Newton_method(f,u_k_1,divide_matrix,N,M,Array,c);
            max(abs(u));
            max(abs(u_k_1));
            u_k=u_k_1-u;
            A_u=Array*u_k;
            sin_h=c*sinh(u_k);
            F=A_u-sin_h-f;
            epsilon=max(abs(F));
            u_k_1=u_k;
            iterate_number=iterate_number+1;
        end
        time2=toc;
end

%計算真實解與數值解差最大的值
error=u_k-u_real;
error=max(abs(error));
%計算真實解與數值解的誤差2-norm
two_norm=norm(u_k-u_real,2);

fprintf('Error( max(abs(u_real-u_k)) ) : %e \n', error);   %顯示最大誤差
fprintf('Two_norm( u_real and u_k ) : %e \n', two_norm);   %顯示2-norm
fprintf('Iterate_number : %d \n', iterate_number);         %顯示迭代次數
fprintf('Iterative time : %f sec\n', time2);               %顯示所花時間

%製造矩陣Array的時間
%time1

%顯示各細項所花時間
if speed_output==1
    profile viewer;	
end

%figure(顯示真實解與數值解的2D比對圖形)
if fig_output==1
    theta = 0:delta_N:(2*pi);%-delta_N);
    phi   = delta_M/2:delta_M:phi_bounded;%-delta_M;
    [PHI,THETA] = meshgrid(phi,theta);
    X = A*cosh(PHI).*cos(THETA);
    Y = A*sinh(PHI).*sin(THETA); 
    u_real_2d = reshape(u_real,N,M);
    u_real_2d = [u_real_2d, u_real_boundary(:)];
    u_real_2d = [u_real_2d; u_real_2d(1,:)];
    u_num_2d = reshape(u_k,N,M);
    u_num_2d =[u_num_2d, u_real_boundary(:)];
    u_num_2d = [u_num_2d; u_num_2d(1,:)];
    
    figure(1); 
    subplot(1,2,1),
        surf(X,Y,u_real_2d);
        shading interp
        view(2);
        colorbar;
        title('u-real');
    subplot(1,2,2)
        surf(X,Y,u_num_2d);
        shading interp
        view(2);
        colorbar;
        title('u-num');
end