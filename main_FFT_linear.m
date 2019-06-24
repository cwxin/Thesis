clear all;
clc;

format long

%linear

%real equation
%第一個測試函數 u:sin(pi*x)sin(2*pi*y);
%第二個測試函數 u:x^2+y^2
%we have equation
%第一個對應到的 f=-5*pi*pi*sin(pi*x)*sin(2*pi*y);
%第二個對應到的 f=4
%因為這個程式是求delta_u=f的解，也就是求出u的數值解 
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
%因此原本要求的delta_u=f，可以轉換成：
%Bu=f
%原本利用matlab可以直接求得u=B^-1*f
%但這個程式是利用傅立葉轉換，去取代用inv函式，進而加速得到結果

%initial conditions
display('The center of elliptic equation is origin point:');
display('(x/a)^2+(y/b)^2=1, and a>b'); %橢圓方程式
a=input('a:');   %輸入橢圓長軸
b=input('b:');   %輸入橢圓短軸
phi_bounded=(atanh(b/a));   %從橢圓座標轉換後成笛卡爾座標(phi,theta)，因為我們domain是一個橢圓，因此會有一個邊界，此邊界會決定轉換後，phi的最大值
A=sqrt(a^2-b^2);            %算出常數A
N=input('Please input the odd number N of block on theta-axis(N>2):');   %輸入座標theta方向被切成幾個點
M=input('Please input the odd number M of block on phi-axis(M>2):');     %輸入座標phi方向被切幾個點
Total_point=M*N;       %計算出總共切成幾個點
fprintf('Total block:%d\n', M*N);   %表現出總共幾個點
Array=[];                           %先把矩陣Array造出來
f=zeros(Total_point,1);             %先把函數陣列f造出來
u_real=zeros(Total_point,1);        %先把真實解陣列u_real造出來
delta_N=2*pi/N;                     %theta方向一個block的長度
delta_M=phi_bounded/(M+1/2);        %phi方向一個block的長度(這個設定的方式是為了不碰觸到phi等於零的地方，原因為在座標轉換的時候會有問題)
%以下是為了要造出矩陣B所需要的係數
coefficient_phi=1/((delta_M)^2);
coefficient_theta=1/((delta_N)^2);
coefficient=(-2)*((1/(delta_M)^2)+1/((delta_N)^2));

%設立選擇控制
case_no = 1;      %選擇要測試的函數f(1 → f=-5*pi*pi*sin(pi*x)*sin(2*pi*y), 0 → f=4)
fig_output = 0;   %選擇要不要開啟2D的值(1→開啟,0→關閉)
speed_output= 0;  %選擇要不要開起profiler(1→開啟,0→關閉) 

[f, u_real, u_real_boundary, divide_matrix] = initial_setup(case_no,N,M,delta_M,delta_N,A,phi_bounded,a,b);
%這個初始設定函式可以把真實解u、初始值f、Eigenvalue matrix算出來

matrix_f=reshape(f,N,M);   %把每一個點對應好在矩陣的相對位置
f=matrix_f;                %把原本的矩陣暫存起來，因為之後還需要原本的矩陣

tic;                       %計算時間(開始)

%fft(theta)
f = FFTtheta(f,1);

%transpose
matrix_f=f';
f=f';

%dst(phi)
f=FFTphi(f,1);

%只選出是sin(Bk*phi)*sin(m*theta)的點
sin_phi_f=f;
sin_phi_f(:,(N+1)/2:1:N)=[];

%dct(phi)
f=matrix_f;
f = FFTphi(f,2);
f(:,1:(N+1)/2-1)=[];

%合併dst&dct的結果
f=[sin_phi_f,f];

%transpose
f=f';

%divide dig
f=f./divide_matrix;

%transpose
matrix_f=f';
f=f';

%idst(phi)
f=FFTphi(f,-1);

%只選出是sin(Bk*phi)*sin(m*theta)的點
sin_phi_f=f;
sin_phi_f(:,(N+1)/2:1:N)=[];

%idct(phi)
f=matrix_f;
f = FFTphi(f,-2);
f(:,1:(N+1)/2-1)=[];

%合併idst&idct的結果
f=[sin_phi_f,f];

%transpose
f=f';

%ifft(theta)
f = FFTtheta(f,-1);

time=toc;     %計算時間(結束)

%計算數值解與真實解的最大誤差
u_num=reshape(f,N*M,1);
error=u_real-u_num;
max_error=max(abs(error));

fprintf('Max error:%e\n', max_error);            %顯示最大誤差
fprintf('Time of (for loop):%f sec\n', time);    %顯示所花時間

%顯示各細項所花時間
if speed_output==1
    profile viewer	
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
    u_num_2d = reshape(u_num,N,M);
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