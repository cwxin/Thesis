clear all;
clc;

format long

%linear

%real equation
%�Ĥ@�Ӵ��ը�� u:sin(pi*x)sin(2*pi*y);
%�ĤG�Ӵ��ը�� u:x^2+y^2
%we have equation
%�Ĥ@�ӹ����쪺 f=-5*pi*pi*sin(pi*x)*sin(2*pi*y);
%�ĤG�ӹ����쪺 f=4
%�]���o�ӵ{���O�Ddelta_u=f���ѡA�]�N�O�D�Xu���ƭȸ� 
%(x/a)^2+(y/b)^2=1, and a>b ����{��(a=���b,b=�u�b)
%���y�Шt���G
%x=A*cosh(phi)*cos(theta)
%y=A*sinh(phi)*sin(theta)
%�bx-y�y�ФU�A���i�H�Q��ܦ�
%x=a*cos(theta)
%y=b*sin(theta)
%�]���i�H���Da=A*cosh(phi), b=A*sinh(phi)
%����i�H�����X(a^2-b^2)=A*(cosh^2(phi)-sinh^2(phi))=A
%�o�̨ϥ�Central-Difference Formulas
%�D�X�ڭ̩Ҥ��ΥX�Ӫ��I���ƭȸ�
%�]���쥻�n�D��delta_u=f�A�i�H�ഫ���G
%Bu=f
%�쥻�Q��matlab�i�H�����D�ou=B^-1*f
%���o�ӵ{���O�Q�γť߸��ഫ�A�h���N��inv�禡�A�i�ӥ[�t�o�쵲�G

%initial conditions
display('The center of elliptic equation is origin point:');
display('(x/a)^2+(y/b)^2=1, and a>b'); %����{��
a=input('a:');   %��J�����b
b=input('b:');   %��J���u�b
phi_bounded=(atanh(b/a));   %�q���y���ഫ�ᦨ�åd���y��(phi,theta)�A�]���ڭ�domain�O�@�Ӿ��A�]���|���@����ɡA����ɷ|�M�w�ഫ��Aphi���̤j��
A=sqrt(a^2-b^2);            %��X�`��A
N=input('Please input the odd number N of block on theta-axis(N>2):');   %��J�y��theta��V�Q�����X���I
M=input('Please input the odd number M of block on phi-axis(M>2):');     %��J�y��phi��V�Q���X���I
Total_point=M*N;       %�p��X�`�@�����X���I
fprintf('Total block:%d\n', M*N);   %��{�X�`�@�X���I
Array=[];                           %����x�}Array�y�X��
f=zeros(Total_point,1);             %�����ư}�Cf�y�X��
u_real=zeros(Total_point,1);        %����u��Ѱ}�Cu_real�y�X��
delta_N=2*pi/N;                     %theta��V�@��block������
delta_M=phi_bounded/(M+1/2);        %phi��V�@��block������(�o�ӳ]�w���覡�O���F���IĲ��phi����s���a��A��]���b�y���ഫ���ɭԷ|�����D)
%�H�U�O���F�n�y�X�x�}B�һݭn���Y��
coefficient_phi=1/((delta_M)^2);
coefficient_theta=1/((delta_N)^2);
coefficient=(-2)*((1/(delta_M)^2)+1/((delta_N)^2));

%�]�߿�ܱ���
case_no = 1;      %��ܭn���ժ����f(1 �� f=-5*pi*pi*sin(pi*x)*sin(2*pi*y), 0 �� f=4)
fig_output = 0;   %��ܭn���n�}��2D����(1���}��,0������)
speed_output= 0;  %��ܭn���n�}�_profiler(1���}��,0������) 

[f, u_real, u_real_boundary, divide_matrix] = initial_setup(case_no,N,M,delta_M,delta_N,A,phi_bounded,a,b);
%�o�Ӫ�l�]�w�禡�i�H��u���u�B��l��f�BEigenvalue matrix��X��

matrix_f=reshape(f,N,M);   %��C�@���I�����n�b�x�}���۹��m
f=matrix_f;                %��쥻���x�}�Ȧs�_�ӡA�]�������ٻݭn�쥻���x�}

tic;                       %�p��ɶ�(�}�l)

%fft(theta)
f = FFTtheta(f,1);

%transpose
matrix_f=f';
f=f';

%dst(phi)
f=FFTphi(f,1);

%�u��X�Osin(Bk*phi)*sin(m*theta)���I
sin_phi_f=f;
sin_phi_f(:,(N+1)/2:1:N)=[];

%dct(phi)
f=matrix_f;
f = FFTphi(f,2);
f(:,1:(N+1)/2-1)=[];

%�X��dst&dct�����G
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

%�u��X�Osin(Bk*phi)*sin(m*theta)���I
sin_phi_f=f;
sin_phi_f(:,(N+1)/2:1:N)=[];

%idct(phi)
f=matrix_f;
f = FFTphi(f,-2);
f(:,1:(N+1)/2-1)=[];

%�X��idst&idct�����G
f=[sin_phi_f,f];

%transpose
f=f';

%ifft(theta)
f = FFTtheta(f,-1);

time=toc;     %�p��ɶ�(����)

%�p��ƭȸѻP�u��Ѫ��̤j�~�t
u_num=reshape(f,N*M,1);
error=u_real-u_num;
max_error=max(abs(error));

fprintf('Max error:%e\n', max_error);            %��̤ܳj�~�t
fprintf('Time of (for loop):%f sec\n', time);    %��ܩҪ�ɶ�

%��ܦU�Ӷ��Ҫ�ɶ�
if speed_output==1
    profile viewer	
end

%figure(��ܯu��ѻP�ƭȸѪ�2D���ϧ�)
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