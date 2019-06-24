clear all;
clc;

format long

%nonlinear

%real equation
%�Ĥ@�Ӵ��ը�� u:sin(pi*x)sin(2*pi*y);
%�ĤG�Ӵ��ը�� u:x^2+y^2
%we have equation
%�Ĥ@�ӹ����쪺 f=-5*pi*pi*sin(pi*x)*sin(2*pi*y);
%�ĤG�ӹ����쪺 f=4
%�]���o�ӵ{���O�Ddelta_u-sinh(u)=f���ѡA�]�N�O�D�Xu���ƭȸ� 
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
%�]���쥻�n�D��delta_u+sinh(u)=f�A�i�H�ഫ���G
%Bu=f+sinh(u)
%�쥻�Q��matlab�i�H�����D�ou=B^-1*(f-sinh(u))
%�ӳo�ӵ{���O�Q�έ��N�k�ӨD�ou���ƭȸ�
%�L�{���Q�γť߸��ഫ�A�h���N��inv�禡�A�i�ӥ[�t�o�쵲�G

display('The center of elliptic equation is origin point:');
display('(x/a)^2+(y/b)^2=1, and a>b');  %����{��
a=input('a:');     %��J�����b
b=input('b:');     %��J���u�b
phi_bounded=(atanh(b/a));   %�q���y���ഫ�ᦨ�åd���y��(phi,theta)�A�]���ڭ�domain�O�@�Ӿ��A�]���|���@����ɡA����ɷ|�M�w�ഫ��Aphi���̤j��
A=sqrt(a^2-b^2);            %��X�`��A
N=input('Please input the odd number N of block on theta-axis(N>2):');  %��J�y��theta��V�Q�����X���I
M=input('Please input the odd number M of block on phi-axis(M>2):');    %��J�y��phi��V�Q���X���I
Total_point=M*N;        %�p��X�`�@�����X���I
fprintf('Total block:%d\n', M*N);   %��{�X�`�@�X���I
delta_N=2*pi/N;                   %theta��V�@��block������
delta_M=phi_bounded/(M+1/2);      %phi��V�@��block������(�o�ӳ]�w���覡�O���F���IĲ��phi����s���a��A��]���b�y���ഫ���ɭԷ|�����D)
%�H�U�O���F�n�y�X�x�}B�һݭn���Y��
coefficient_phi=1/((delta_M)^2);  
coefficient_theta=1/((delta_N)^2);
coefficient=(-2)*((1/(delta_M)^2)+1/((delta_N)^2));

%�]�߿�ܱ���
case_no = 1;      %��ܭn���ժ����f(1 �� f=-5*pi*pi*sin(pi*x)*sin(2*pi*y), 0 �� f=4)
case_iterate = 1; %��ܭ��N��k(1 �� �T�w�I���N, 0 �� ���y�k���N)
fig_output = 0;   %��ܭn���n�}��2D�Ϥ�����(1���}��,0������)
speed_output= 0;  %��ܭn���n�}�_profiler(1���}��,0������) 
c=1;              %c*sinh�e���Y��
[f, u_real, divide_matrix, u_real_boundary] = nonlinear_initial_setup(case_no,N,M,delta_M,delta_N,A,phi_bounded,c,a,b); 
%�o�Ӫ�l�]�w�禡�i�H��u���u�B��l��f�BEigenvalue matrix��X�ӡA���F�e�X�ϧΡA�ݭn�u��Ѫ����

[u_k_1] = eliminate_intial(case_no,N,M,delta_M,delta_N,A,a,b);
u_k_1=Fastpoisson(u_k_1,divide_matrix,N,M); 
u_k_1=reshape(u_k_1,N*M,1);                      
%��X���N���Ĥ@�Ӫ�l�ȡA�令u_k_1=(A^-1*f)

epsilon=1;          %�w�]���N�~�t�ϥH�U�禡�i�H�B�@
iterate_number=0;   %�w�]���N���Ƭ�0

%��ܭn�ϥέ��@�ح��N�覡
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
        [Array]=Make_matrix(N,M,coefficient,coefficient_theta,coefficient_phi);   %�s�yA�x�}
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

%�p��u��ѻP�ƭȸѮt�̤j����
error=u_k-u_real;
error=max(abs(error));
%�p��u��ѻP�ƭȸѪ��~�t2-norm
two_norm=norm(u_k-u_real,2);

fprintf('Error( max(abs(u_real-u_k)) ) : %e \n', error);   %��̤ܳj�~�t
fprintf('Two_norm( u_real and u_k ) : %e \n', two_norm);   %���2-norm
fprintf('Iterate_number : %d \n', iterate_number);         %��ܭ��N����
fprintf('Iterative time : %f sec\n', time2);               %��ܩҪ�ɶ�

%�s�y�x�}Array���ɶ�
%time1

%��ܦU�Ӷ��Ҫ�ɶ�
if speed_output==1
    profile viewer;	
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