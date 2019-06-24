function [f]=eliminate_intial(case_no,N,M,delta_M,delta_N,A,a,b)

switch case_no
    case 1
        %make initial f
        %u = sin(pi*x)*sin(2*pi*y)

        %�o�O�令�V�q�x�}����k�h���B�z�����
        x1=0:(M-1);
        y1=ones(1,N);
        z1=kron(x1,y1);
        x2=0:(N-1);
        y2=ones(1,M);
        z2=kron(y2,x2);
        u_real=(sin(pi.*A.*cosh(delta_M/2+z1.*delta_M).*cos(z2.*delta_N))).*(sin(2.*pi.*A.*sinh(delta_M/2+z1.*delta_M).*sin(z2.*delta_N)));
        u_real=u_real';
        denominator=A.*A.*((sinh(delta_M/2+z1.*delta_M)).^2+(sin(z2.*delta_N)).^2);
        f=(-5.0).*(pi*pi).*(sin(pi.*A.*(cosh(delta_M/2+z1.*delta_M).*cos(z2.*delta_N))).*sin(2.*pi.*A.*(sinh(delta_M/2+z1.*delta_M).*sin(z2.*delta_N)))).*(denominator);
        f=f';        
        
        %�B�z�o��O�]�����p�ɻݭn���Dphi=phi_bounded���ȡA���ڭ̭n�D���ƭȸѨä��|�]�t��ɡA�ҥH���e��ݭn��ɦ��p���ȥ�����
        i = M;
        u_real_boundary=zeros(N,1);
        for j=1:N
            u_real_boundary(j) = (sin(pi*A*cosh(delta_M/2+i*delta_M)*cos((j-1)*delta_N)))*(sin(2*pi*A*sinh(delta_M/2+i*delta_M)*sin((j-1)*delta_N)));
            f(N*M-N+j,1)=f(N*M-N+j,1)-1/(delta_M)^2*u_real_boundary(j);
        end
        
    otherwise
        %make initial f
        % u = x^2+y^2
        
        %�o�O�令�V�q�x�}����k�h���B�z�����
        x1=0:(M-1);
        y1=ones(1,N);
        z1=kron(x1,y1);
        x2=0:(N-1);
        y2=ones(1,M);
        z2=kron(y2,x2);
        u_real=(A.*cosh(delta_M/2+z1.*delta_M).*cos(z2.*delta_N)).^2.+(A.*sinh(delta_M/2+z1.*delta_M).*sin(z2.*delta_N)).^2;
        u_real=u_real'./(a^2+b^2);
        denominator=A.*A.*((sinh(delta_M/2+z1.*delta_M)).^2+(sin(z2.*delta_N)).^2);
        f=4*(denominator)./(a^2+b^2);
        f=f';
        
        %�B�z�o��O�]�����p�ɻݭn���Dphi=phi_bounded���ȡA���ڭ̭n�D���ƭȸѨä��|�]�t��ɡA�ҥH���e��ݭn��ɦ��p���ȥ�����
        i = M;
        u_real_boundary=zeros(N,1);
        for j=1:N
            u_real_boundary(j) = (A*cosh(delta_M/2+i*delta_M)*cos((j-1)*delta_N))^2+(A*sinh(delta_M/2+i*delta_M)*sin((j-1)*delta_N))^2;
            f(N*M-N+j,1)=f(N*M-N+j,1)-1/(delta_M)^2*u_real_boundary(j)./(a^2+b^2);
        end      
        
end