function [f]=eliminate_intial(case_no,N,M,delta_M,delta_N,A,a,b)

switch case_no
    case 1
        %make initial f
        %u = sin(pi*x)*sin(2*pi*y)

        %這是改成向量矩陣的方法去做處理比較快
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
        
        %處理這邊是因為估計時需要取道phi=phi_bounded的值，但我們要求的數值解並不會包含邊界，所以提前把需要邊界估計的值先扣掉
        i = M;
        u_real_boundary=zeros(N,1);
        for j=1:N
            u_real_boundary(j) = (sin(pi*A*cosh(delta_M/2+i*delta_M)*cos((j-1)*delta_N)))*(sin(2*pi*A*sinh(delta_M/2+i*delta_M)*sin((j-1)*delta_N)));
            f(N*M-N+j,1)=f(N*M-N+j,1)-1/(delta_M)^2*u_real_boundary(j);
        end
        
    otherwise
        %make initial f
        % u = x^2+y^2
        
        %這是改成向量矩陣的方法去做處理比較快
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
        
        %處理這邊是因為估計時需要取道phi=phi_bounded的值，但我們要求的數值解並不會包含邊界，所以提前把需要邊界估計的值先扣掉
        i = M;
        u_real_boundary=zeros(N,1);
        for j=1:N
            u_real_boundary(j) = (A*cosh(delta_M/2+i*delta_M)*cos((j-1)*delta_N))^2+(A*sinh(delta_M/2+i*delta_M)*sin((j-1)*delta_N))^2;
            f(N*M-N+j,1)=f(N*M-N+j,1)-1/(delta_M)^2*u_real_boundary(j)./(a^2+b^2);
        end      
        
end