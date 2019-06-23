function [u]=Newton_method(f,u,divide_matrix,N,M,A,c)
    
    %解(A+D)x=F=A*u+sinh-f → x_k=(A^-1)*(F-D*x_k_1) 
    x_k_1=zeros(N*M,1);    %給初始的x_k_1
    A_u=A*u;
    sin_h=c*sinh(u);       %算出sinh(u)
    F=A_u-sin_h-f;         %算出F
    epsilon=1;             %預設迭代誤差使以下函式可以運作
    %固定點迭代
    while epsilon>=10^-8 && isnan(epsilon)==0 && isinf(epsilon)==0
        Di=c*cosh(u).*x_k_1;                   %算出D*x_k_1
        b=F+Di;                                %F+D
        x_k=Fastpoisson(b,divide_matrix,N,M);  %算出(A^-1)*(F-D*x_k_1)
        x_k=reshape(x_k,N*M,1);
        error=x_k-x_k_1;
        epsilon=max(abs(error));
        x_k_1=x_k;
        %pause
    end
    u=x_k;
    
end