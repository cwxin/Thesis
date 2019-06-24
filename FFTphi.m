function P = FFTphi(F, s)

[M,N] = size(F);

if s == 1
    
    % dst(phi)
    fm=zeros((M*2),N);
    for i=1:N
        fm1=fm(:,i);
        fm2=F(:,i);
        fm1(1:2:2*M)=fm2;
        fm(:,i)=fm1;
    end
    fm=dst(fm);
    P=zeros(M,N);
    for i=1:M
        P(i,:)=fm(M-i+1,:);
    end
    
elseif s == -1
    
    % idst(phi)
    F=flipud(F);
    fm=[F;zeros(M,N)];
    fm=dst(fm);
    P=zeros(M,N);
    for i=1:M
        P(i,:)=fm(i*2-1,:);
    end
    P=P.*4/(M*2+1);
    
elseif s == 2
    
    % dct(phi)
    
    %{
    這個需要比較大的空間
    fm=zeros(2*M+2,N);
    for i=1:M
        fm(2*i,:)=F(i,:);
    end
    fm1=-flipud(fm);
    fm1(1,:)=[];
    fm=[fm;fm1];
    fm=[fm;zeros((2*M+1)*2-1,N)];
    fm=real(fft(fm))./2;
    P=zeros(M,N);
    for i=1:M
        P(i,:)=fm((M-i+1)*2,:);
    end
    %}
    
    fm=zeros(M+1,N);
    fm1=[F;fm];
    fm1=(dct(fm1))./(2/(M*2+1))^0.5;
    P=zeros(M,N);
    for i=1:M
        P(i,:)=fm1(i*2,:);
    end
    P=flipud(P);
       
elseif s == -2
    
    % idct(phi)
    
    %{
    這個需要比較大的空間
    f=flipud(F);
    fm=zeros(2*M+2,N);
    for i=1:M
        fm(2*i,:)=f(i,:);
    end
    fm1=-flipud(fm);
    fm1(1,:)=[];
    fm=[fm;fm1];
    fm=[fm;zeros((2*M+1)*2-1,N)];
    fm=real(fft(fm))./2;
    P=zeros(M,N);
    for i=1:M
        P(i,:)=fm((M-i+1)*2,:);
    end
    P=P*4/(M*2+1);
    P=flipud(P);
    %}
    
    F=flipud(F);
    fm=zeros(M+1,N);
    fm1=[F;fm];
    fm1=fm1.*(2/(2*M+1))^0.5;
    fm1=(dct(fm1)).*2;
    P=zeros(M,N);
    for i=1:M
        P(i,:)=fm1(i*2,:);
    end
    
end