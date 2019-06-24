function P = FFTtheta(F, s)

P = F;
[N,M] = size(F);

if s == 1
    
    %fft(theta)
    fm=zeros(2*N,M);
    for i=1:N
        fm(i*2-1,:)=F(i,:);
    end
    fm=fft(fm);

    %����Hsin(Bk*phi)*sin(m*theta)��ܪ��I
    j=0;
    for i=N-(N-1)/2:N-1
        j=1+j;
        P(j,:)=-imag(fm(i+1,:));
    end

    %����Hcos(Ak*phi)*cos(m*theta)��ܪ��I
    for i=N-(N-1)/2:N
        P(i,:)=real(fm(i+1,:));
    end
    
elseif s == -1
    
    %ifft(theta)(���bfft�����A�N�b����)
    f_real=zeros(N,M);
    f_imag=zeros(N,M);
    f_real(1,:)=F(N,:);
    for i=1:(N-1)/2
        f_real(i+1,:)=F(N-i,:);
        f_imag(i+1,:)=F((N+1)/2-i,:);
    end
    j=1;
    for i=(N-1)/2+1:N-1
        f_real(i+1,:)=F(i,:);
        f_imag(i+1,:)=-F(j,:);
        j=j+1;
    end
    f_real=kron([1;1],f_real);
    f_imag=kron([1;1],f_imag);
    i=sqrt(-1);
    f_imag=f_imag.*i;
    fm=f_real+f_imag;
    fm=ifft(fm);
    for i=1:N
        P(i,:)=fm(i*2-1,:);
    end 
    
else
    
    disp('s=1(FFT) or -1(iFFT)');
    
end
        