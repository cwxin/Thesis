function [u]=Fastpoisson(f,divide_matrix,N,M)
    
            matrix_f=reshape(f,N,M);    %把每一個點對應好在矩陣的相對位置
            f=matrix_f;                 %把原本的矩陣暫存起來，因為之後還需要原本的矩陣

            %fft(theta)
            f = FFTtheta(f,1);

            %transpose，以便處理phi方向
            matrix_f=f';
            f=f';

            %dst(phi)
            f=FFTphi(f,1);

            %只選出是sin(Bk*phi)*sin(m*theta)的點
            sin_phi_f=f;
            sin_phi_f(:,(N+1)/2:1:N)=[];

            %要選出是cos(Ak*phi)*cos(m*theta)的點
            %dct(phi)
            f=matrix_f;
            f = FFTphi(f,2);
            f(:,1:(N+1)/2-1)=[];

            %合併dst(phi)&dct(phi)的結果
            f=[sin_phi_f,f];

            %transpose
            f=f';

            %divide dig
            f=f./divide_matrix;

            %transpose，以便處理phi方向
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

            %合併idst(phi)&idct(phi)的結果
            f=[sin_phi_f,f];

            %transpose
            f=f';

            %ifft(theta)
            u = FFTtheta(f,-1);
        
end