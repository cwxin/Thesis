function [u]=Fastpoisson(f,divide_matrix,N,M)
    
            matrix_f=reshape(f,N,M);    %��C�@���I�����n�b�x�}���۹��m
            f=matrix_f;                 %��쥻���x�}�Ȧs�_�ӡA�]�������ٻݭn�쥻���x�}

            %fft(theta)
            f = FFTtheta(f,1);

            %transpose�A�H�K�B�zphi��V
            matrix_f=f';
            f=f';

            %dst(phi)
            f=FFTphi(f,1);

            %�u��X�Osin(Bk*phi)*sin(m*theta)���I
            sin_phi_f=f;
            sin_phi_f(:,(N+1)/2:1:N)=[];

            %�n��X�Ocos(Ak*phi)*cos(m*theta)���I
            %dct(phi)
            f=matrix_f;
            f = FFTphi(f,2);
            f(:,1:(N+1)/2-1)=[];

            %�X��dst(phi)&dct(phi)�����G
            f=[sin_phi_f,f];

            %transpose
            f=f';

            %divide dig
            f=f./divide_matrix;

            %transpose�A�H�K�B�zphi��V
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

            %�X��idst(phi)&idct(phi)�����G
            f=[sin_phi_f,f];

            %transpose
            f=f';

            %ifft(theta)
            u = FFTtheta(f,-1);
        
end