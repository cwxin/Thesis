function [A]= Make_matrix(N,M,coefficient,coefficient_theta,coefficient_phi)
    %Make Array

    %Sparse
            
    %Processing Ai
    Total_point=N*M;
    Xi=1:Total_point;
    temp=ones(1,6)'*Xi;
    temp1=kron(ones(1,M),[2:N+1]);
    temp(6,:)=temp1;
    temp(3,:)=temp(3,:)+1;
    temp(5,1:N)=temp(5,1:N)-Total_point;
    temp(4,Total_point-N+1:Total_point)=temp(4,Total_point-N+1:Total_point)-Total_point;
    for i=1:M
        temp(3,i*N)=temp(3,i*N)-N;
    end
    Ai=reshape(temp,1,6*Total_point);
    %Processing Aj
    temp=ones(1,6)'*Xi;
    temp(2,:)=temp(2,:)+1;
    temp(4,:)=temp(4,:)+N;
    temp(5,:)=temp(5,:)-N;
    for i=1:M
        temp(2,i*N)=temp(2,i*N)-N;
    end
    temp1=kron(ones(1,M),[N:-1:1]);
    temp(6,:)=temp1;
     Aj=reshape(temp,1,6*Total_point);
    %Processing Av
    Yi=[coefficient,coefficient_theta,coefficient_theta,coefficient_phi,coefficient_phi,(coefficient_phi)/N];
    temp=(ones(1,Total_point)'*Yi)';
    Av=reshape(temp,1,6*Total_point);

    NID = (Aj <= 0 | Aj > Total_point);
    Ai(NID) = [];
    Aj(NID) = [];
    Av(NID) = [];
    Array = sparse(Ai,Aj,Av,Total_point,Total_point);

    %ideal with special points
    if rem(N,2)==0
        Array((N/2)+1,(N/2)+1)=coefficient+(coefficient_phi);
        Array(1,1)=coefficient+(coefficient_phi);
        Array(N+1,1)=coefficient_phi;
    else
        Array(1,1)=coefficient+(coefficient_phi);
        Array((N+1)/2,(N+1)/2+1)=(coefficient_theta)+(coefficient_phi);
        Array((N+1)/2+1,(N+1)/2)=(coefficient_theta)+(coefficient_phi);
        Array(N+1,1)=coefficient_phi;
    end
    
    A=Array;
    
end