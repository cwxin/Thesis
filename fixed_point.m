function [u_k] = fixed_point(f,u_k_1,divide_matrix,N,M,c) 

        tilde_f=f+c*sinh(u_k_1);
        u_k=Fastpoisson(tilde_f,divide_matrix,N,M);
        
end