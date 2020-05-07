function [L,D,R,M]=MPGLRAM__(original_data, DATA,L,D,R,n,k,d,ite)
%% input:
%original_data(m*m,n) m*m:features n:samples.
%DATA{ij}: each cell contains an (m,m) image.
%L{j},R{j} intitalized by L0(1:d,1:d)=eye(d)and R0(1:d,1:d)=eye(d) respectively.
%D{i}: each cell contains a (d,d) reduced data from given by initilizing stage.
%n: the number of samples.
%k: the number of transformation pairs.
%d: the desired reduced dimension.
%ite: the number of the function's iteration.
%% output:
%L{j},R{j} intitalized by L0(1:d,1:d)=eye(d)and R0(1:d,1:d)=eye(d) respectively.
%D{i}: each cell contains a (d,d) reduced data from given by initilizing stage.
%M (d*d,n): the final reduced data 
%% Proposed method

for ITE=1:ite
    Sigma_kron=0;
    for K=1:k
        Sigma_A_TMM=0;Sigma_MMTMM=0;Sigma_A_MM2=0;Sigma_MM2TMM2=0;
        for i=1:n
            Sigma_LDR=0;
            for jj=1:k
                Sigma_LDR=Sigma_LDR+L{jj}*D{i}*R{jj}';
            end
            A_bar{i}=DATA{i}-(Sigma_LDR-L{K}*D{i}*R{K}');
        end
        for i=1:n % fix L, calculate R
            MM{i}=L{K}*D{i};
            Sigma_A_TMM=Sigma_A_TMM+A_bar{i}'*MM{i};
            Sigma_MMTMM=Sigma_MMTMM+MM{i}'*MM{i};
        end
        R{K}=Sigma_A_TMM/Sigma_MMTMM;
        
        for i=1:n % fix R, calculate L
            MM2{i}=R{K}*D{i}';
            Sigma_A_MM2=Sigma_A_MM2+A_bar{i}*MM2{i};
            Sigma_MM2TMM2=Sigma_MM2TMM2+MM2{i}'*MM2{i};
        end
        L{K}=Sigma_A_MM2/Sigma_MM2TMM2;
    end
    
    for KK=1:k
        Sigma_kron=Sigma_kron+kron(R{KK},L{KK});
    end
    M=Sigma_kron\original_data ;
    for ij=1:n
        D{ij}=reshape(M(:,ij),d,d);
    end
end