function error = ERROR(L_, M, R_, A, n,k)
%  SUBROUTINE TO COMPUTE THE RMSRE VALUE

error = 0;
A_const=cell(1,n);
for  i=1:n
    A_const{i}=0;
end
for i =1:n
    for K=1:k
        A_const{i}=A_const{i}+L_{K}*M{i}*R_{K}';
    end
    TEMPT       = A{i} - A_const{i};
    error      = error + norm(TEMPT,'fro')^2;
end

error = sqrt(error/n);
