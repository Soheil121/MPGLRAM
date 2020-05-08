function psnr = PSNR(L_, M, R_, A, n,k,msize)
%  SUBROUTINE TO COMPUTE THE RMSRE VALUE
     M_PSNR=msize;
     N_PSNR=msize;
A_const=cell(1,n);
for  i=1:n
    A_const{i}=0;
end
for i =1:n
    for K=1:k
        A_const{i}=A_const{i}+L_{K}*M{i}*R_{K}';
    end
end

MSE = sum(sum((A{i}-A_const{i}).^2))/(M_PSNR*N_PSNR);
    psnr = 10*log10(256*256/MSE);
%     fprintf('\nMSE: %7.2f ', MSE);
%     fprintf('\nPSNR: %9.7f dB', psnr);

