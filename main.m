close all;
clear;
clc;
%% Data
load()% load the data
% We assume that our data stored in "fea". 
%fea=your data name;
n=size(fea,1);
m=sqrt(size(fea,2));
original_data=fea'; % m*m:features n:samples

%% Parameters

ite=; % the number of MPGLRAM's iteration
K_=; % k-pairs of transformations
d=; % the desired reduced dimension
%%
DATA=cell(1,n); % input matrix to n separate matrices with m*m dimension

% Mat2Img
for ij=1:n
    DATA{ij}=reshape(original_data(:,ij),m,m);
end
%% MPGLRAM


%% initialize

D   =cell(1,n); % reduced matrix to n separate matrices with d*d dimension
L   =cell(K_,1); % left transformations matrices m*d
R   =cell(K_,1); % right transformations matrices m*d
Sigma_kron   =0; % sum of kronecker product of R and L
% L0
L0=zeros(m,d);
L0(1:d,1:d)=eye(d);
% R0
R0=zeros(m,d);
R0(1:d,1:d)=eye(d);

%% The first level of multiple k-pairs of transportaion
tic
for j=1:K_
    
    L{j}=L0;
    R{j}=R0;
    
    Sigma_kron=Sigma_kron+kron(R{j},L{j});
end
M=Sigma_kron\original_data; % solving least square problem

for ij=1:n
    
    D{ij}=reshape(M(:,ij),d,d);
    
end

%% MPGLRAM function
[L,D,R,M]=MPGLRAM__(original_data, DATA, L,D,R,n,K_,d,ite);
toc


