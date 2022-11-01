function [CC,R2,CC_node,R2_node,VecP] = MY_gpuArray_GLM_shuffle(V,DesignMatrix,lamda,shuf_column)
Vec = gpuArray(V);
DM =  gpuArray(DesignMatrix);
fa =  gpuArray(lamda);

% for i=1:1000
S.v = Vec;
S.dm = DM;
S.fa = fa;
S.sc = shuf_column;
% S{i} = {Vec;DM;fa;shuf_column};
S.mI = 1000;

clear fa DM Vec

[CCgpu,R2gpu,CCgpu_node,R2gpu_node,VecPgpu] = arrayfun(@MY_GLM_CC_R2_estimate,S,'UniformOutput',false);

C = cell2mat(CCgpu)/S.mI;
R =  cell2mat(R2gpu)/S.mI;
clear CCgpu R2gpu 
C_node = cell2mat(CCgpu_node)/S.mI;
R_node =  cell2mat(R2gpu_node)/S.mI;
clear CCgpu_node R2gpu_node 
VecPx = cell2mat(VecPgpu)/S.mI;
clear VecPgpu 



CC = gather(C);
R2 = gather(R);
CC_node = gather(C_node);
R2_node = gather(R_node);
VecP = gather(VecPx);
clear C R VecPx


end


function [CCgpu,R2gpu,CCgpu_node,R2gpu_node,VecPgpu] = MY_GLM_CC_R2_estimate(S)

% Vec = S{1};
% DM = S{2};
% fa =S{3};
% shuf_column = S{4};



Vec = S.v;
DM = S.dm;
fa =S.fa;
shuf_column = S.sc;
maxIter = S.mI;
[Nx,Ny] = size(Vec);
R2gpu = gpuArray(zeros(size(Vec,1),1));
CCgpu = gpuArray(zeros(size(Vec,1),1));
R2gpu_node = gpuArray(zeros(size(Vec,1),sqrt(Ny)));
CCgpu_node = gpuArray(zeros(size(Vec,1),sqrt(Ny)));
VecPgpu = gpuArray(zeros(size(Vec,2),size(Vec,1)));

for iter = 1:maxIter
    
for shf = shuf_column
    X = randperm(size(DM,1));
    DM(:,shf) = DM(X,shf);
end

A = DM'*DM+fa*eye(size(DM,2));
b = DM'*Vec';
beta = A\b; clear A b
VecP = (beta'*DM');

[C,R]=MY_estimate_CC_R2_perRow(Vec,VecP);
[Nx,Ny] = size(Vec);
A0=squeeze(reshape(Vec,Nx,sqrt(Ny),sqrt(Ny)));
B0=squeeze(reshape(VecP,Nx,sqrt(Ny),sqrt(Ny)));
clear C_node R_node
for nl=1:sqrt(Ny)
[C_node(nl,:),R_node(nl,:)]=MY_estimate_CC_R2_perRow(A0(:,:,nl),B0(:,:,nl));
end

CCgpu = CCgpu+C(:);%
R2gpu = R2gpu+R(:);%
CCgpu_node = CCgpu_node+C_node';%
R2gpu_node = R2gpu_node+R_node';%
VecPgpu = VecPgpu+VecP';

clear sigma* Cov* C
clear VecP Vecf SS* m* beta b R
end

end

function [CC,R2]=MY_estimate_CC_R2_perRow(Vec,VecP)
% Vec: raw data; point*time
% VecP: prediected data; point*time
%
CovXY = sum((Vec-mean(Vec,2)).*(VecP-mean(VecP,2)),2);
sigmaX = sqrt(sum((Vec-mean(Vec,2)).^2,2));
sigmaY = sqrt(sum((VecP-mean(VecP,2)).^2,2));
CC = CovXY./sigmaX./sigmaY;
%
mu = mean(VecP,2).*mean(Vec,2)-mean(Vec.*VecP,2);
md = mean(VecP,2).^2 - mean(VecP.^2,2);
b = mean(Vec,2) - mu./md.*mean(VecP,2);
Vecf = mu./md.*VecP+b;
SSR = sum((Vecf-mean(Vec,2)).^2,2);
SST = sum((Vec-mean(Vec,2)).^2,2);
R2 = SSR./SST;

end
