function [CC,R2] = MY_cpuArray_GLM_shuffle(V,DesignMatrix,lamda,shuf_column)
% Vec = gpuArray(V);
% DM =  gpuArray(DesignMatrix);
% fa =  gpuArray(lamda);
Vec = (V);
DM =  (DesignMatrix);
fa =  (lamda);
% for i=1:1000
S.v = Vec;
S.dm = DM;
S.fa = fa;
S.sc = shuf_column;
% S{i} = {Vec;DM;fa;shuf_column};
S.mI = 1000;

clear fa DM Vec

[CCgpu,R2gpu] = arrayfun(@MY_GLM_CC_R2_estimate,S,'UniformOutput',false);
clear S

C = cell2mat(CCgpu)/1000;
R =  cell2mat(R2gpu)/1000;
clear CCgpu R2gpu 
% CC = gather(C);
% R2 = gather(R);
CC = (C);
R2 = (R);
clear C R

end


function [CCgpu,R2gpu] = MY_GLM_CC_R2_estimate(S)

% Vec = S{1};
% DM = S{2};
% fa =S{3};
% shuf_column = S{4};



Vec = S.v;
DM = S.dm;
fa =S.fa;
shuf_column = S.sc;
maxIter = S.mI;

R2gpu = (zeros(size(Vec,1),1));
CCgpu = (zeros(size(Vec,1),1));
for iter = 1:maxIter
    
for shf = shuf_column
    X = randperm(size(DM,1));
    DM(:,shf) = DM(X,shf);
end

A = DM'*DM+fa*eye(size(DM,2));
b = DM'*Vec';
beta = A\b; clear A b
VecP = beta'*DM';

CovXY = sum((Vec-mean(Vec,2)).*(VecP-mean(VecP,2)),2);
sigmaX = sqrt(sum((Vec-mean(Vec,2)).^2,2));
sigmaY = sqrt(sum((VecP-mean(VecP,2)).^2,2));
C = CovXY./sigmaX./sigmaY;
CCgpu = CCgpu+C(:);%
clear sigma* Cov* C
%
mu = mean(VecP,2).*mean(Vec,2)-mean(Vec.*VecP,2);
md = mean(VecP,2).^2 - mean(VecP.^2,2);
b = mean(Vec,2) - mu./md.*mean(VecP,2);
Vecf = mu./md.*VecP+b;
SSR = sum((Vecf-mean(Vec,2)).^2,2);
SST = sum((Vec-mean(Vec,2)).^2,2);
R = SSR./SST;
R2gpu = R2gpu+R(:);%
clear VecP Vecf SS* m* beta b R
end

end
