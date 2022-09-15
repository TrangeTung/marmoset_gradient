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
