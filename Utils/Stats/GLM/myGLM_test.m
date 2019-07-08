clear

load(['/Users/sorooshafyouni/Home/GitClone/UKBUtils/Analysis/T1MRI/20K/AtAsCent/R_Confounds/UKB-Confounds-AtAsCent-VBM-20K-eid34077.mat'])
load(['/Users/sorooshafyouni/Home/GitClone/UKBUtils/Analysis/T1MRI/20K/AtAsCent/R_GMVols/UKB-HvrdOxGM-AtAsCent-VBM-20K-eid34077.mat'])
% 
% %{'onesample','ttest','ftest'}
% 
% N = size(RegMat_Ind,1);
% 
%load(['/Volumes/Extra/UKB/GEO/VBM-20K-eid34077/UKB-GeoDMat-inverse-K20-AtAsCent-Bin0-VBM-20K-eid34077.mat'],'Lat20K','Lon20K')
%X = [RegMat_IndT1Cont];
%y = HOAGMV;
%P = size(X,2)+1; 
%N = size(X,1); 
%MDL = ['Intrcpt', MDL_IndT1Cont];
%mdl = fitglm(X,HOAGMV(:,1),'VarNames',[MDL_IndT1Cont,'HOGM']);
%ii = 20;
%C_tmp = zeros(1,P);
%C_tmp(ii) = 1;
%y = log(y); 
%X(:,ii-1) = X(:,ii-1).^3;

Interest = 1;
X = [RegMat_IndT1 RegMat_Cont(:,Interest)];
MDL_Cont(Interest)

P = size(X,2)+1; 
C = zeros(1,P);
C(end) = 1;

GLM_tmp = myGLM_MakeGLM(X,HOAGMV,'ttest',C,'fwe',5000);  
GG_mce  = myGLM_perm(GLM_tmp);

if sum(GG_mce.h)
    UKBparse_Labels(GG_mce.h)
end

%[GG_mce.Betas(1),GG_mce.test_stat(1,1),GG_mce.unadj_pvals(1)]

%GLM.se = sqrt(mse* (GLM_tmp.contrast * pinv(GLM_tmp.X'*GLM_tmp.X) * GLM_tmp.contrast') );

ROI = 1; 
% WhichOnes = [1:20];
% X = RegMat_IndT1Cont;
% X = X-mean(X); 
% X = X(:,WhichOnes);

Y = HOAGMV;
Y = Y - mean(Y); 
Y = Y(:,ROI);
% ,'VarName',[MDL_IndT1Cont,'HOGMA']
mdl_tmp = fitglm(X,Y); mdl_tmp
