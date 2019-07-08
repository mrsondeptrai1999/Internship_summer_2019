function GLM = myGLM_MakeGLM(X,y,testmethod,C,mcemethod,Nperm)

% demean the regressors and y:
X = X-mean(X); 
y = y-mean(y); 

% also add a column of ones and check whether the size of the constrats
% matches the new regressor matrix. 
O = size(X,1);
X = [ones(O,1) X];
P = size(X,2);
N = size(y,2); 

if numel(C)~=P
    error('length of contrast doesnt match the RegMat -- have you considered a intercept?')
end

if ~exist('Nperm','var'); Nperm = 2000; end; 

GLM.P        = P; 
GLM.O        = O;
GLM.N        = N;
GLM.perms    = Nperm; 
GLM.X        = X;
GLM.y        = y;
GLM.contrast = C;
GLM.test     = testmethod;
GLM.alpha    = 0.05;
GLM.mce      = mcemethod;

disp(['Number of predictors: ' num2str(P) ' (including intercept)'])
disp(['Number of observations/subjects: ' num2str(O)])
disp(['Number of ROIs: ' num2str(N)])


