function vif = myGLM_VIF(X)

R0      = corr(X); % correlation matrix
vif = diag(pinv(R0))';