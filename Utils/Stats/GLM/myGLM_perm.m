function GLM=myGLM_perm(GLM)

%Number of predictors (including intercept)
p = length(GLM.contrast);

%Number of independent GLM's to fit
M = size(GLM.y,2);

%Number of observations
n = size(GLM.y,1);

% if isfield(GLM,'exchange')
%     %Set up exchange blocks
%     blks=unique(GLM.exchange); 
%     %Number of blocks
%     n_blks=length(blks);
%     %Number of observations per block
%     sz_blk=n/n_blks; 
%     blk_ind=zeros(n_blks,sz_blk);
%     for nRlz=1:length(blks)
%         blk_ind(nRlz,:)=find(blks(nRlz)==GLM.exchange);
%     end
% end


%Variance Inflation Factor:
% VIFs are also the diagonal elements of the inverse of the correlation matrix [1], 
% a convenient result that eliminates the need to set up the various regressions
% [1] Belsley, D. A., E. Kuh, and R. E. Welsch. Regression Diagnostics. Hoboken, NJ: John Wiley & Sons, 1980.
R0      = corr(GLM.X(2:end,2:end)); % correlation matrix
GLM.vif = diag(pinv(R0))';

%Determine nuisance predictors not in contrast
ind_nuisance = find(~GLM.contrast);
ind_interest = find(GLM.contrast);

if ~isempty(ind_nuisance)
    disp(['Removing nausances...'])
    %Regress out nuisance predictors and compute residual
    Betas       = zeros(length(ind_nuisance),M);
    resid_y     = zeros(n,M);
    
    Betas       = GLM.X(:,ind_nuisance) \ GLM.y;
    yHat        = GLM.X(:,ind_nuisance) * Betas;
    resid_y     = GLM.y - yHat; 
end

%%% Some basic info from the mode -- betas and residuals 
% Note that these estimates are Beta and Res of *deconfounded* model
GLM.Betas  = GLM.X(:,ind_interest) \ resid_y;
resid_yhat =(GLM.X(:,ind_interest) * GLM.Betas);
GLM.Res    = resid_y - resid_yhat;
GLM.yhat   = resid_yhat;

%%% Model Diagnostics; AIC and BIC
GLM.AIC   = 2*(p+1) +n*log(sum(GLM.Res.^2)/n);
GLM.BIC   = log(n)*(p+1) +n*log(sum(GLM.Res.^2)/n);

GLM.RawBetas  = GLM.X \ GLM.y;

%%% Testing

test_stat=zeros(GLM.perms+1,M);

disp(['Permutation will start now...'])
for nRlz=1:GLM.perms+1
    
    if ~mod(nRlz,1000); disp(['myGLM_perm:: realisation: ' num2str(nRlz)]); end;
    
    y_perm=zeros(n,M);
    % Permute
    
    if nRlz==1
        % Don't permute the first run
        y_perm=GLM.y; 
    else
        if isempty(ind_nuisance) 
            %Permute signal 
            y_perm  = GLM.y(randperm(n)',:);
        else
            resid_y = resid_y(randperm(n)',:);
        end
    end
    % Freedman & Lane
    % Add permuted residual back to nuisance signal, giving a realisation 
    % of the data under the null hypothesis 
    
    if ~isempty(ind_nuisance)
        y_perm=resid_y+[GLM.X(:,ind_nuisance)]*Betas;
    end
       
    Betas_perm = zeros(p,M);
    Betas_perm = GLM.X \ y_perm;
        
    %Compute statistic of interest
    if strcmp(GLM.test,'onesample')
        %Flip signs
        if nRlz==1
            %Don't permute first run
            % SA: Oh, this is THE test stat
            % How the math would work here?! Why does he mean over y_perm?!
            % Because we test on the MEAN?!
            test_stat(nRlz,:) = mean(y_perm); 
        else
            test_stat(nRlz,:) = mean(y_perm.*repmat(sign(rand(n,1)-0.5),1,M)); 
        end
        
    elseif strcmp(GLM.test,'ttest')
        resid             = zeros(n,M);
        mse               = zeros(n,M);
        
        % SA: Can you please write the math down somewhere?
        resid             = y_perm-GLM.X*Betas_perm;
        mse               = sum(resid.^2)/(n-p);
        se                = sqrt(mse* (GLM.contrast * pinv(GLM.X'*GLM.X) * GLM.contrast') );
        test_stat(nRlz,:) = (GLM.contrast * Betas_perm)./se;
        
    elseif strcmp(GLM.test,'ftest')
        sse               = zeros(1,M);
        ssr               = zeros(1,M);
        %Sum of squares due to error
        sse               = sum((y_perm-GLM.X*Betas_perm).^2);
        %Sum of square due to regression
        ssr               = sum((GLM.X*Betas_perm-repmat(mean(y_perm),n,1)).^2);
        if isempty(ind_nuisance)
            test_stat(nRlz,:) = (ssr/(p-1))./(sse/(n-p));
        else
            %Get reduced model
            %Column of ones will be added to the reduced model unless the
            %resulting matrix is rank deficient
            X_new           = [ones(n,1),GLM.X(:,ind_nuisance)];
            %+1 because a column of 1's will be added to the reduced model
            b_red           = zeros(length(ind_nuisance)+1,M);
            %Number of remaining variables
            v               = length(find(GLM.contrast))-1;
            [n,ncolx]       = size(X_new);
            [Q,R,perm]      = qr(X_new,0);
            rankx           = sum(abs(diag(R)) > abs(R(1))*max(n,ncolx)*eps(class(R)));
            
            if rankx < ncolx
                %Rank deficient, remove column of ones
                X_new   = GLM.X(:,ind_nuisance);
                b_red   = zeros(length(ind_nuisance),M);
                v       = length(find(GLM.contrast));
            end
            sse_red           = zeros(1,M);
            ssr_red           = zeros(1,M);
            
            b_red             = X_new\y_perm;
            sse_red           = sum((y_perm-X_new*b_red).^2);
            ssr_red           = sum((X_new*b_red-repmat(mean(y_perm),n,1)).^2);
            test_stat(nRlz,:) = ((ssr-ssr_red)/v)./(sse/(n-p));
        end
    end
end

%Added to v1.1.2
%Covers the case where the dependent variable is identically zero for all
%observations. The test statistic in this case in NaN. Therefore, force any
%NaN elements to zero. This case is typical of connectivity matrices
%populated using streamline counts, in which case some regional pairs are
%not interconnected by any streamlines for all subjects. 
test_stat(isnan(test_stat))=0; 

GLM.test_stat = test_stat;

%%% Get the pvalues

%K       = size(test_stat,1)-1;
%pvals   = zeros(1,M);
%for i=1:M
%        pvals(i,:) = mean(test_stat(1,:)<=test_stat(2:end,:));
%end
%pvals
%pvals = pvals/K;

pvals = mean(abs(test_stat(1,:))<=test_stat(2:end,:));
GLM.unadj_pvals = pvals; 

if strcmpi(GLM.mce,'fdr')

    [GLM.h,GLM.crit_p,GLM.adj_p] = myGLM_fdr(GLM);
    
elseif strcmpi(GLM.mce,'fwe')
    
    max_test_stat = max(test_stat(2:end,:),[],2);
    GLM.adj_p  = mean(abs(test_stat(1,:))<=max_test_stat); 
    GLM.crit_p = []; 
    GLM.h      = GLM.adj_p<GLM.alpha;
    
else
    error('Unrecongnised MCE method; choose either fwe or fdr')
end

%%%% END OF FUNCTION 
end

function [h,crit_p,adj_p] = myGLM_fdr(GLM)
% Outputs:
%   h       - A binary vector or matrix of the same size as the input "pvals."
%             If the ith element of h is 1, then the test that produced the 
%             ith p-value in pvals is significant (i.e., the null hypothesis
%             of the test is rejected).
%   crit_p  - All uncorrected p-values less than or equal to crit_p are 
%             significant (i.e., their null hypotheses are rejected).  If 
%             no p-values are significant, crit_p=0.
%   adj_p   - All adjusted p-values less than or equal to q are significant
%             (i.e., their null hypotheses are rejected). Note, adjusted 
%             p-values can be greater than 1.
%
%
% References:
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.

pvals = GLM.unadj_pvals; 
s     = size(pvals);
q     = GLM.alpha; 

%p-values are already a row vector
[p_sorted, sort_ids] = sort(pvals);
[~, unsort_ids]      = sort(sort_ids); %indexes to return p_sorted to pvals order
m                    = length(p_sorted); %number of tests

%BH procedure for independence or positive dependence
thresh = (1:m)*q/m;
wtd_p  = m*p_sorted./(1:m);

%compute adjusted p-values
adj_p = zeros(1,m)*NaN;
[wtd_p_sorted, wtd_p_sindex] = sort(wtd_p);
nextfill = 1;
for k = 1 : m
    if wtd_p_sindex(k)>=nextfill
        adj_p(nextfill:wtd_p_sindex(k)) = wtd_p_sorted(k);
        nextfill = wtd_p_sindex(k)+1;
        if nextfill>m
            break;
        end
    end
end

adj_p  = reshape(adj_p(unsort_ids),s);
rej    = p_sorted<=thresh;
max_id = find(rej,1,'last'); %find greatest significant pvalue

if isempty(max_id)
    crit_p = 0;
    h      = pvals*0;
else
    crit_p = p_sorted(max_id);
    h      = pvals<=crit_p;
end

%%%% END OF FUNCTION 
end
    