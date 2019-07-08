function [h,crit_p,adj_p] = myGLM_fdr(GLM)

pvals = GLM.pvals; s=size(pvals);
q     = GLM.alpha; 

%p-values are already a row vector
[p_sorted, sort_ids] = sort(pvals);
[~, unsort_ids]      = sort(sort_ids); %indexes to return p_sorted to pvals order
m=length(p_sorted); %number of tests

%BH procedure for independence or positive dependence
thresh = (1:m)*q/m;
wtd_p  = m*p_sorted./(1:m);

%compute adjusted p-values
adj_p = zeros(1,m)*NaN;
[wtd_p_sorted, wtd_p_sindex] = sort( wtd_p );
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


end
