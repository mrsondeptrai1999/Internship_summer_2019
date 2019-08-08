#Testing my BIC and AIC functions
my_rss = sum((y - beta %*% x)^2)
my_BIC(df4,my_rss) #510346.3
my_AIC(df4,my_rss) #510154.4
BIC(lm.fit2) #532201.3
AIC(lm.fit2) #532009.4
# Explanation for the slight difference is that we assume the data is Gaussian distribution and 
# calculte log-likelihood accordingly,
# While the package calculate log-likelihood using logLik, and I assume it does not use Gaussian as
# straightforward as I do