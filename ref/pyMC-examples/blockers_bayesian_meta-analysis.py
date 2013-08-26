# -*- coding: utf-8 -*-
# <nbformat>3</nbformat>

# <markdowncell>

# Blocker: random effects meta-analysis of clinical trials
# ========================================================
# 
# Ported to PyMC by Abraham Flaxman, 4/18/2012, from http://www.openbugs.info/Examples/Blockers.html
# 
# Carlin (1992) considers a Bayesian approach to meta-analysis, and includes the following examples of 22 trials of beta-blockers to prevent mortality after myocardial infarction.
# 
# In a random effects meta-analysis we assume the true effect (on a log-odds scale) $d_i$ in a trial $i$
# is drawn from some population distribution. Let $r^C_i$ denote number of events in the control group in trial $i$,
# and $r^T_i$ denote events under active treatment in trial $i$. Our model is:
# 
# $$
# r^C_i \sim \text{Binomial}\left(p^C_i, n^C_i\right),
# $$
# $$
# r^T_i \sim \text{Binomial}\left(p^T_i, n^T_i\right),
# $$
# $$
# \text{logit}\left(p^C_i\right) = \mu_i,
# $$
# $$
# \text{logit}\left(p^T_i\right) = \mu_i + \delta_i,
# $$
# $$
# \delta_i \sim \text{Normal}(d, t).
# $$
# 
# "Noninformative" priors are given for $\mu_i$, $t$, and $d$. We want to make inferences about the population effect $d$,
# and the predictive distribution for the effect $\delta_{\text{new}}$ in a new trial. Empirical Bayes methods estimate
# $d$ and $t$ by maximum likelihood and use these estimates to form the predictive distribution $p( \delta_\text{new} | \hat{d}, \hat{t})$. Full Bayes allows for the uncertainty concerning $d$ and $t$ . 

# <codecell>

import pylab as pl
import pymc as mc

# <codecell>

### data
r_t_obs = [3, 7, 5, 102, 28, 4, 98, 60, 25, 138, 64, 45, 9, 57, 25, 33, 28, 8, 6, 32, 27, 22]
n_t_obs = [38, 114, 69, 1533, 355, 59, 945, 632, 278,1916, 873, 263, 291, 858, 154, 207, 251, 151, 174, 209, 391, 680]
r_c_obs = [3, 14, 11, 127, 27, 6, 152, 48, 37, 188, 52, 47, 16, 45, 31, 38, 12, 6, 3, 40, 43, 39]
n_c_obs = [39, 116, 93, 1520, 365, 52, 939, 471, 282, 1921, 583, 266, 293, 883, 147, 213, 122, 154, 134, 218, 364, 674]
N = len(n_c_obs)

# <markdowncell>

# Here is the BUGS model:
# 
#     model
#     {
#        for( i in 1 : Num ) {
#           rc[i] ~ dbin(pc[i], nc[i])
#           rt[i] ~ dbin(pt[i], nt[i])
#           logit(pc[i]) <- mu[i]
#           logit(pt[i]) <- mu[i] + delta[i]
#           mu[i] ~ dnorm(0.0,1.0E-5)
#           delta[i] ~ dnorm(d, tau)
#        }
#        d ~ dnorm(0.0,1.0E-6)
#        tau ~ dgamma(0.001,0.001)
#        delta.new ~ dnorm(d, tau)
#        sigma <- 1 / sqrt(tau)
#     }

# <codecell>

### hyperpriors
d = mc.Normal('d', 0., 1.e-6, value=0.)
tau = mc.Gamma('tau', 1.e-3, 1.e-3, value=1.)

sigma = mc.Lambda('sigma', lambda tau=tau: tau**-.5)
delta_new = mc.Normal('delta_new', d, tau, value=0.)


### priors
mu = [mc.Normal('mu_%d'%i, 0., 1.e-5, value=0.) for i in range(N)]
delta = [mc.Normal('delta_%d'%i, d, tau, value=0.) for i in range(N)]

p_c = mc.Lambda('p_c', lambda mu=mu: mc.invlogit(mu))
p_t = mc.Lambda('p_t', lambda mu=mu, delta=delta: mc.invlogit(array(mu)+delta))


### likelihood
r_c = mc.Binomial('r_c', n_c_obs, p_c, value=r_c_obs, observed=True)
r_t = mc.Binomial('r_t', n_t_obs, p_t, value=r_t_obs, observed=True)

# <markdowncell>

# BUGS uses Gibbs steps automatically, so it only takes 10000 steps of MCMC after a 1000 step burn in for this model in their example.
# 
# PyMC only uses Gibbs steps if you set them up yourself, and it uses Metropolis steps by default.  So 10000 steps
# go by more quickly, but the chain takes longer to converge to the stationary distribution.

# <codecell>

m = mc.MCMC([d, tau, sigma, delta_new, mu, delta, p_c, p_t, r_c, r_t])

# not long enough to mix
%time m.sample(11000, 1000, progress_bar=False)

# <codecell>

mc.Matplot.plot(sigma)
mc.Matplot.plot(d)
mc.Matplot.plot(delta_new)

# <markdowncell>

# Using some fancy step methods, and thinning the chain a little gives more reliable results

# <codecell>

%time mc.MAP([d, tau, sigma, delta_new, mu, delta, p_c, p_t, r_c, r_t]).fit(method='fmin_powell')
m = mc.MCMC([d, tau, sigma, delta_new, mu, delta, p_c, p_t, r_c, r_t])

%time m.sample(110000, 10000, 10, progress_bar=False)

# <codecell>

mc.Matplot.plot(sigma)
mc.Matplot.plot(d)
mc.Matplot.plot(delta_new)

# <markdowncell>

# Results
# -------
# 
# The OpenBUGS results table shows the following:
# 
#                  mean     sd
#     d           -0.2492   0.06422
#     delta_new   -0.2499   0.1509
#     sigma        0.1189   0.07
# 
# To be compared with the following:

# <codecell>

for node in [d, delta_new, sigma]:
    stats = node.stats()
    print '%10s\t%1.3f \t%.3f' % (node.__name__, stats['mean'], stats['standard deviation'])

# <codecell>

# use the following to generate plots of all stochastics
# (good for checking convergence in detail)
# mc.Matplot.plot(m)

# <markdowncell>

# Extensions
# ----------
# 
# In some circumstances it might be reasonable to assume that the population distribution has heavier tails, for example a $t$-distribution with low degrees of freedom. This is easily accomplished in BUGS by using the <code>dt</code> distribution function instead of <code>dnorm</code> for $d$ and $d_\text{new}$.
# 
# In PyMC, the relevant stochastic is <code>mc.NoncentralT</code>

# <codecell>

### hyperpriors
d = mc.Normal('d', 0., 1.e-6, value=0.)
tau = mc.Gamma('tau', 1.e-3, 1.e-3, value=1.)

sigma = mc.Lambda('sigma', lambda tau=tau: tau**-.5)
delta_new = mc.NoncentralT('delta_new', d, tau, 4., value=0.)


### priors
mu = [mc.Normal('mu_%d'%i, 0., 1.e-5, value=0.) for i in range(N)]
delta = [mc.NoncentralT('delta_%d'%i, d, tau, 4., value=0.) for i in range(N)]

p_c = mc.Lambda('p_c', lambda mu=mu: mc.invlogit(mu))
p_t = mc.Lambda('p_t', lambda mu=mu, delta=delta: mc.invlogit(array(mu)+delta))


### likelihood
r_c = mc.Binomial('r_c', n_c_obs, p_c, value=r_c_obs, observed=True)
r_t = mc.Binomial('r_t', n_t_obs, p_t, value=r_t_obs, observed=True)

# <markdowncell>

# Our estimates are lower and with tighter precision - in fact similar to the values obtained by Carlin for the empirical Bayes estimator. The discrepancy appears to be due to Carlin's use of a uniform prior for s 2 in his analysis, which will lead to increased posterior mean and standard deviation for d, as compared to our (approximate) use of p( s 2 )  ~ 1 /  s 2 (see his Figure 1).
# 
# It is easy to use a uniform prior on $s^2$ to compare results.

# <markdowncell>

# And how can we use PyMC to do empirical bayes?  That requires going to Carlin's paper to see what is done there precisely.

# <codecell>


