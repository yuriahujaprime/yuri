# LASSO Gauss Copula implementation script
# Infers regression coefficients maximizing the LASSO Gauss Copula likelihood function
# Yuri Ahuja
# Last modified 8/21/2019


import xalglib
import numpy as np
import scipy as sc
import scipy.stats as st
import math
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score


# User sets input datasets and array of lambdas to validate
lambs = 10**np.array([-2.5,-2.0,-1.5,-1.0,-0.5])
X = np.loadtxt('Xtrain.csv',skiprows=1,delimiter=',')			# Training feature set
y = np.loadtxt('Ytrain.csv',skiprows=1,delimiter=',')			# Training label set
Xtest = np.loadtxt('Xtest.csv',skiprows=1,delimiter=',')		# Test feature set
ytest = np.loadtxt('Ytest.csv',skiprows=1,delimiter=',')		# Test label set
beta_actual = np.loadtxt('betas.csv',skiprows=1,delimiter=',')	# True beta (assuming data is simulated)
n,p = np.shape(X)												# n = number of individuals; p = number of features

# Computes data correlation matrix (with diagonal offset to avoid singularity) to be used for beta covariance matrix
# CorrMat = np.corrcoef(X.transpose()) + 1e-6*np.identity(p)	# Uncomment to compute matrix in python
CorrMat = np.loadtxt('Sigma.csv',skiprows=1,delimiter=',')		# Computed corr(X) in R and added 1e-6 to diagonal
CorrMatInv = np.linalg.inv(CorrMat)								# Inverse of correlation matrix

# Makes X design matrix with intercept column of 1s
X = np.column_stack((np.ones(n),X))

## Set Sig to empiric covariance matrix in X vs. identity matrix
# Sig = np.identity(p)
# for i in np.arange(0,p,2):
# 	Sig[i,i+1] = 0.9
# 	Sig[i+1,i] = 0.9


# Maps p-dimensional input vector to p-dimensional output vector
# See Eric's pdf on Lasso Gauss Copula for defintion of q()
def q(x):
	return np.array([st.norm.ppf(st.laplace.cdf(x[i],scale=1/(n*lamb)), scale=1/(n*lamb)) for i in range(len(x))])


# Maps p-dimensionl input vector to p-dimensional output vector corresponding to diag(Jacobian(q))
# Note that Jacobian(q) here is a pxp diagnoal matrix; this function JUST returns the diagonal
def Jq(x):
	G = st.laplace.cdf(x,scale=1/(n*lamb))
	Gprime = st.laplace.pdf(x,scale=1/(n*lamb))
	return np.array([Gprime[i]/st.norm.pdf(st.norm.ppf(G[i], scale=1/(n*lamb)), scale=1/(n*lamb)) for i in range(len(G))])


# Defines objective function (fi[0]), equality constraint (fi[1] = 0), inequality constraints (fi[2:p+1] < 0)
# and Jacobian matrix (J[i][j] = [dfi/dxj]) (see Eric's pdf for derivations of these terms)
# See http://www.alglib.net/translator/man/manual.cpython.html#example_minnlc_d_mixed for nlcfunc2_jac specification
# Inputs:
# betapm = (2p+1)-dimensional vector consisting of [beta0,beta+,beta-] (see Eric's pdf for definition of beta+ and beta-)
# fi, jac, and param do not have to be predefined; XALGLIB will set these
def nlcfunc2_jac(betapm,fi,jac,param):
	#
	# functions:
	# f0(beta): y^TX(beta+ - beta-) - Sum[log(exp((beta+ - beta-)^Txi)+1)] - 0.5q(beta+ - beta-)^T(Sigma^-1 - Sigmadiag^-1)q(beta+ - beta-) - lambda * Sum[beta+ + beta-]
	# f1(beta): (beta+)^T(beta-) = 0
	# f2(beta): beta+,beta- >= 0
	#
	beta0 = betapm[0]						# Intercept term beta0
	betaP = np.array(betapm[1:(p+1)])		# beta+
	betaM = np.array(betapm[(p+1):])		# beta-
	betaR = np.append(beta0,betaP-betaM)	# beta = [beta0, beta+ - beta-]
	YtX = np.dot(y,X)
	#
	# Set fi
	fi[0] = float(-1.0*np.dot(YtX,betaR) + 1.0*np.sum(np.log(1.0+np.exp(np.dot(X,betaR)))) + 0.5 * np.dot(np.transpose(q(betaR[1:])), np.dot(SigDiff, q(betaR[1:]))) + n*lamb * np.sum(betaP+betaM))
	fi[1] = float(np.dot(betaP,betaM))		# Sigdiff = (Sigma^-1 = Sigmadiff^-1) defined below; see Eric's pdf for details
	for i in range(2*p):
		fi[i+2] = -1.0*float(betapm[i+1])
	#
	# Set jac
	exponent = np.exp(np.dot(X,betaR))
	term1 = 1.0*YtX												# Log likelihood component of gradient part 1
	term2 = 1.0*np.dot(np.transpose(X), exponent/(exponent+1))	# Log likelihood component of gradient part 2
	term3 = Jq(betaR[1:]) * np.dot(SigDiff, q(betaR[1:]))		# LASSO Gauss Copula component of gradient (note that Jq is the Jacobian DIAGONAL)
	term4 = n * lamb * np.ones(p)								# Marginal LASSO component of gradient
	jac[0] = [float(-1.0*term1[0]+term2[0])] + list(map(float,-1.0*term1[1:]+term2[1:]+term3+term4)) + list(map(float,term1[1:]-term2[1:]-term3+term4))
	jac[1] = [0.0] + list(map(float,betaM)) + list(map(float,betaP))
	for i in range(2*p):
		jac[i+2] = [0.0 for j in range(2*p+1)]
		jac[i+2][i+1] = -1.0


# Runs copula on set of lambdas and datasets specified at top of page
for lamb in lambs:
	print("Lambda = " + str(lamb))
	#
	SigDiff = (n*lamb)**2 * CorrMatInv		# SigDiff = Sigma^-1 - Sigmadiag^-1
	for k in range(p):
		SigDiff[k,k] = 0.0
	#
	# Fits LASSO (note normalization parameter must be scaled by n) to use as initial beta
	lasso = LogisticRegression(penalty='l1', C=1/(n*lamb), solver='liblinear')
	lasso.fit(X[:,1:],y)
	betaInit = list(lasso.intercept_)+[max(c,0) for c in lasso.coef_[0]]+[abs(min(c,0)) for c in lasso.coef_[0]]
	#
	# Sets parameter scales (XALGLIB requirement) to abs(maximum parameter from LASSO fit)
	s = list(np.max(np.abs(betaInit)) * np.ones(len(betaInit)))
	epsx = 1e-6		# Tolerance parameter (descent halts when maximum parameter step falls below tolerance)
	maxits = 0		# Max number of iterations (0 = unlimited)
	rho = 1000.0	# Some kind of tolerance thing only used in Augmented Lagrangian (not used here)
	outerits = 5	# Another parameter only used in Augmented Lagrangian (not used here)
	#
	state = xalglib.minnlccreate(len(betaInit), betaInit)	# Creates minnlc (min-nonlinearconstrained) object
	xalglib.minnlcsetcond(state,epsx,maxits)		# Sets stopping conditions (epsx and maxits)
	xalglib.minnlcsetscale(state,s)					# Sets scales of coefficients (s)
	xalglib.minnlcsetstpmax(state,10.0)				# Sets maximum step length (after scaling by s)
	xalglib.minnlcsetalgoaul(state,rho,outerits)	# Activates Augmented Lagrangian solver (AUL)
	xalglib.minnlcsetalgoslp(state)					# Activates Successive Linear Programming (SLP), comment out to stick with AUL
	xalglib.minnlcsetnlc(state,1,2*p)				# Sets number of equality (1) and inequality (2p) constraints
	xalglib.minnlcoptguardsmoothness(state)			# Sets error flag to 1 if detects dicontinuity/non-smoothness
	xalglib.minnlcoptguardgradient(state,0.001)		# Sets error flag to 1 if occasional numerical differentiation of fi significantly differs from user-supplied jac
	xalglib.minnlcoptimize_j(state,nlcfunc2_jac)	# Performs optimization according to functions defined in nlcfunc2_jac above
	x1, rep = xalglib.minnlcresults(state)			# Saves (beta+/beta-) estimator in x1
	#
	# Automatically sets any negative (beta+/beta-) estimators to 0 (should be all nonnegative)
	for i in range(1,len(x1)):
		x1[i] = max(x1[i],0.0)
	#
	beta0 = x1[0]							# Intercept term beta0
	betaP = np.array(x1[1:(p+1)])			# beta+
	betaM = np.array(x1[(p+1):])			# beta-
	betaR = np.append(beta0,betaP-betaM)	# beta = [beta0, beta+ - beta-]
	np.savetxt("pascal_copula_beta_lambda-"+str(lamb)+".csv",betaR,delimiter=',')	# Saves beta estimator
	#
	# Computes cosine similarities with true beta (assuming simulated beta is provided) and predictive AUCs
	print("Cosine similarity between LASSO and true beta (lambda=" + str(lamb) + "): " + str(np.dot(lasso.coef_[0],beta_actual)/(np.linalg.norm(lasso.coef_[0])*np.linalg.norm(beta_actual))))
	print("Cosine similarity between Copula and true beta (lambda=" + str(lamb) + "): " + str(np.dot(betaR[1:],beta_actual)/(np.linalg.norm(betaR[1:])*np.linalg.norm(beta_actual))))
	print("LASSO AUC (lambda=" + str(lamb) + "): " + str(roc_auc_score(ytest,np.dot(Xtest,lasso.coef_[0]))))
	print("Copula AUC (lambda=" + str(lamb) + "): " + str(roc_auc_score(ytest,np.dot(Xtest,betaR[1:]))))
	#
	# Quality control flags (all 0s indicates no issues; 1s indicate issues detected during optimization)
	ogrep = xalglib.minnlcoptguardresults(state)
	print("Quality control (all 0s indicates no issues with solver):")
	print(ogrep.badgradsuspected, ogrep.nonc0suspected, ogrep.nonc1suspected)


# Reads in copula beta estimates saved above and retrospectively computes cosine similarity/predictive AUC
# No need to uncomment unless you need to recompute accuracy metrics 
for lamb in lambs:
	print("Lambda = " + str(lamb))
	lasso = LogisticRegression(penalty='l1', C=1/(n*lamb), fit_intercept=True, solver='liblinear')
	lasso.fit(X[:,1:],y)
	betaR = np.loadtxt("beta_learned_lambda-"+str(lamb)+".csv",delimiter=',')
	print("Cosine similarity between LASSO and true beta (lambda=" + str(lamb) + "): " + str(np.dot(lasso.coef_[0],beta_actual)/(np.linalg.norm(lasso.coef_[0])*np.linalg.norm(beta_actual))))
	print("Cosine similarity between Copula and true beta (lambda=" + str(lamb) + "): " + str(np.dot(betaR[1:],beta_actual)/(np.linalg.norm(betaR[1:])*np.linalg.norm(beta_actual))))
	print("LASSO AUC (lambda=" + str(lamb) + "): " + str(roc_auc_score(ytest,np.dot(Xtest,lasso.coef_[0]))))
	print("Copula AUC (lambda=" + str(lamb) + "): " + str(roc_auc_score(ytest,np.dot(Xtest,betaR[1:]))))

