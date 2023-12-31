# check if all following packages are installed, if not then install
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table,Rcpp,RcppArmadillo,optparse) 

library(Rcpp)
library(RcppArmadillo)
library(data.table)
library(optparse)

cat('\n### All required packages are successfully loaded.\n')

option_list = list(
  make_option(c("-g", "--gene"), type="character", default=NULL, 
              help="gene expression file name", metavar="character"),
  make_option(c("-p", "--pheno"), type="character", default=NULL, 
              help="phenotype file name", metavar="character"),
  make_option(c("-w", "--w-step"), type="numeric", default=6000, 
              help="number of warm-up steps in MCMC", metavar="numeric"),
  make_option(c("-s", "--s-step"), type="numeric", default=20000, 
              help="number of sampling steps in MCMC", metavar="numeric"),
  make_option(c("-o", "--out"), type="character", default="FABIO_PIP.txt", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$gene)){
  print_help(opt_parser)
  stop("### At least one argument must be supplied (input file)", call.=FALSE)
}

source('FABIO_utility.R')
sourceCpp('FABIO_utility.cpp')

X = fread(opt$gene)
X = as.data.frame(X)
gene_names = as.vector(X[,1])
X = scale(t(X[,-1]))
y = fread(opt$pheno)
y = y$pheno

w_step = opt$`w-step`
s_step = opt$`s-step`
cat('\n### Input data are successfully loaded.')

start_time = Sys.time()
cat('\n### Analysis starts at',format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat('\n')

## pre-define basic parameters 
# number of analyzed individuals
ni_test = nrow(X)
# number of analyzed genes
ng_test = ncol(X)
pve_null = 0
n_accept = 0
# -hmax, -hmin, h for PVE
h_min = 0
h_max = 1
# -hscale
h_scale = min(1, 10 / sqrt(ni_test))
# -pmin, -pmax, log(pi)
logp_min = -log(ng_test)
logp_max = 0
# -pscale
logp_scale = min(1, 5 / sqrt(ni_test))
# -gmin, -gmax, range of n_gamma, n_gamma: number of gamma
g_min = 0
g_max = 300
# -rpace, record one state in every [num] steps (default 10)
r_pace = 10
# -wpace, write values down in every [num] recorded steps (default 1000)
w_pace = 1000
# -mh, number of M-H steps per iteration
n_mh = 10
# -gmean, the given mean of expected number of genes in final model
geo_mean = 2000
trace_G = 0
randseed = -1

if (ng_test < 300) {g_max = ng_test}
# pre-allocate
Xb_new = rep(0, ni_test)
Xb_old = rep(0, ni_test)
z_hat = rep(0, ni_test)
z = rep(0, ni_test)

Xgamma_old = matrix(0, ni_test, g_max)
XtX_old = matrix(0, g_max, g_max)
Xtz_old = rep(0, g_max)
beta_old = rep(0, g_max)

Xgamma_new = matrix(0, ni_test, g_max)
XtX_new = matrix(0, g_max, g_max)
Xtz_new = rep(0, g_max)
beta_new = rep(0, g_max)

z = y

# center z
mean_z = mean(z)
z = z - mean_z
ztz = t(z) %*% z

pheno_mean = 0

# calculate and rank the log likelihood of linear model for each gene
pos_loglr = cbind(1:ncol(X), rep(0, ncol(X)))
for (i in 1:ncol(X)) {
  X_col = X[,i]
  xtx = t(X_col) %*% X_col
  xtz = t(X_col) %*% z
  
  log_lr = 0.5 * length(z) * (log(ztz) - log(ztz - xtz * xtz / xtx))
  pos_loglr[i,2] = log_lr
}
pos_loglr = pos_loglr[order(pos_loglr[,2], decreasing = T),]
# map rank to gene index
rank = 1:ncol(X)
pos = pos_loglr[,1]
mapRank2pos = cbind(rank, pos)
# create a pos vector for rcpp function
pos_vec = mapRank2pos[,2] - 1

if (randseed == -1) {
  c_time = unclass(as.POSIXlt(Sys.time()))
  randseed = round(c_time$hour * 3600 + c_time$min * 60 + c_time$sec)
}
set.seed(randseed)

# Calculate proposal distribution for gamma (unnormalized)
p_gamma = CalcPgamma(geo_mean, ng_test)

# Initial parameters.
result_IM = InitialMCMC(X, z, pos_loglr)
cHyp_old = result_IM$cHyp
rank_old = result_IM$rank

Xgamma_old[,1:length(rank_old)] = X[,mapRank2pos[rank_old,2]]
CalcXtX(Xgamma_old, z, length(rank_old), XtX_old, Xtz_old)

logPost_old = CalcPosterior(cHyp_old,ng_test,ztz,Xgamma_old,length(rank_old),XtX_old,Xtz_old,Xb_old,beta_old)
# Calculate centered z_hat, and pve.
CalcCC_PVEnZ(Xb_old, ni_test, cHyp_old, z_hat)

# start MCMC
total_step = w_step + s_step
Result_hyp = matrix(0, w_pace, 5)
Result_gamma = matrix(0, w_pace, g_max)
beta_g = matrix(0, ng_test, 2)

output = mcmc_iter(total_step, w_step, r_pace, w_pace, n_mh, ng_test, ni_test,
                   h_max, h_min, h_scale,
                   logp_max, logp_min, logp_scale,
                   y, z_hat, z, rank_old, beta_old, beta_new,
                   Xtz_old, Xtz_new, Xb_old, Xb_new, p_gamma,
                   pos_vec, cHyp_old, X, Xgamma_old, Xgamma_new, XtX_old, XtX_new,
                   Result_hyp, Result_gamma, beta_g)

end_time = Sys.time()
cat('\n### Analysis ends at',format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat(c("\n### Total time used for the analysis:",paste0(round(as.numeric(difftime(end_time,start_time,units=c("mins"))),2)," mins\n")))

w = output$w
beta_g = output$beta_g
beta = rep(NA, ng_test)
gamma = rep(NA, ng_test)
for (i in 1:ng_test) {
  if (beta_g[i,2] != 0) {
    beta[i] = beta_g[i,1] / beta_g[i,2]
    gamma[i] = beta_g[i,2] / w
  } else {
    beta[i] = 0
    gamma[i] = 0
  }
}
param = cbind(beta,gamma)
param = as.data.frame(cbind(gene_names,param))
param = param[,c(1,3)]
colnames(param) = c("Gene","PIP")
param = param[order(param$PIP,decreasing = T),]
write.table(param,opt$out,row.names = F,col.names = T,quote = F)
cat("\n### Data analysis is done, and results are saved.")
