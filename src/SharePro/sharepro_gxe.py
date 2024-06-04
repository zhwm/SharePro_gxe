import argparse
import logging
import pandas as pd
import numpy as np
from scipy.special import softmax, expit
from scipy.stats import chi2

np.set_printoptions(precision=4, linewidth=200)

def title():
    print('**********************************************************************')
    print('* SharePro for joint fine-mapping and GxE analysis                   *')
    print('* Version 2.0.0                                                      *')
    print('* (C) Wenmin Zhang (wenmin.zhang@mail.mcgill.ca)                     *')
    print('**********************************************************************')


class SparseReg(object):
    def __init__(self, P, K, N, sigma, varb):
        """initialize and set parameters"""
        self.p = P
        self.k = K
        self.n = N
        self.sigma = sigma
        self.y_tau = 1.0 / (1 - varb)
        self.priorb = 1.0 / varb
        self.postn = self.n * self.y_tau + self.priorb
        self.lamb = self.priorb / self.postn
        self.beta_mu = np.zeros([self.p, self.k])
        self.delta = np.zeros((self.p, self.k))
        self.u1 = np.log(self.sigma / (1 - self.sigma)) + 0.5 * np.log(self.lamb)

    def infer_q_beta(self, b_hat_s, R_s, GAM, k):
        """perform variational updates for the k-th effect with beta and c"""
        idxall = list(range(self.k))
        idxall.remove(k)
        beta_all_k = (GAM[:, idxall] * self.beta_mu[:, idxall] * self.delta[:, idxall]).sum(axis=1)
        res_beta = b_hat_s - np.dot(R_s, beta_all_k)
        self.beta_mu[:, k] = (1 - self.lamb) * res_beta
        u = self.u1 + 0.5 * self.beta_mu[:, k] ** 2 * self.postn
        self.delta[:, k] = expit(u)
        return u

    def get_elbo(self, b_hat_s, R_s, GAM):
        """get elbo with beta and c"""
        beta_all = (GAM * self.delta * self.beta_mu).sum(axis=1)
        ll1 = self.y_tau * np.dot(beta_all, self.n * b_hat_s)
        ll2 = - 0.5 * self.y_tau * (
            (((GAM * self.delta * (self.beta_mu ** 2 + 1 / self.postn)).sum(axis=1) * self.n).sum()))
        W = GAM * self.delta * self.beta_mu
        WtRW = np.dot(np.dot(W.transpose(), self.n * R_s), W)
        ll3 = - 0.5 * self.y_tau * (WtRW.sum() - np.diag(WtRW).sum())
        ll = ll1 + ll2 + ll3
        mklbeta = - 0.5 * (GAM * self.delta * (
                self.priorb * (self.beta_mu ** 2 + 1 / self.postn) +
                np.log(self.postn / self.priorb) - 1)).sum()
        nozerodelta = self.delta[self.delta != 0]
        nozeroGAM = GAM[self.delta != 0]
        noonedelta = self.delta[self.delta != 1]
        nooneGAM = GAM[self.delta != 1]
        mkldelta = (nozeroGAM * nozerodelta * np.log(self.sigma / nozerodelta)).sum() + (
                nooneGAM * (1 - noonedelta) * np.log((1 - self.sigma) / (1 - noonedelta))).sum()
        return ll, mklbeta, mkldelta


class SharePro(object):
    def __init__(self, P, K, num, sigma, Nlst, varblst):
        """initialize and set parameters"""
        self.p = P
        self.k = K
        self.num = num
        self.sigma = sigma
        self.prior_pi = np.ones((self.p,)) * (1 / self.p)
        self.gamma = np.zeros((self.p, self.k))
        self.usum1 = self.num * np.log(1 - self.sigma**(1/self.num)) + np.log(self.prior_pi)
        self.SR = [SparseReg(P, K, Nlst[i], self.sigma**(1/self.num), varblst[i]) for i in range(self.num)]

    def infer_q_s(self, bhatlst, ld):
        """perform variational updates for s"""
        for k in range(self.k):
            idxall = list(range(self.k))
            idxall.remove(k)
            unum = np.array([j.infer_q_beta(bhatlst[i], ld, self.gamma, k) for (i,j) in enumerate(self.SR)])
            unumexp = np.log(1 + np.exp(-unum))
            uall = unum.sum(axis=0) + unumexp.sum(axis=0) + self.usum1
            self.gamma[:, k] = softmax(uall)

    def get_PIP(self):
        PIP_overall = self.gamma.max(axis=1).round(4) # causal in either
        return PIP_overall

    def get_effect(self, ld, cthres, pthres):
        vidx = np.argsort(-self.gamma, axis=1)
        matidx = np.argsort(-self.gamma, axis=0)
        mat_eff = np.zeros((self.p, self.k))  # effective gamma
        for p in range(self.p):
            mat_eff[p, vidx[p, 0]] = self.gamma[p, vidx[p, 0]]
        csum = mat_eff.sum(axis=0).round(2)
        logging.info(f"Attainable coverage for effect groups: {csum}")
        eff = {}
        eff_gamma = {}
        eff_purity = {}
        for k in range(self.k):
            if csum[k] > cthres:
                p = 0
                while np.sum(mat_eff[matidx[0:p, k], k]) < cthres * csum[k]:
                    p = p + 1
                cidx = matidx[0:p, k].tolist()
                purity = abs(ld[np.ix_(cidx, cidx)]).min()
                if purity > pthres:
                    eff[k] = cidx
                    eff_gamma[k] = mat_eff[eff[k], k].round(4).tolist()
                    eff_purity[k] = purity
        return eff, eff_gamma, eff_purity

    def update_pi(self, new_pi):
        self.prior_pi = new_pi

    def get_elbo(self, bhatlst, ld):
        """get elbo with s"""
        llsum, mklbetasum, mkldeltasum = zip(*[j.get_elbo(bhatlst[i], ld, self.gamma) for (i,j) in enumerate(self.SR)])
        gammaterm1 = (self.gamma * np.tile(self.prior_pi.reshape(-1, 1), (1, self.k))).sum()
        gammaterm2 = (self.gamma[self.gamma != 0] * np.log(self.gamma[self.gamma != 0])).sum()
        mklgamma = gammaterm1 - gammaterm2
        ll = sum(llsum)
        mklbeta = sum(mklbetasum)
        mkldelta = sum(mkldeltasum)
        elbo = ll + mklbeta + mkldelta + mklgamma
        return ll, mklbeta, mkldelta, mklgamma, elbo

    def train(self, bhatlst, ld, maxite, eps, ubound):
        loss = 100000.0
        for ite in range(maxite):
            self.infer_q_s(bhatlst, ld)
            ll, mklbeta, mkldelta, mklgamma, elbo = self.get_elbo(bhatlst, ld)
            logging.info('*' * 70)
            logging.info('Iteration-->{} . Likelihood: {:.1f} . KL_b: {:.1f} . KL_c: {:.1f} . KL_s: {:.1f} . ELBO: {:.1f}'.format(ite, ll, mklbeta, mkldelta, mklgamma, elbo))
            if abs(elbo - loss) < eps:
                converged = True
                break
            if ite == (maxite - 1) or elbo > ubound:
                print("Detected mismatch between summary statistics and LD matrix.")
                converged = False
                break
            loss = elbo
        return converged

def get_eff_maxld(eff, ld):
    idx = [i[0] for i in eff.values()]
    if len(eff)>1:
        maxld = np.abs(np.tril(ld[np.ix_(idx,idx)],-1)).max()
    else:
        maxld = 0
    return maxld

def adaptive_train(bhatlst, ld, Nlst, sigma, K, maxite, eps, ubound, cthres, pthres):
    P = bhatlst[0].shape[0]
    num = len(Nlst)
    varb = [(i**2).max() for i in bhatlst]
    Nmaxlst = [i.max() for i in Nlst]
    model = SharePro(P, K, num, sigma, Nmaxlst, varb)
    mc = model.train(bhatlst, ld, maxite, eps, ubound)
    while not mc:
        varb = [i*0.1 for i in varb]
        model = SharePro(P, K, num, sigma, Nmaxlst, varb)
        mc = model.train(bhatlst, ld, maxite, eps, ubound)
    eff, eff_gamma, eff_purity = model.get_effect(ld, cthres, pthres)
    PIP_overall = model.get_PIP()
    return eff, eff_gamma, eff_purity, PIP_overall

def get_bhat(beta, se, n):
    zscore = beta / se
    bhat = zscore / np.sqrt(zscore**2 + n)
    return bhat.round(4)

def parse_args():
    parser = argparse.ArgumentParser(description='SharePro Commands:')
    parser.add_argument('--z', type=str, default=None, nargs='+', help='path to matched summary statisticis', required=True)
    parser.add_argument('--ld', type=str, default=None, help='path to matched ld matrix', required=True)
    parser.add_argument('--save', type=str, default=None, help='path to save results', required=True)
    parser.add_argument('--sigma', type=float, default=0.99, help='prior shared prob')
    parser.add_argument('--K', type=int, default=10, help='largest number of causal signals')
    parser.add_argument('--maxite', type=int, default=100, help='max number of iterations')
    parser.add_argument('--eps', type=float, default=1e-2, help='convergence criterion')
    parser.add_argument('--ubound', type=int, default=1e10, help='upper bound for inconvergence')
    parser.add_argument('--cthres', type=float, default=0.95, help='attainable coverage threshold for effect groups')
    parser.add_argument('--pthres', type=float, default=0.8, help='purity threshold for effect groups')
    args = parser.parse_args()
    return args

def print_args(args):
    for arg in vars(args):
        logging.info(f"{arg}: {getattr(args, arg)}")

def main(args):
    zfile = [pd.read_csv(i, sep='\s+') for i in args.z]
    ld = pd.read_csv(args.ld, sep='\s+', header=None).values
    bhatlst = [get_bhat(i['BETA'].values, i['SE'].values, i['N'].values) for i in zfile]
    Nlst = [i['N'].values for i in zfile]
    eff, eff_gamma, eff_purity, PIP_overall = adaptive_train(bhatlst, ld, Nlst, args.sigma, args.K, args.maxite, args.eps, args.ubound, args.cthres, args.pthres)
    df_res = pd.DataFrame({'SNP': zfile[0]['SNP']})
    df_res['PIP'] = PIP_overall
    df_res['cs'] = 0
    df_res['cs_variants'] = 'NA'
    df_res['cs_purity'] = 0.0
    cstop = []
    for m in range(len(zfile)):
        for n in range(m + 1, len(zfile)):
            pdiff = np.nan_to_num(chi2.logsf((zfile[m]['BETA'] - zfile[n]['BETA'])**2 / (zfile[m]['SE']**2 + zfile[n]['SE']**2), 1) / -np.log(10))
            df_res['nlog10pGxE_{}_{}'.format(m,n)] = ['{:.4f}'.format(val) for val in pdiff]
            df_res['nlog10pdiff_{}_{}'.format(m,n)] = 'NA'
    for e,val in eff.items():
        mcs_idx = [zfile[0]['SNP'][j] for j in val]
        logging.info(f'The {e}-th effect group contains effective variants:')
        logging.info(f'causal variants: {mcs_idx}')
        logging.info(f'variant probabilities for this effect group: {eff_gamma[e]}')
        logging.info(f'variant purity for this effect group: {eff_purity[e]}')
        df_res.loc[val[0], 'cs'] = e+1
        cstop.append(val[0])
        df_res.loc[val[0], 'cs_variants'] = '/'.join(mcs_idx)
        df_res.loc[val[0], 'cs_purity'] = eff_purity[e]
        for m in range(len(zfile)):
            for n in range(m + 1, len(zfile)):
                df_res.loc[val[0], 'nlog10pdiff_{}_{}'.format(m,n)] = df_res.loc[val[0], 'nlog10pGxE_{}_{}'.format(m,n)]
                #logging.info('negative log GxE p-value for this effect group: {}'.format(df_res.loc[val[0], 'nlog10pGxE_{}_{}'.format(m,n)].values))
    for i in range(len(zfile)):
        pval = np.nan_to_num(chi2.logsf((zfile[i]['BETA']/zfile[i]['SE'])**2, 1)/-np.log(10))
        df_res['nlog10p{}'.format(i)] = ['{:.4f}'.format(val) for val in pval]
    df_res.to_csv('{}.sharepro.gxe.txt'.format(args.save), sep='\t', header=True, index=False)

if __name__ == '__main__':
    args = parse_args()
    logging.basicConfig(filename='{}.sharepro.gxe.log'.format(args.save), level=logging.INFO, filemode='w', format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    print_args(args)
    main(args)
