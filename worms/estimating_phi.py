import numpy as np
import pandas as pd
from scipy.special import gamma, gammaln
from scipy.optimize import minimize

def safe_ln(x, minval=0.0000000001):
    return np.log(x.clip(min=minval))

def log_nb(count, mu, k):
    term1 = gammaln(k + count) - gammaln(k) - gammaln(count+1)
    term2 = k * np.log((k) / (k + mu))
    term3 = count * np.log(mu / (mu + k))
    return term1 + term2 + term3

groups=[]

with open('expset.txt') as handle:
    for line in handle.readlines():
        groups.append(line.split())
        

data = pd.read_csv('spectrum.csv')

for group in groups:
    foo = data.loc[data['name'].isin(group)].values
    #foo=foo[:,1:].sum(axis=0)
    print foo
    break
    
M = np.random.uniform(size=[645*11])
K = [100]
M = np.concatenate([M, K], axis=0)

const=0
for group in groups:
    const+=len(group)
print const*11

def demon_loss(M, args):
    data = args[0]
    groups= args[1]
    loss = 0.0
    K = M[-1]
    M = np.reshape(M[:-1], newshape=(645, 11))
    for group, i in zip(groups, xrange(len(groups))):
        counts = data.loc[data['name'].isin(group)].values
        counts=np.asarray(counts[:,1:].sum(axis=0), dtype=np.float64)
        partial_loss = log_nb(counts, M[i], K)
        loss+=np.sum(partial_loss)
    
    return loss

result = minimize(demon_loss, M, args=[data, groups], bounds=[(0.00001, None),]*(645*11+1), options={'maxiter' : 100, 'disp'
 : True})


