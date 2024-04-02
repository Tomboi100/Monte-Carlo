import numpy as np
import numpy.random as rnd
import math
import matplotlib .pyplot as plt

# def MS(r):
#     endpoint = r.size()
#     startpoint = r[0]
#     MS = np.sum((endpoint - startpoint)**2)
#     return MS
#
# def ROG(r):
#     N = r.size()
#     rcmsum = 0  # rcmsum is a vector
#     for i in range(0, N):
#         rcmsum += r[i,:]  # addition of vectors
#     rcm = rcmsum / (N+1)  # rcm is a vector
#     rg2sum = 0  # rg2sum is a number
#     for i in range(0, N):
#         rg2sum += np.sum((r[i,:] - rcm) ** 2)  # length of difference vector squared
#     rg2 = rg2sum / (N + 1)
#     return rg2

def get_re2(r):
    return ((r[-1,:]-r[0,:])**2).sum()

def get_rcm(r):
    return r.mean(axis=0)

def get_rgt(r):
    rcm = get_rcm(r)
    return ((r-rcm)**2).sum()/r.shape[0]

def GenRWstep(d):
    u = np.zeros(d)
    drnd = rnd.randint(d)
    u[drnd] = 2 * rnd.randint(2)-1
    return u

def GenRW(r, niter):
    # niter means n
    N, d = r.shape
    N -= 1
    sum_re2 = 0
    sum_rg2 = 0
    for i in range(0, niter+1):
        r[0, :] = 0  # r[0] is a vector
        for i in range(0, N):
            r[i+1,:] = r[i,:] + GenRWstep(d)  # r[i] and u are vectors
        #calculate re2  # end point distance squared
        sum_re2 += get_re2(r)
        #calculate rg2  # radius of gyration squared
        sum_rg2 += get_rgt(r)

        inverseniter = 1/niter
        av_re2 = sum_re2 *inverseniter
        av_rg2 = sum_rg2 *inverseniter
        return av_re2, av_rg2

if __name__ == "__main__":
    print("start main")
    N = 10
    n = 1000
    d = 2
    r = np.empty([N+1, d])
    av_re2, av_rg2 = GenRW(r, n)
    print(GenRW(r, n))

    # fig, ax = plt.subplots()
    # ax.plot(av_re2, av_rg2)
    # ax.set(xlabel='time (s)', ylabel='av_rg2', title='PLOT')
    # ax.grid()
    # plt.show()