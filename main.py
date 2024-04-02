import numpy as np
import numpy.random as rnd
import math
import matplotlib .pyplot as plt

def firstq():
    sm = 0
    N = 1000000
    tau = 1
    n = 3
    for i in range(1, N):
        r = rnd.random()
        t = -tau * math.log(1 - r)
        sm += t ** n
    moment = sm/N
    print(moment)

def secondq(N = 10000):
    xvalues = np.empty(N)
    yvalues = np.empty(N)
    sigma = 1
    for i in range(1, N):
        s = rnd.random()
        t = rnd.random()
        r = sigma * np.sqrt(-2 * np.log(1 - t))
        phi = s * 2 * np.pi
        xvalues[i] = r * np.cos(phi)
        yvalues[i] = r * np.sin(phi)

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.hist(xvalues, bins = 50, range = [-4.0, 4.0], density = True)
    ax2 = fig.add_subplot(122)
    ax2.hist(yvalues, bins = 50, range = [-4.0, 4.0], density = True)
    plt.show()

def f(x):
    return np.exp(x**2)-1
def W(x):
    return 3*x**2

def thirdq(N = 10000):
    sum_f = 0.
    sum_f2 = 0.
    for i in range(1, N):
        r = rnd.random()
        sum_f += f(r)
        sum_f2 += f(r)**2
    av_f = sum_f/N
    av_f2 = sum_f2/N
    result = av_f
    error = np.sqrt((av_f2-av_f**2)/N)
    return result

def forthq(N = 10000):
    x = rnd.random()
    Wx = 0
    sum_f = 0.
    step = 0.4

    for i in range(1, N):
        r = rnd.random()
        xtry = x + 2 * step * (r - 0.5)
        if 0.0 <= xtry and xtry < 1.0:
            Wxtry = W(xtry)
            if Wxtry > Wx:
                x = xtry
                Wx = Wxtry
            else:
                r = rnd.random()
                if Wxtry > Wx * r:
                    x = xtry
                    Wx = Wxtry
        sum_f += (f(x)/W(x))
    av_f = sum_f/N
    result = av_f
    return result

def ran():
    Nmin =1
    Nmax =1e4
    Nsamples = 9
    Ns = np.geomspace(Nmin, Nmax, Nsamples)
    print(Ns)
    Ns = np.array(sorted(set(Ns.astype(int))))
    #print(Ns)
    fig, axs = plt.subplots(1, 2)
    axs[0].plot(Ns, '.-')
    axs[1].plot(Ns, '.-')
    axs[1].set_yscale('log')
    plt.show()

if __name__ == "__main__":
    #firstq()
    #secondq()
    #ran()
    print(thirdq())
    print(forthq())