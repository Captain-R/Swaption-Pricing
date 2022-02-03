from scipy.stats import norm
import numpy as np
import math
from scipy.optimize imprt

fsolve
from scipy.interpolate import interp1d
from scipy.integrate import quad

rv = norm()


def PV(f, t, T):
    p = math.exp(-f(t, T) * (T - t))  # price at time t of zero coupon bond with unit face value and maturity T
    return p


def Annuity(f, m, n, t, T):
    sum1 = 0
    for i in range(m * n):
        k = PV(f, t, T + ((i + 1) / m))
        sum1 += k
    sum1 = sum1 / m
    return sum1


def B(a, t, T):
    v = (1 / a) * (1 - math.exp(-a * (T - t)))
    return v


def A(a, s, t, T, Pm, fm):
    v = (Pm(T) / Pm(t)) * math.exp(
        B(a, t, T) * fm(0) - (s * s / (4 * a)) * (1 - math.exp(-2 * a * t)) * (B(a, t, T) * B(a, t, T)))
    return v


def P(model, a, s, t, T, Pm, fm):
    v = A(a, s, t, T, Pm, fm) * math.exp(-B(a, t, T) * model.rate(t))
    return v


def ZBP(model, a, sigma, t, T, Pm, fm, S, X):
    d = (1 - math.exp(-2 * a * (T - t))) / (2 * a)
    sigma_p = sigma * math.sqrt(d) * B(a, T, S)
    d1 = P(model, a, sigma, t, S, Pm, fm) / P(model, a, sigma, t, T, Pm, fm)
    h = (1 / sigma_p) * math.log(d1 / X) + sigma_p / 2
    v = X * P(model, a, sigma, t, T, Pm, fm) * rv.cdf(-h + sigma_p) - P(model, a, sigma, t, S, Pm, fm) * rv.cdf(-h)
    return v


def Blackswaprate(f, T, n, m, LIBOR):
    fixed = 0
    floating = 0
    for i in range(m * n):
        fixed += PV(f, T, T + (i + 1) / m)
        floating += LIBOR(T, T + (i + 1) / m) * PV(f, T, T + (i + 1) / m)
    r = floating / fixed
    return r


# f(t,T) is the interest rate for period starting at t and ending at T
# L is the notional of the swap
# r is the fixed rate of the swap
# n is the duration of swap in years
# m is the number of yeaarly payments of the swap
# rt is forward swap rate with maturity T  at time=t(if at-the-money then rt=r0=r)
# T is the time in years when swap starts
# s is the sigma of the swaption
# t is the time at which value is desired
def BlackPayerValue(f, L, r, n, m, rt, T, s, t, LIBOR):
    A = Annuity(f, m, n, t, T)
    d1 = (math.log(rt / r) + s * s * (T - t) / 2) / (s * math.sqrt(T - t))
    d2 = d1 - s * math.sqrt(T - t)

    V = L * A * (rt * rv.cdf(d1) - r * rv.cdf(d2))

    R = Blackswaprate(f, T, n, m, LIBOR)

    Vega = L * A * R * math.sqrt(T) * rv.pdf(d1)

    return V, Vega


def BachelierPayerValue(f, L, r, n, m, rt, T, s, t)
    A = Annuity(f, m, n, t, T)
    d = (rt - r) / (s * math.sqrt(T - t))

    V = L * A * s * math.sqrt(T - t) * (d * rv.cdf(d) + rv.pdf(d))
    return V


# r is the rate function (yet to be modelled)
# L is the notional of the swap
# k is the fixed rate of the swap
# n is the duration of swap in years
# m is the number of yeaarly payments of the swap
# Pm is the zero coupon curve obtained from the market
# fm is the market instantaneous forward rates curve
# T is the time in years when swap starts
# sigma is the sigma of the swaption
# a is the constant according to HW model
# t is the time at which value is desired
class HWJamshidian(a, L, k, n, m, Pm, fm, T, sigma):
    def __init__(self):
        self.m = m
        self.a = a
        self.n = n
        self.k = k
        self.T = T
        self.L = L
        self.Pm = Pm
        self.fm = fm
        self.sigma = sigma
        self.c = []
        for i in range(m * n - 1):
            self.c.append(k / m)
        self.c.append(1 + k / m)

    def rate(self, t):
        size = 1000
        alpha = self.fm(t) + ((self.sigma * self.sigma) / (2 * self.a * self.a)) * (1 - math.exp(-self.a * t)) * (
                    1 - math.exp(-self.a * t))
        delta_t = t / size
        N = []
        W = [0] * (size + 1)
        x = [0] * (size + 1)
        for i in range(size):
            N.append(np.random.normal(0, 1).item())
        for i in range(len(W) - 1):
            W[i + 1] = W[i] + math.sqrt(delta_t) * N[i]
        for i in range(size):
            x[i + 1] = x[i] - self.a * x[i] * delta_t + self.sigma * (W[i + 1] - W[i])
        rt = x[size] + alpha
        return rt

    def find(self, r):
        sum1 = 0
        for i in range(self.m * self.n):
            ti = self.T + (i + 1) / self.m
            sum1 += c[i] * A(self.a, self.s, self.T, ti, self.Pm, self.fm) * math.exp(-B(self.a, self.T, ti) * r)

        return sum1 - 1


def HWJam(t, a, L, k, n, m, Pm, fm, T, sigma):
    model = HWJamshidian(a, L, k, n, m, Pm, fm, T, sigma)
    r_j = fsolve(model.find(), 1)
    X = []
    for i in range(m * n):
        ti = T + (i + 1) / m
        X.append(model.c[i] * A(a, s, T, ti, Pm, fm) * math.exp(-B(a, T, ti) * r_j))

    V = 0
    for i in range(m * n):
        ti = T + (i + 1) / m
        V += model.c[i] * ZBP(model, a, sigma, t, T, Pm, fm, ti, X[i])

    Value = V * L
    return Value


class MonteCarlo(a, L, k, n, m, Pm, fm, T, sigma, r0):
    def __init__(self):
        self.m = m
        self.a = a
        self.n = n
        self.k = k
        self.T = T
        self.L = L
        self.Pm = Pm
        self.fm = fm
        self.sigma = sigma
        self.r0 = r0

    def simulate(self):
        dt = 0.01 / self.m
        dt1 = 0.01 / self.m
        end = self.T + self.n
        rates = [self.fm(0)]
        x_current = 0
        times = [0]
        while (end - dt > 0):
            dWt = np.random.normal(0, math.sqrt(dt1)).item()
            dxt = -self.a * x_current * dt1 + self.sigma * dWt
            x_current += dxt
            alpha = self.fm(dt) + ((self.sigma * self.sigma) / (2 * self.a * self.a)) * (1 - math.exp(-self.a * dt)) * (
                        1 - math.exp(-self.a * dt))
            rates.append(x_current + alpha)
            times.append(dt)
            dt += dt1
        r = interp1d(times, rates)
        return r

    def discount(self, r, t, T):
        I = quad(r, t, T)
        d = math.exp(-I)
        return d

    def swaprate(self, r, L):  # L is the forward LIBOR rate curve
        fixed = 0
        floating = 0
        for i in range(self.m * self.n):
            fixed += self.discount(r, self.T, self.T + (i + 1) / self.m)
            floating += L(self.T, self.T + (i + 1) / self.m) * self.discount(r, self.T, self.T + (i + 1) / self.m)
        R = floating / fixed
        return R

    def payoff(self, t, r, L):
        payoff = 0
        R = self.swaprate(r, L)
        if R <= self.k:
            return payoff
        else:
            for i in range(self.m * self.n):
                payoff += self.L * (R - self.k) * self.discount(r, self.T, self.T + (i + 1) / self.m)
            return payoff * self.discount(r, t, self.T)


def HWMonteCarlo(a, L, k, n, m, Pm, fm, t, T, sigma, r0, iterations, LIBOR):
    model = MonteCarlo(a, L, k, n, m, Pm, fm, T, sigma, r0)
    avg = 0
    for i in range(iterations):
        r = model.simulate()
        avg += model.payoff(t, r, LIBOR)
    V = avg / iterations
    return V


#########################################################################################
a =  # constant in HW Model
sigma =  # constant in HW model
L =  # notional amount for underlying swap
k =  # strike price for swaption
n =  # duration of swap in years
m =  # number of yearly payments in the swap
Pm =  # a function such that Pm(T) is the market discount factor for maturity T
fm =  # a function such that fm(T) is the market instantaneous forward rate at time 0 for maturity T
T =  # maturity of swaption in years
t =  # time for which value of swaption is required in years (t < T)
r0 =  # current rate of interest (required for starting value for MC simulation)
iterations =  # number of iterations to run for MC simulation
LIBOR =  # function such that LIBOR(t,T) gives the forward LIBOR rate for period from t to T (needed to get swaprate for payoff calc. and Vega calc.) (i did not know how to get this from rate dynamics)
f =  # function such that f(t,T) gives the interest rate for period starting at t and ending at T (for Black/Bachelier Models)
rt =  # swaprate for maturity T at time t (if t=0 and at-the-money then rt=k)

BlackPrice, BlackVega = BlackPayerValue(f, L, k, n, m, rt, T, sigma, t, LIBOR)

BachelierPrice = BachelierPayerValue(f, L, k, n, m, rt, T, sigma, t)

HWJamPrice = HWJam(t, a, L, k, n, m, Pm, fm, T, sigma)

HWMCPrice = HWMonteCarlo(a, L, k, n, m, Pm, fm, t, T, sigma, r0, iterations, LIBOR)

print("Black Price is", BlackPrice)
print("Black Model Vega is", BlackVega)
print("Bachelier Price is", BachelierPrice)
print("HW Jamshidian Price is", HWJamPrice)
print("HW MC Price is", HWMCPrice)
