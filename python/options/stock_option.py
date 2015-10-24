__author__ = 'christopheralexander_nmp'

import math

class StockOption(object):
    def __init__(self, S0, K, r, T, N, params):
        self.S0 = S0
        self.K = K
        self.r = r
        self.T = T
        self.N = max(1, N) # N needs at least one time step
        self.STs = None # stock prices tree

        """ Optional params used by child classes"""
        self.pu = params.get("pu", 0) # Prob of up state
        self.pd = params.get("pd", 0) # prob of down state
        self.div = params.get("div", 0) # dividend yield
        self.sigma = params.get("sigma", 0) # volatility
        self.is_call = params.get("is_call", True) # call or put
        self.is_european = params.get("is_eu", True) # euro or american

        """ Computed values """
        self.dt = T / float(N) # single time step, in years
        self.df = math.exp(-(r - self.div) * self.dt) # discount factor: e^-rt