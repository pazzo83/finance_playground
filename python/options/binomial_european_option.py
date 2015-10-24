__author__ = 'christopheralexander_nmp'
from stock_option import StockOption

import math
import numpy as np

class BinomialEuropeanOption(StockOption):

    def _setup_parameters_(self):
        """ Required calculations for the model """
        self.M = self.N + 1 # number of terminal nodes in the tree
        self.u = 1 + self.pu # expected value in the up state
        self.d = 1 - self.pd # Expected value in the down state
        self.qu = (math.exp((self.r - self.div) * self.dt) - self.d) / (self.u - self.d)
        self.qd = 1 - self.qu

    def _initialize_stock_price_tree_(self):
        # Initialize terminal price nodes to zero
        self.STs = np.zeros(self.M)

        # Calculate expected stock prices for each node
        for i in range(self.M):
            self.STs[i] = self.S0 * (self.u ** (self.N - i)) * (self.d ** i)

    def _initialize_payoffs_tress_(self):
        # get payoffs when the option expires at terminal nodes
        payoffs = np.maximum(0, (self.STs - self.K) if self.is_call else (self.K - self.STs))
        return payoffs

    def _traverse_tree_(self, payoffs):
        # starting from time the option expires, traverse backwards and calculate the discounted payoffs at each node
        for _ in range(self.N):
            payoffs = (payoffs[:-1] * self.qu + payoffs[1:] * self.qd) * self.df

        return payoffs

    def __begin_tree_traversal__(self):
        payoffs = self._initialize_payoffs_tress_()
        return self._traverse_tree_(payoffs)

    def price(self):
        """
         Pricing implementation
        :return first node payoff
        """
        self._setup_parameters_()
        self._initialize_stock_price_tree_()
        payoffs = self.__begin_tree_traversal__()

        return payoffs[0]