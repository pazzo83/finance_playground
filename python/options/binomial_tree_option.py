__author__ = 'christopheralexander_nmp'

from stock_option import StockOption
import math
import numpy as np

class BinomialTreeOption(StockOption):

    def _setup_parameters_(self):
        self.u = 1 + self.pu # Expected value in the up state
        self.d = 1 - self.pd # Expected value in the down state
        self.qu = (math.exp((self.r - self.div) * self.dt) - self.d) / (self.u - self.d)
        self.qd = 1 - self.qu

    def _initialize_stock_price_tree_(self):
        # initialize a 2D tree at T=0
        self.STs = [np.array([self.S0])]

        print("initialize STs")
        print(self.STs)

        # Simulate the possible stock prices path
        for i in range(self.N):
            prev_branches = self.STs[-1]
            print('prev branches ', prev_branches)
            st = np.concatenate((prev_branches * self.u, [prev_branches[-1] * self.d]))
            print('st ', st)
            self.STs.append(st)
            print("STs now ", self.STs)

    def _initialize_payoffs_tree_(self):
        return np.maximum(0, (self.STs[self.N] - self.K) if self.is_call else (self.K - self.STs[self.N]))

    def __check_early_exercise__(self, payoffs, node):
        early_ex_payoff = (self.STs[node] - self.K) if self.is_call else (self.K - self.STs[node])
        print('early ex ', early_ex_payoff)

        return np.maximum(payoffs, early_ex_payoff)

    def _traverse_tree_(self, payoffs):
        print('payoffs initial ', payoffs)
        for i in reversed(range(self.N)):
            # The payoffs from NOT exercising the option
            payoffs = (payoffs[:-1] * self.qu + payoffs[1:] * self.qd) * self.df
            # print('payoffs now ', payoffs)

            if not self.is_european:
                payoffs = self.__check_early_exercise__(payoffs, i)
                # print('payoffs american ', payoffs)

        return payoffs

    def __begin_tree_traversal__(self):
        payoffs = self._initialize_payoffs_tree_()
        return self._traverse_tree_(payoffs)

    def price(self):
        self._setup_parameters_()
        self._initialize_stock_price_tree_()
        payoffs = self.__begin_tree_traversal__()

        return payoffs[0]