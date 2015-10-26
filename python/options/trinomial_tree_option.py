__author__ = 'christopheralexander'
import math
import numpy as np

from binomial_tree_option import BinomialTreeOption

class TrinomialTreeOption(BinomialTreeOption):

    def _setup_parameters_(self):
        """ Required calculations for model """
        self.u = math.exp(self.sigma * math.sqrt(2.0 * self.dt))
        self.d = 1.0 / self.u
        self.m = 1.0

        self.qu = ((math.exp((self.r - self.div) * self.dt / 2.0) - math.exp(-self.sigma * math.sqrt(self.dt / 2.0))) /
                   (math.exp(self.sigma * math.sqrt(self.dt / 2.0)) - math.exp(-self.sigma * math.sqrt(self.dt / 2.0)))) ** 2

        self.qd = ((math.exp(self.sigma * math.sqrt(self.dt / 2.0)) - math.exp((self.r - self.div) * self.dt / 2.0)) /
                   (math.exp(self.sigma * math.sqrt(self.dt / 2.0)) - math.exp(-self.sigma * math.sqrt(self.dt / 2.0)))) ** 2

        self.qm = 1 - self.qu - self.qd

    def _initialize_stock_price_tree_(self):
        """ Initialize a 2-D tree at t=0 """
        self.STs = [np.array([self.S0])]

        print("initial STs ", self.STs)

        for i in range(self.N):
            prev_nodes = self.STs[-1]
            self.ST = np.concatenate((prev_nodes * self.u, [prev_nodes[-1] * self.m, prev_nodes[-1] * self.d]))
            self.STs.append(self.ST)
            print("now STs are ", self.STs)

    def _traverse_tree_(self, payoffs):
        """ Traverse tree backwards """
        for i in reversed(range(self.N)):
            payoffs = (payoffs[:-2] * self.qu + payoffs[1:-1] * self.qm + payoffs[2:] * self.qd) * self.df

            if not self.is_european:
                payoffs = self.__check_early_exercise__(payoffs, i)

        return payoffs