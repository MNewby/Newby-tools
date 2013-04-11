#    Copyright (C) 2004, Maxime Biais <maxime@biais.org>
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# $Id: pso.py,v 1.2 2004/10/07 12:30:42 max Exp $

# Adapted by willeb to use 2 dimensional functions

from random import seed 
from random import uniform
import subprocess

class PSO:
    def __init__(self, pop_size, min, max, phi, phi2, lr, max_iter, func):
        self.func = func
        self.pop = []
        # 0,1,2,3,4,5,6: position, 7,8,9,10,11,12,13: velocity, 14: fitness
        self.min = min
        self.max = max
	seed()

        for i in xrange(pop_size):
            self.pop.append([uniform(5,50), uniform(-300,300), uniform(-300,300), uniform(-300,300), uniform(100, 150), uniform(0.5,1.5), uniform(5,20), uniform(-0.1, 0.1), uniform(-1, 1), uniform(-1, 1), uniform(-1, 1), uniform(-1, 1), uniform(-0.01,0.01), uniform(-0.1,0.1), 0])
        self.evaluate()
        self.gdest = self.pop[0]
	self.pop = []
        for i in xrange(pop_size):
            self.pop.append([uniform(5,50), uniform(-300,300), uniform(-300,300), uniform(-300,300), uniform(100, 150), uniform(0.5,1.5), uniform(5,20), uniform(-0.1, 0.1), uniform(-1, 1), uniform(-1, 1), uniform(-1, 1), uniform(-1, 1), uniform(-0.01,0.01), uniform(-0.1,0.1), 0])
        self.evaluate()
        self.pdest = self.pop[0]
        self.phi = phi
        self.phi2 = phi2
        self.lr = lr
        self.max_iter = max_iter
    
    def update_velocity(self):
        for i in self.pop:
            i[7] = self.lr * i[7] + uniform(0, self.phi) * (self.pdest[0] \
                    - i[0]) + uniform(0, self.phi2) * (self.gdest[0] - i[0])
            i[8] = self.lr * i[8] + uniform(0, self.phi) * (self.pdest[1] \
                    - i[1]) + uniform(0, self.phi2) * (self.gdest[1] - i[1])
            i[9] = self.lr * i[9] + uniform(0, self.phi) * (self.pdest[2] \
                    - i[2]) + uniform(0, self.phi2) * (self.gdest[2] - i[2])
            i[10] = self.lr * i[10] + uniform(0, self.phi) * (self.pdest[3] \
                    - i[3]) + uniform(0, self.phi2) * (self.gdest[3] - i[3])
            i[11] = self.lr * i[11] + uniform(0, self.phi) * (self.pdest[4] \
                    - i[4]) + uniform(0, self.phi2) * (self.gdest[4] - i[4])
            i[12] = self.lr * i[12] + uniform(0, self.phi) * (self.pdest[5] \
                    - i[5]) + uniform(0, self.phi2) * (self.gdest[5] - i[5])
            i[13] = self.lr * i[13] + uniform(0, self.phi) * (self.pdest[6] \
                    - i[6]) + uniform(0, self.phi2) * (self.gdest[6] - i[6])

    def evaluate(self):
        for i in self.pop:
            i[14] = self.func(i[0],i[1],i[2],i[3],i[4],i[5],i[6])

    def move(self):
        for i in self.pop:
            i[0] += i[7]
            i[1] += i[8]
            i[2] += i[9]
            i[3] += i[10]
            i[4] += i[11]
            i[5] += i[12]
            i[6] += i[13]

            if (i[0] < 5 or i[0] > 50) or (i[1] < -300 or i[1] > 300) or (i[2] < -300 or i[2] > 300) or (i[3] < -300 or i[3] > 300) or (i[4] < 100 or i[4] > 150) or (i[5] < 0.5 or i[5] > 1.5) or (i[6] < 5 or i[6] > 15):
            		i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], i[11], i[12], i[13], i[14] = uniform(5,50), uniform(-300,300), uniform(-300,300), uniform(-300,300), uniform(100, 150), uniform(0.5, 1.5), uniform(5, 20), uniform(-0.1, 0.1), uniform(-1, 1), uniform(-1, 1), uniform(-1, 1), uniform(-1, 1), uniform(-0.01, 0.01), uniform(-0.1, 0.1), 0
    def __cmp_by_fitness(self, a, b):
        return cmp(a[14], b[14])
    
    def run(self, update_func=False):
        for i in xrange(self.max_iter):
            if update_func:
                update_func()
            self.update_velocity()
            self.move()
            self.evaluate()
            self.pop.sort(self.__cmp_by_fitness, reverse=0)
            self.pdest = self.pop[0]
            print self.pdest[14], self.gdest[14]
            if self.pdest[14] < self.gdest[14]:
                self.gdest[0] = self.pdest[0]
                self.gdest[1] = self.pdest[1]
                self.gdest[2] = self.pdest[2]
                self.gdest[3] = self.pdest[3]
                self.gdest[4] = self.pdest[4]
                self.gdest[5] = self.pdest[5]
                self.gdest[6] = self.pdest[6]
                self.gdest[7] = self.pdest[7]
                self.gdest[8] = self.pdest[8]
                self.gdest[9] = self.pdest[9]
                self.gdest[10] = self.pdest[10]
                self.gdest[11] = self.pdest[11]
                self.gdest[12] = self.pdest[12]
                self.gdest[13] = self.pdest[13]
                self.gdest[14] = self.pdest[14]
            

    def __str__(self):
        ret = ""
        for i in self.pop:
            ret += str(i) + "\n"
        return ret
                                 
def test():
    import math
    # func = lambda x:math.cos(x*math.sin(x*0.3)-x) / 1.5
    #func = lambda x:math.cos(x) * math.exp(math.sin(x)) * math.sin(x)  / 1.5
    func = lambda R1, vx1, vy1, vz1, vhalo, q, d: float(subprocess.Popen([r"./func.pl", str(R1), str(vx1), str(vy1), str(vz1), str(vhalo), str(q), str(d)], stdout=subprocess.PIPE).communicate()[0])
    p = PSO(10, -4.5, 4.5, phi=0.3, phi2=0.3, lr=0.9, max_iter=1000000, func=func)

    #printer = PygamePrinter(p)
    p.run(update_func=False)
    print p
    
if __name__ == "__main__":
    test()
