import sys

""" display:  PREFIX [#########=     ] 50% SUFFIX """


"""http://thelivingpearl.com/2012/12/31/creating-progress-bars-with-python/"""



class Progressbar:
    def __init__(self, steps=100, prefix="Progress:", suffix=None, symbol="#", 
                active="=", brackets="[]", percent=True, size=30):
        self.steps = steps
        self.prefix = prefix
        self.suffix = suffix
        self.symbol = symbol
        self.active = active
        self.brackets = brackets
        self.percent = percent
        self.size = size  #Visual size in spaces
        self.state = 0
        self.printbar()
    def printbar(self, end=1):
        out1 = "{0} {1}".format(self.prefix,self.brackets[0])
        out2 = (self.symbol)*self.state + " "*(self.size-self.state)
        out3 = "{0}".format(self.brackets[1])
        if self.percent==True:  out3 = out3+" {0}% ".format(int(100.0*(float(self.state)/float(self.size))))
        if self.suffix != None:  out2 = out2 + suffix
        out = out1+out2+out3
        width = len(out)+2
        print out,
        if end==1:  print '\b'*width,
        sys.stdout.flush()
    def updatebar(self, prog, suff=None):
        if prog >= 1.0:  self.endbar(); return
        complete = int(self.size*prog)
        if complete > self.state:  
            self.state = complete
            self.printbar()
    def updatesuffix(self):
        return -1
    def endbar(self):
        self.state = self.size
        self.printbar(0)
        print "  Done!"
    
        

if __name__=="__main__":
    import time
    testbar = Progressbar()
    #testbar.displaybar()
    for i in range(101):
        time.sleep(0.1)
        testbar.updatebar(i/100.0)
    print "Test Successful"
