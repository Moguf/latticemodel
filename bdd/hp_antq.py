#!/usr/bin/env python

class HpAntQ:
    def __init__(self):
        self.seq=[]
        pass
    
    def read(self,seq_str):
        self.seq=self._convertSeq2Binary(seq_str)
        print self.seq


    def _convertSeq2Binary(self,seq_str):
        seq=[]
        for i in seq_str:
            if i=="H":
                seq.append(1)
            else:
                seq.append(0)
        return seq
        
    
if __name__=="__main__":
    tclass=HpAntQ()
    #tclass.read("HPPHPPHPHPPHPHPHHH")
    tclass.read("HPPH") 
            

