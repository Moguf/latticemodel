#! /usr/bin/evn python

class binary:
    def __init__(self):
        self.state=0
        self.left=0
        self.right=0


class SolveHPwithBDD:
    def __init__(self):
        self.seq=[]
        self.root=binary()

        
    def makeBDD(self,binary):
        binary.left=binary()
        binary.right=binary()

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
        
    def show(self,binary):
        print binary.state
        if binary.left==0 and binary.right==0:
            return 0
        self.show(binary.left)
        self.show(binary.right)

    
if __name__=="__main__":
    tclass=SolveHPwithBDD()
    #tclass.read("HPPHPPHPHPPHPHPHHH")
    tclass.read("HPPH") 
    tclass.makeBDD()
    tclass.show(tclass.root)
