
import sys

class fastaReader:
    def __init__(self, fn):
        self.done=False
        self.ifile=open(fn)
        self.hdr=self.ifile.readline().rstrip()
    def __iter__(self):
        return self
    def next(self):
        if self.done:
            raise StopIteration
        body=''
        while True:
            l=self.ifile.readline().rstrip()
            if not l:
                self.done=True
            if not l or l[0]=='>':
                hdr=self.hdr
                self.hdr=l
                return (hdr, body)
            else:
                body+=l
                
class fastaWriter:
    def __init__(self, fn, linelen=60):
        self.ofile=open(fn, 'w')
        self.linelen=linelen
    def close(self):
        self.ofile.close()
    def writeFA(self, hdr, body):
        pos=0
        stop=len(body)
        self.ofile.write(hdr)
        self.ofile.write('\n')
        while pos<stop:
            self.ofile.write(body[pos:pos+self.linelen])
            self.ofile.write('\n')
            pos+=self.linelen

if __name__=='__main__':
    rdr=fastaReader(sys.argv[1])
    wrter=fastaWriter(sys.argv[2], 10)
    for hdr, body in rdr:
        wrter.writeFA(hdr, body)
    wrter.close()
