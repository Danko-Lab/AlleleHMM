
import glob, sys, os.path
import MergeBowtie
import Mapping2

# this is needed to account for the fact that the maternal, paternal, and ref genomes all have different coordinates now
# that indels are handled

def makeMappers(maptmplt):
    mappers={}
    cs=['chr%s' % str(c) for c in range(1,23)] + ['chrX', 'chrY', 'chrM']
 
    for c in cs:
        f=maptmplt % c
        if os.path.exists(f):
            mappers[c] = Mapping2.Mapping(f)
    return mappers

if __name__=='__main__':

    Fathers=sorted(glob.glob("*_AltRefFather.txt"))
    Mothers=sorted(glob.glob("*_AltRefMother.txt"))
    maptmplt=sys.argv[1]
    outfile=sys.argv[2]
    
    mappers=makeMappers(maptmplt)
    log=open('Merge.log', 'w')
    ofp=open(outfile, 'w')
    
    for f, m in zip(Fathers, Mothers):
        print >>log, 'Merging %s %s' % (f, m)
        log.flush()
        assert (f.replace('AltRefFather', 'AltRefMother') == m)
        MergeBowtie.process(open(f), open(m), mappers, ofp, log)

