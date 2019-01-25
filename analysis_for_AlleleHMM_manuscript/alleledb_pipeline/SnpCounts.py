
import gc, os, sys, string, re, pdb, scipy.stats
import getNewer1000GSNPAnnotations, Bowtie, cPickle, utils

USAGE="Usage"

if __name__=='__main__':

    if len(sys.argv) != 5:
        print USAGE % sys.argv[0]
        sys.exit(-1)

    snpfile=sys.argv[1]
    readfile=sys.argv[2]
    maptmplt=sys.argv[3]
    OUTFP = open(sys.argv[4], 'w')

    gc.disable()

    pat=re.compile('(\w+)_([mp]aternal)')

    mappers=utils.makeMappers(maptmplt)

    ref_1000G=getNewer1000GSNPAnnotations.Handler(snpfile, None, 'PAT', hasHeader=True, onlyHets=True)
    pat_1000G=getNewer1000GSNPAnnotations.Handler(snpfile, mappers, 'PAT', hasHeader=True, onlyHets=True)
    mat_1000G=getNewer1000GSNPAnnotations.Handler(snpfile, mappers, 'MAT', hasHeader=True, onlyHets=True)

    mp={'paternal': pat_1000G, 'maternal': mat_1000G}

    counts={}
    ksvals={}
    
    for ident, strand, chrom, start, seq, qual, dummy in Bowtie.parse(readfile):
        start=int(start)+1 # bowtie is 0 based, everything else here is 1 based.
        chrm, parent = utils.parseChrom(chrom)
        readlen=len(seq)
        poss, snps=mp[parent].getAnnotations(chrm, start, start+readlen-1)
        for pos, snp in zip(poss, snps): # these are the snps that are in parent coordinates
            ref_pos, mat_genotype, pat_genotype, child_genotype, mat_allele, pat_allele, typ, ref, hetSNP = snp
            base=seq[pos-start].lower()
            if chrm not in counts:
                counts[chrm]={}
                ksvals[chrm]={}
            if ref_pos not in counts[chrm]:
                counts[chrm][ref_pos]={'a':0, 'c':0, 'g':0, 't':0, 'n':0}
                ksvals[chrm][ref_pos]=[]
            counts[chrm][ref_pos][base]+=1
            ksvals[chrm][ref_pos].append(float(pos-start)/readlen)

    cPickle.dump(counts, OUTFP)

