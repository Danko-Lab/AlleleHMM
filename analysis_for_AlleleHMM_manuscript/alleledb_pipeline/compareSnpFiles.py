
import sys, pdb, get1000GSNPAnnotations, getNew1000GSNPAnnotations
c='22'

slice=(0,1,2,6,7) # skip phasing for now

if __name__=='__main__':
    oldfile=sys.argv[1]
    newfile=sys.argv[2]

    h1=get1000GSNPAnnotations.Handler(oldfile)
    h2=getNew1000GSNPAnnotations.Handler(newfile, True)

    s1=set(h1.getClassPositions(c))
    s2=set(h2.getClassPositions(c))
    onlyin1=len(s1.difference(s2))
    onlyin2=len(s2.difference(s1))
    allboth=s1.intersection(s2)

    pdb.set_trace()

    print 'considering all'
    print 'only s1 %d' % onlyin1
    print 'only s2 %d' % onlyin2
    print 'common %d' % len(allboth)

    hets1=[pos for c, pos, snp in h1.getAnnotationsGenerator(c) if snp[7]]
    hets2=[pos for c, pos, snp in h2.getAnnotationsGenerator(c) if snp[7]]
    
    s1=set(hets1)
    s2=set(hets2)
    onlyin1=len(s1.difference(s2))
    onlyin2=len(s2.difference(s1))
    hetboth=s1.intersection(s2)

    print 'considering hets'
    print 'only s1 %d' % onlyin1
    print 'only s2 %d' % onlyin2
    print 'common %d' % len(hetboth)

    motherdiff=fatherdiff=childdiff=phasediff=same=0

    for pos in hetboth:
        flag=False
        pos1, snp1=h1.getAnnotations(c, pos, pos)
        pos2, snp2=h2.getAnnotations(c, pos, pos)
        snp1=snp1[0]
        snp2=snp2[0]
        if snp1[0]!=snp2[0]: 
            motherdiff+=1
            flag=True
        if snp1[1]!=snp2[1]: 
            fatherdiff+=1
            flag=True
        if snp1[2]!=snp2[2]: 
            childdiff+=1
            flag=True
        
        if not flag and snp1[3]!=None and (snp1[3]!=snp2[3] or snp1[4] != snp2[4]):
            phasediff+=1
            flag=True

        if not flag:
            same+=1


#        if flag:
#            print c, pos1, snp1
#            print c, pos2, snp2
#            print '----------------'

    print "same %d" % same
    print "mother diff %d" %motherdiff
    print "father diff %d" %fatherdiff
    print "child diff %d" %childdiff    
    print "phase diff %d" %phasediff
