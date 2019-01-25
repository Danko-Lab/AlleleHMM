import sys

f = sys.argv[1]
#f_input=open(sys.argv[2], 'r')
f_input=sys.stdin
#outp=sys.stdout


with open (f+'_plus', 'w') as p_out:
    with open (f+'_minus', 'w') as m_out:
        with open (f+'_question', 'w') as q_out:
            for l in f_input:
                s = l.split('\t')[1]
                #print s
                if s == '+':
                    p_out.write(l)
                elif s == '-':
                    m_out.write(l)
                else:
                    q_out.write(l)
