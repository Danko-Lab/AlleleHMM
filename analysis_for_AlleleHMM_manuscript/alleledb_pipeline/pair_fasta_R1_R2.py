# python pair_fasta_R1_R2.py DIPLO_CGATA_R1_filtered.fasta

import sys

#######################################################
#
# Change the following settings to suit your needs

file_head= sys.argv[1]

input_forward_filename = file_head+"_R1.fasta"
input_reverse_filename = file_head+"_R2.fasta"

#output_pairs_filename = "out_interleaved_pairs.fastq"
output_orphan_filename = file_head+"_singleton.fasta"
output_paired_forward_filename = file_head+"_paired_R1.fasta"
output_paired_reverse_filename = file_head+"_paired_R2.fasta"

print output_paired_forward_filename


#######################################################
print "Scanning reverse file to build list of names..."    
reverse_ids = set()
paired_ids = set()
with open(input_reverse_filename) as rf:
    while True:
        header = rf.readline().strip().split(' ')[0]
        if header=="": break
        assert(header[0]=='>')
        seq=rf.readline().rstrip()
        reverse_ids.add(header)


print "Processing forward file..."
forward_handle = open(output_paired_forward_filename, "w")
orphan_handle = open(output_orphan_filename, "w")
with open(input_forward_filename) as ff:
    while True:
        header = ff.readline().strip().split(' ')[0]
        if header=="": break
        assert(header[0]=='>')
        seq=ff.readline().rstrip()
        if header in reverse_ids:
            #Paired
            paired_ids.add(header)
            reverse_ids.remove(header) #frees a little memory
            forward_handle.write("%s_1\n%s\n" % (header, seq))
        else:
            #Orphan
            orphan_handle.write("%s_1\n%s\n" % (header, seq))
forward_handle.close()

del reverse_ids #frees memory, although we won't need more now
print "Processing reverse file..."
reverse_handle = open(output_paired_reverse_filename, "w")
with open(input_reverse_filename) as rf:
    while True:
        header = rf.readline().strip().split(' ')[0]
        if header=="": break
        assert(header[0]=='>')
        seq=rf.readline().rstrip()
        if header in paired_ids:
            #Paired
            reverse_handle.write("%s_2\n%s\n" % (header, seq))
        else:
            #Orphan
            orphan_handle.write("%s_2\n%s\n" % (header, seq))
orphan_handle.close()
reverse_handle.close()
print "Done"