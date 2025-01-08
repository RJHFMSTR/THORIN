
ibd_p = "ibd.txt"
phase_p = "unphased_vs_phased.txt"

#------------------------------------------------------------------------------#
# Import IBD states
#------------------------------------------------------------------------------#

ibd_state = {'A': [], 'B': [], 'C':[], 'D':[]}
with open(ibd_p, 'r') as ibd_s:
    ibd_s.readline()
    for line in ibd_s:
        splt_line = line.rstrip().split('\t')
        ibd_state[splt_line[3]].append((int(splt_line[1]), int(splt_line[2])))

#------------------------------------------------------------------------------#
# Read  
#------------------------------------------------------------------------------#

with open(phase_p, 'r') as phase_s:
    for line in phase_s:
        pos, phasing = line.rstrip().split()
        pos = int(pos)
        check_pos_in = False
        if phasing == "|":
            for state in ('A', 'B'):
                for start, end in ibd_state[state]:
                    if pos <= end and pos >= start:
                        check_pos_in = True
                        break
        else:
            assert phasing == "/"
            for state in ('C', 'D'):
                for start, end in ibd_state[state]:
                    if pos <= end and pos >= start:
                        check_pos_in = True
                        break
        
        if not check_pos_in: raise Exception("Position", pos, "is incorrect")

