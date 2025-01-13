
ibd_p = "ibd.txt"
phase_p = "unphased_vs_phased.3cM.txt"
cM_size = 3

#------------------------------------------------------------------------------#
# Import IBD states
#------------------------------------------------------------------------------#

ibd_state = {'A': [], 'B': [], 'C':[], 'D':[]}
with open(ibd_p, 'r') as ibd_s:
    ibd_s.readline()
    for line in ibd_s:
        splt_line = line.rstrip().split('\t')
        letter = splt_line[3]
        ibd_state[letter].append((int(splt_line[1]), int(splt_line[2]), float(splt_line[4])))

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
                for start, end, cM in ibd_state[state]:
                    if pos <= end and pos >= start and cM >= cM_size:
                        check_pos_in = True
                        break
        else:
            assert phasing == "/"
            # if we have "/" it means either it was in c or d or it was a A, B with a size < 3
            # first check if it is c or d
            for state in ('C', 'D'):
                for start, end, cM in ibd_state[state]:
                    if pos <= end and pos >= start:
                        check_pos_in = True
                        break
            # if not, check if it is a small segment
            if not check_pos_in:
                for state in ('A', 'B'):
                    for start, end, cM in ibd_state[state]:
                        if pos <= end and pos >= start:
                            check_pos_in = True
                            break
        
        if not check_pos_in: raise Exception("Position", pos, "is incorrect")

