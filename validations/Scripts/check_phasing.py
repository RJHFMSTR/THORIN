
ibd_p = "ibd.txt"
phase_p = "changed_phase.txt"

#------------------------------------------------------------------------------#
# Import IBD states
#------------------------------------------------------------------------------#

ibd_b_pos = []
with open(ibd_p, 'r') as ibd_s:
    ibd_s.readline()
    for line in ibd_s:
        splt_line = line.rstrip().split('\t')
        if splt_line[3] == 'B': ibd_b_pos.append((int(splt_line[1]), int(splt_line[2])))

#------------------------------------------------------------------------------#
# Read  
#------------------------------------------------------------------------------#

with open(phase_p, 'r') as phase_s:
    for line in phase_s:
        pos = int(line.rstrip())
        check_pos_in = False
        for start, end in ibd_b_pos:
            if pos <= end and pos >= start: check_pos_in = True
        if not check_pos_in:
            raise Exception("Position", pos, "should not have changed phase")
