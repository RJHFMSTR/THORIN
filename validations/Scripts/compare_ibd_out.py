import sys

#------------------------------------------------------------------------------#
# Params
#------------------------------------------------------------------------------#
ibd_v1_p = sys.argv[1]
ibd_v2_p = sys.argv[2]

#------------------------------------------------------------------------------#
# Compare files
#------------------------------------------------------------------------------#

ibd_v1_s= open(ibd_v1_p, 'r')
ibd_v2_s= open(ibd_v2_p, 'r')
assert ibd_v1_s.readline() == ibd_v2_s.readline() # check header

for line1, line2 in zip(ibd_v1_s, ibd_v2_s):
    splt_line1 = line1.split()
    splt_line2 = line2.split()
    # check CHR, target and UPD
    assert splt_line1[0] == splt_line2[0] and splt_line1[5] == splt_line2[5] and splt_line1[6] == splt_line2[6]

    # check length in cM
    assert (float(splt_line1[4]) - float(splt_line2[4]) < 0.1)
    
    # compare positions
    pos_1 = list(map(int, splt_line1[1:3]))
    pos_2 = list(map(int, splt_line2[1:3]))
    if not (abs(pos_1[0] - pos_2[0]) < 1000 and abs(pos_1[1] - pos_2[1]) < 1000):
        print("------------------------------------------------------------------------------")
        print("This lines have a difference > 1000 bp in one of their position")
        print(line1)
        print(line2)

ibd_v1_s.close()
ibd_v2_s.close()