# 20220731

def diff_tad_bed(diff_tad_file, out_tad_file):
    f = open(out_tad_file, "w")
    with open(diff_tad_file, "r") as f1:
        diff_number = 1
        for line in f1:
            if line.startswith("#"):
                pass
            else:
                line_list = line.strip().split("\t")
                f.write(line_list[0] + "\t" + line_list[1] + "\t" + line_list[2] + "\t" + f"diff_tad_{diff_number}" + "\n")
                diff_number += 1

    f.close()
    return

def test_tad_bed(test_tad_file, out_tad_file):
    f = open(out_tad_file, "w")
    with open(test_tad_file, "r") as f1:
        diff_number = 1
        for line in f1:
            if line.startswith("#"):
                pass
            else:
                line_list = line.strip().split("\t")
                f.write(line_list[0] + "\t" + line_list[1] + "\t" + line_list[2] + "\t" + f"test_tad_{diff_number}" + "\n")
                diff_number += 1

    f.close()
    return

def control_tad_bed(control_tad_file, out_tad_file):
    f = open(out_tad_file, "w")
    with open(control_tad_file, "r") as f1:
        diff_number = 1
        for line in f1:
            if line.startswith("#"):
                pass
            else:
                line_list = line.strip().split("\t")
                f.write(line_list[0] + "\t" + line_list[1] + "\t" + line_list[2] + "\t" + f"control_tad_{diff_number}" + "\n")
                diff_number += 1

    f.close()
    return

import subprocess as subp

def Subp_call(cmd):
    print(cmd)
    subp.check_call(cmd, shell=True)
    return


# out_dir = ./path/path
# coefficient(0,1)
def main(diff_tad_file, test_tad_file, control_tad_file, out_dir, coefficient=0.1):
    
    # tad -> bed
    diff_tad_bed_file = f"{out_dir}/diff_tad.bed"
    test_tad_bed_file = f"{out_dir}/test_tad.bed"
    control_tad_bed_file = f"{out_dir}/control_tad.bed"
    diff_tad_bed(diff_tad_file, diff_tad_bed_file)
    test_tad_bed(test_tad_file, test_tad_bed_file)
    control_tad_bed(control_tad_file, control_tad_bed_file)
    
    # bedtools intersect
    bedtools_diff = f"{out_dir}/differential_tad_control_tad_intersect.tsv"
    Subp_call(f"bedtools intersect -a {diff_tad_bed_file} -b {control_tad_bed_file} > {bedtools_diff}")
    
    bedtools_countrol = f"{out_dir}/control_tad_differential_tad_intersect.tsv"
    Subp_call(f"bedtools intersect -a {control_tad_bed_file} -b {diff_tad_bed_file} > {bedtools_countrol}")
    
    # overlap and location
    dict_control = {}
    control_location = {}
    with open(bedtools_countrol, "r") as f:
        for line in f:
            line_list = line.strip().split("\t")
            diff_tad = line_list[7]
            control_tad = line_list[3]
            start = eval(line_list[1])
            end = eval(line_list[2])
            chrn = line_list[0]
            location_list = []
            location_list.append(start)
            location_list.append(end)
            location_list.append(chrn)
            if control_tad in dict_control.keys():
                over_list = dict_control[control_tad]
            else:
                over_list = []
                control_location[control_tad] = location_list
            over_list.append(diff_tad)
            dict_control[control_tad] = over_list
    
    dict_diff = {}
    diff_location = {}
    with open(bedtools_diff, "r") as f:
        for line in f:
            line_list = line.strip().split("\t")
            start = eval(line_list[1])
            end = eval(line_list[2])
            chrn = line_list[0]
            location_list = []
            location_list.append(start)
            location_list.append(end)
            location_list.append(chrn)
            diff_tad = line_list[3]
            control_tad = line_list[7]
            if diff_tad in dict_diff.keys():
                over_list = dict_diff[diff_tad]
            else:
                over_list = []
                diff_location[diff_tad] = location_list
            over_list.append(control_tad)
            dict_diff[diff_tad] = over_list
    
    # separate
    separate_tad = {}
    for diff_tad, control_list in dict_diff.items():
        # c 
        c_list = []
        if len(control_list) >= 2:
            c1 = [control_location[control_list[0]][0], control_location[control_list[-1]][1]]
            c_list.append(c1)
        if len(control_list) >= 3:
            c2 = [control_location[control_list[1]][0], control_location[control_list[-1]][1]]
            c3 = [control_location[control_list[0]][0], control_location[control_list[-2]][1]]
            c_list.append(c2)
            c_list.append(c3)
        if len(control_list) >= 4:
            c4 = [control_location[control_list[1]][0], control_location[control_list[-2]][1]]
            c_list.append(c4)
        
        location_a = diff_location[diff_tad]
        a0 = location_a[0]
        a2 = location_a[1]
        a1 = (a0 + a2) / 2
        d = (a1 - a0)*(coefficient*2)
        c_n = 0 # typle_c
        for c in c_list:
            c_n = c_n +1
            c0 = c[0]
            c2 = c[1]
            c1 = (c0 + c2) / 2
            left_d = abs(c0-a0)
            right_d = abs(c2-a2)
            if (a0 < c1 < a2) and (c0 < a1 < c2):
                if (left_d <= d) and (right_d <= d):
                    list_tuple = []
                    list_tuple.append(c)
                    list_tuple.append(c_n)
                    separate_tad[diff_tad]=list_tuple
                    #break
    
    # fusion
    fusion_tad = {}
    for control_tad, diff_list in dict_control.items():
        # c 
        c_list = []
        if len(diff_list) >= 2:
            c1 = [diff_location[diff_list[0]][0], diff_location[diff_list[-1]][1]]
            c_list.append(c1)
        if len(diff_list) >= 3:
            c2 = [diff_location[diff_list[1]][0], diff_location[diff_list[-1]][1]]
            c3 = [diff_location[diff_list[0]][0], diff_location[diff_list[-2]][1]]
            c_list.append(c2)
            c_list.append(c3)
        if len(diff_list) >= 4:
            c4 = [diff_location[diff_list[1]][0], diff_location[diff_list[-2]][1]]
            c_list.append(c4)
        
        location_a = control_location[control_tad]
        a0 = location_a[0]
        a2 = location_a[1]
        a1 = (a0 + a2) / 2
        d = (a1 - a0)*(coefficient*2)
        c_n = 0 # type_c
        for c in c_list:
            c_n = c_n + 1
            c0 = c[0]
            c2 = c[1]
            c1 = (c0 + c2) / 2
            left_d = abs(c0-a0)
            right_d = abs(c2-a2)
            if (a0 < c1 < a2) and (c0 < a1 < c2):
                if (left_d <= d) and (right_d <= d):
                    list_tuple = []
                    list_tuple.append(c)
                    list_tuple.append(c_n)
                    fusion_tad[control_tad]=list_tuple
                    #break
    
    # shifted unfilter
    shift_tad = {}
    for diff_tad, control_list in dict_diff.items():
        location_a = diff_location[diff_tad]
        a0 = location_a[0]
        a2 = location_a[1]
        a1 = (a0 + a2) / 2
        d = (a1 - a0)*(coefficient * 2)
        for control in control_list:
            location_b = control_location[control]
            b0 = location_b[0]
            b2 = location_b[1]
            b1 = (b0 + b2) / 2
            left_d = abs(b0-a0)
            right_d = abs(b2-a2)
            if (a0 < b1 < a2) and (b0 < a1 < b2):
                if ((left_d <= d) and (right_d > d)):
                    c_n = 2 # right change
                    list_tuple = []
                    list_tuple.append(control)
                    list_tuple.append(c_n)
                    shift_tad[diff_tad] = list_tuple
                elif ((left_d > d) and (right_d <= d)):
                    c_n = 1 # left  change
                    list_tuple = []
                    list_tuple.append(control)
                    list_tuple.append(c_n)
                    shift_tad[diff_tad] = list_tuple
    
    # output tad type file
    boundary_list = []        
    ## separate
    separate_tad_file = f"{out_dir}/test_to_control_separate.csv"
    f = open(separate_tad_file, "w")
    for tad in separate_tad.keys():
        location_a = diff_location[tad]
        chrn = location_a[2]
        type_ = separate_tad[tad][-1]
        if type_ == 1:
            tad_list = dict_diff[tad]            
        elif type_ == 2:
            tad_list = dict_diff[tad][1:]
        elif type_ == 3:
            tad_list = dict_diff[tad][:-1]
        elif type_ == 4:
            tad_list = dict_diff[tad][1:-1]
        else:
            pass
        
        # boundary location
        tmp_list = []# 将所有边界信息存好
        for b in tad_list:
            location_b = control_location[b]
            f.write(tad + "," + location_a[2] + "," + str(location_a[0]) + "," + str(location_a[1]) + "," + b + "," + location_b[2] + "," + str(location_b[0]) + "," + str(location_b[1]) + "," + str(type_) + "\n")
            bou_list = []
            bou_list.append(location_b[2])
            bou_list.append(location_b[0])
            tmp_list.append(bou_list)
            bou_list = []
            bou_list.append(location_b[2])
            bou_list.append(location_b[1])
            tmp_list.append(bou_list)
        for bou_real in tmp_list[1:-1]:# 过滤开头和结尾
            boundary_list.append(bou_real)   
    f.close()
    ## fusion : test_to_control_fusion.csv
    fusion_tad_file = f"{out_dir}/test_to_control_fusion.csv"
    f = open(fusion_tad_file, "w")
    fusion_both = []
    for tad in fusion_tad.keys():
        location_b = control_location[tad]
        type_ = fusion_tad[tad][-1]
        if type_ == 1:
            tad_list = dict_control[tad]
        elif type_ == 2:
            tad_list = dict_control[tad][1:]
        elif type_ == 3:
            tad_list = dict_control[tad][:-1]
        elif type_ == 4:
            tad_list = dict_control[tad][1:-1]
        else:
            pass
        tmp_list = []
        for a in tad_list:
            fusion_both.append(a)
            location_a = diff_location[a]
            f.write(a + "," + location_a[2] + "," + str(location_a[0]) + "," + str(location_a[1]) + "," + tad + "," + location_b[2] + "," + str(location_b[0]) + "," + str(location_b[1]) + "," + str(type_) + "\n")
            bou_list = []
            bou_list.append(location_a[2])
            bou_list.append(location_a[0])
            tmp_list.append(bou_list)
            bou_list = []
            bou_list.append(location_a[2])
            bou_list.append(location_a[1])
            tmp_list.append(bou_list)
        for bou_real in tmp_list[1:-1]:# 过滤开头和结尾
            boundary_list.append(bou_real)
    f.close()
    
    ## shilt(test_to_control_shift.csv)
    shift_tad_file = f"{out_dir}/test_to_control_shift.csv"
    f = open(shift_tad_file, "w")
    for tad in shift_tad.keys():
        if tad in separate_tad.keys():
            pass
        elif tad in fusion_both:
            pass
        else:
            location_a = diff_location[tad]
            b = shift_tad[tad][0]
            location_b = control_location[b]
            type_ = shift_tad[tad][1]
            f.write(tad + "," + location_a[2] + "," + str(location_a[0]) + "," + str(location_a[1]) + "," + b + "," + location_b[2] + "," + str(location_b[0]) + "," + str(location_b[1]) + "," + str(type_) + "\n")
            # boundary location
            bou_list = []
            bou_list.append(location_a[2])
            if type_ == 1:
                bou_list.append(max(location_a[0], location_b[0]))
            else:
                bou_list.append(min(location_a[1], location_b[1]))
            boundary_list.append(bou_list)
    f.close()
    
    ## boundary file: tad_type_boundary.csv
    boundary_file = f"{out_dir}/tad_type_boundary.csv"
    with open(boundary_file, "w") as f:
        f.write("chrn" + "," + "change_position" + "\n")
        for bou_list in boundary_list:
            f.write(bou_list[0] + "," + str(bou_list[1]) + "\n")
    
    return