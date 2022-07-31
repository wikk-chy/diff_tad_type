# diff_tad_type
![avatar](./tad_type_all.jpg)  
Analysis the type of diff_tad which from HicExplorer tools HicDifferentialTAD  
more infomation see the PPT file(diff_tad_type.pptx)  

## input file  
test_differential_rejected.diff_tad(diff tad file)  
test_tad_domains.bed(test tad file)  
control_tad_domains.bed(control tad file)

## output file
differential_tad.bed    
test_tad.bed    
control_tad.bed  
differential_tad_control_tad_intersect.tsv  
control_tad_differential_tad_intersect.tsv  
test_to_control_shift.csv  
test_to_control_separate.csv  
test_to_control_fusion.csv  
tad_type_boundary.csv  


## useage
```python
python diff_tad_type.py diff_tad_file test_tad_file control_tad_file out_dir percent  
```

## more
percent = percent_conserve = percent_change = (0,1), defaults is 0.1
