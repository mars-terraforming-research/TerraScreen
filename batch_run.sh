#!/usr/bin/env bash


filelist=(
"QEXT_96IR_84_VIS_Al_8um_60um"
"QEXT_96IR_84_VIS_fosterite_100nm"
"QEXT_96IR_84_VIS_nanoring_9um7um"
"QEXT_96IR_84_VIS_coreshell_9p1um9um"
"QEXT_96IR_84_VIS_fosterite_1um"
"QEXT_96IR_84_VIS_dust1.5"
"QEXT_96IR_84_VIS_graph_1000nm"
"QEXT_96IR_84_VIS_nanoribbon_9um1um0p1um"
"QEXT_96IR_84_VIS_nanoring_9um8p8um"
"QEXT_96IR_84_VIS_SiO2_Rod_5um_1um"
"QEXT_96IR_84_VIS_SiO2_Rod_7um_1um"
"QEXT_96IR_84_VIS_SiO2_Rod_9um_1um"
"QEXT_96IR_84_VIS_Al2O3_0p5um_poro_90pct"
"QEXT_96IR_84_VIS_SiO2_0p5um_poro_90pct"
"QEXT_96IR_84_VIS_nanoribbon_9um1um0p1um"
"QEXT_96IR_84_VIS_graph_mix"
#"QEXT_96IR_84_VIS_graph_mix_Qext_cst_w0"
#"QEXT_96IR_84_VIS_graph_mix_Qext_cst_w1_g1"
#"QEXT_96IR_84_VIS_graph_mix_Qext_cst_w1_g0"
#"QEXT_96IR_84_VIS_graph_mix_Qext_cst_w1_gm1"

)

for f in "${filelist[@]}"; do
    ./TerraScreen data/$f      
done
