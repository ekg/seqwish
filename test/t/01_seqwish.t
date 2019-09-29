#!/usr/bin/env bash

BASH_TAP_ROOT=bash-tap
. bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for seqwish

plan tests 29

is $(seqwish -h 2>&1 | grep "seqwish: a variation graph inducer" | wc -l) 1 "seqwish prints its help"

is $( seqwish -s HLA/A-3105.fa.gz -p HLA/A-3105.paf.gz -b HLA/A-3105.fa.gz.work -g HLA/A-3105.fa.gz.gfa && md5sum HLA/A-3105.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/A-3105.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for A-3105"
is $( seqwish -s HLA/B-3106.fa.gz -p HLA/B-3106.paf.gz -b HLA/B-3106.fa.gz.work -g HLA/B-3106.fa.gz.gfa && md5sum HLA/B-3106.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/B-3106.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for B-3106"
is $( seqwish -s HLA/C-3107.fa.gz -p HLA/C-3107.paf.gz -b HLA/C-3107.fa.gz.work -g HLA/C-3107.fa.gz.gfa && md5sum HLA/C-3107.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/C-3107.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for C-3107"
is $( seqwish -s HLA/DMA-3108.fa.gz -p HLA/DMA-3108.paf.gz -b HLA/DMA-3108.fa.gz.work -g HLA/DMA-3108.fa.gz.gfa && md5sum HLA/DMA-3108.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DMA-3108.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DMA-3108"
is $( seqwish -s HLA/DMB-3109.fa.gz -p HLA/DMB-3109.paf.gz -b HLA/DMB-3109.fa.gz.work -g HLA/DMB-3109.fa.gz.gfa && md5sum HLA/DMB-3109.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DMB-3109.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DMB-3109"
is $( seqwish -s HLA/DOA-3111.fa.gz -p HLA/DOA-3111.paf.gz -b HLA/DOA-3111.fa.gz.work -g HLA/DOA-3111.fa.gz.gfa && md5sum HLA/DOA-3111.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DOA-3111.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DOA-3111"
is $( seqwish -s HLA/DOB-3112.fa.gz -p HLA/DOB-3112.paf.gz -b HLA/DOB-3112.fa.gz.work -g HLA/DOB-3112.fa.gz.gfa && md5sum HLA/DOB-3112.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DOB-3112.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DOB-3112"
is $( seqwish -s HLA/DPA1-3113.fa.gz -p HLA/DPA1-3113.paf.gz -b HLA/DPA1-3113.fa.gz.work -g HLA/DPA1-3113.fa.gz.gfa && md5sum HLA/DPA1-3113.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DPA1-3113.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DPA1-3113"
is $( seqwish -s HLA/DPB1-3115.fa.gz -p HLA/DPB1-3115.paf.gz -b HLA/DPB1-3115.fa.gz.work -g HLA/DPB1-3115.fa.gz.gfa && md5sum HLA/DPB1-3115.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DPB1-3115.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DPB1-3115"
is $( seqwish -s HLA/DQA1-3117.fa.gz -p HLA/DQA1-3117.paf.gz -b HLA/DQA1-3117.fa.gz.work -g HLA/DQA1-3117.fa.gz.gfa && md5sum HLA/DQA1-3117.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DQA1-3117.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DQA1-3117"
is $( seqwish -s HLA/DQB1-3119.fa.gz -p HLA/DQB1-3119.paf.gz -b HLA/DQB1-3119.fa.gz.work -g HLA/DQB1-3119.fa.gz.gfa && md5sum HLA/DQB1-3119.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DQB1-3119.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DQB1-3119"
is $( seqwish -s HLA/DRA-3122.fa.gz -p HLA/DRA-3122.paf.gz -b HLA/DRA-3122.fa.gz.work -g HLA/DRA-3122.fa.gz.gfa && md5sum HLA/DRA-3122.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DRA-3122.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DRA-3122"
is $( seqwish -s HLA/DRB1-3123.fa.gz -p HLA/DRB1-3123.paf.gz -b HLA/DRB1-3123.fa.gz.work -g HLA/DRB1-3123.fa.gz.gfa && md5sum HLA/DRB1-3123.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DRB1-3123.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DRB1-3123"
is $( seqwish -s HLA/DRB3-3125.fa.gz -p HLA/DRB3-3125.paf.gz -b HLA/DRB3-3125.fa.gz.work -g HLA/DRB3-3125.fa.gz.gfa && md5sum HLA/DRB3-3125.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DRB3-3125.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DRB3-3125"
is $( seqwish -s HLA/DRB4-3126.fa.gz -p HLA/DRB4-3126.paf.gz -b HLA/DRB4-3126.fa.gz.work -g HLA/DRB4-3126.fa.gz.gfa && md5sum HLA/DRB4-3126.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DRB4-3126.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DRB4-3126"
is $( seqwish -s HLA/DRB5-3127.fa.gz -p HLA/DRB5-3127.paf.gz -b HLA/DRB5-3127.fa.gz.work -g HLA/DRB5-3127.fa.gz.gfa && md5sum HLA/DRB5-3127.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/DRB5-3127.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for DRB5-3127"
is $( seqwish -s HLA/E-3133.fa.gz -p HLA/E-3133.paf.gz -b HLA/E-3133.fa.gz.work -g HLA/E-3133.fa.gz.gfa && md5sum HLA/E-3133.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/E-3133.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for E-3133"
is $( seqwish -s HLA/F-3134.fa.gz -p HLA/F-3134.paf.gz -b HLA/F-3134.fa.gz.work -g HLA/F-3134.fa.gz.gfa && md5sum HLA/F-3134.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/F-3134.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for F-3134"
is $( seqwish -s HLA/G-3135.fa.gz -p HLA/G-3135.paf.gz -b HLA/G-3135.fa.gz.work -g HLA/G-3135.fa.gz.gfa && md5sum HLA/G-3135.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/G-3135.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for G-3135"
is $( seqwish -s HLA/H-3136.fa.gz -p HLA/H-3136.paf.gz -b HLA/H-3136.fa.gz.work -g HLA/H-3136.fa.gz.gfa && md5sum HLA/H-3136.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/H-3136.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for H-3136"
is $( seqwish -s HLA/J-3137.fa.gz -p HLA/J-3137.paf.gz -b HLA/J-3137.fa.gz.work -g HLA/J-3137.fa.gz.gfa && md5sum HLA/J-3137.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/J-3137.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for J-3137"
is $( seqwish -s HLA/K-3138.fa.gz -p HLA/K-3138.paf.gz -b HLA/K-3138.fa.gz.work -g HLA/K-3138.fa.gz.gfa && md5sum HLA/K-3138.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/K-3138.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for K-3138"
is $( seqwish -s HLA/L-3139.fa.gz -p HLA/L-3139.paf.gz -b HLA/L-3139.fa.gz.work -g HLA/L-3139.fa.gz.gfa && md5sum HLA/L-3139.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/L-3139.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for L-3139"
is $( seqwish -s HLA/MICA-100507436.fa.gz -p HLA/MICA-100507436.paf.gz -b HLA/MICA-100507436.fa.gz.work -g HLA/MICA-100507436.fa.gz.gfa && md5sum HLA/MICA-100507436.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/MICA-100507436.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for MICA-100507436"
is $( seqwish -s HLA/MICB-4277.fa.gz -p HLA/MICB-4277.paf.gz -b HLA/MICB-4277.fa.gz.work -g HLA/MICB-4277.fa.gz.gfa && md5sum HLA/MICB-4277.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/MICB-4277.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for MICB-4277"
is $( seqwish -s HLA/TAP1-6890.fa.gz -p HLA/TAP1-6890.paf.gz -b HLA/TAP1-6890.fa.gz.work -g HLA/TAP1-6890.fa.gz.gfa && md5sum HLA/TAP1-6890.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/TAP1-6890.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for TAP1-6890"
is $( seqwish -s HLA/TAP2-6891.fa.gz -p HLA/TAP2-6891.paf.gz -b HLA/TAP2-6891.fa.gz.work -g HLA/TAP2-6891.fa.gz.gfa && md5sum HLA/TAP2-6891.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/TAP2-6891.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for TAP2-6891"
is $( seqwish -s HLA/V-352962.fa.gz -p HLA/V-352962.paf.gz -b HLA/V-352962.fa.gz.work -g HLA/V-352962.fa.gz.gfa && md5sum HLA/V-352962.fa.gz.gfa | cut -f 1 -d\ ) $( cat HLA/V-352962.fa.gz.gfa.md5 ) "seqwish correctly builds the graph for V-352962"

rm -f HLA/*gfa
