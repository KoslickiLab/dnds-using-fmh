#!/bin/bash

# run kaks_calculator on axt file
nohup KaKs_Calculator -i negative_selection_queries_10002_0.01.fna.axt -o negative_selection_queries_10002_0.01.axt.kaks -m GNG -m GY -m LPB -m LWL -m NG -m YN > .kaks_calc.log 2>&1 &


