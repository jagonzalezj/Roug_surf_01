#! /bin/bash
# Here are some examples on how to use pyrough.py

python pyrough.py atoms intercept Cu1.lmp wire
python pyrough.py atoms intercept Cu1.lmp slab
python pyrough.py atoms intercept Cu1.lmp sphere
python pyrough.py atoms intercept Cu1.lmp wire RS b 1.8 C 10 N 30 M 30 graph
python pyrough.py atoms intercept Cu1.lmp wire graph
python pyrough.py atoms intercept Cu1.lmp slab RS b 2.8
