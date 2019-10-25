python INITIAL_STEP.py
make
./execslice
tot=$(ls -d SLICED_* | wc -l)
python Python_main.py $tot
python Final_step.py $tot
python ANALYSIS.py $tot
