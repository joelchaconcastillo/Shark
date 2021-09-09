
PATH1=../SMS-EMOA/

DI=0.4;
for px in 0.2 0.4 0.6 0.8 1.0
#for px in 0.2
do

  for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF1 UF2 UF3 UF4 UF5 UF6 UF7 IMB1 IMB2 IMB3 IMB7 IMB8 IMB9
  do
     for seed in {1..35}
     do
     tail -100 $PATH1/POF_${instance}_nobj_2_*_seed_${seed}_px_${px}*_2500000 | cut -f1,2 -d' '  > POF/${instance}_2_${seed}_${px}
     done
  done
  for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF8 UF9 UF10 IMB4 IMB5 IMB6 IMB10
  do
     for seed in {1..35}
     do
     tail -100 $PATH1/POF_${instance}_nobj_3_*_seed_${seed}_px_${px}*_2500000 | cut -f1,2,3 -d' '  > POF/${instance}_3_${seed}_${px}
     done
  done
done
