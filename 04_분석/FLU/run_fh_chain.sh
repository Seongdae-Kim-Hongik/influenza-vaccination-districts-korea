#!/bin/bash
#  무거운 단계는 5분 부하 14 미만에서만 시작(lw)
S="${ST_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}/01_스크립트"; cd $S; T=FLU; L="$S/../04_분석/$T"
export ST_MV_N=127 ST_GRAPH=$L/graph_593.graph ST_THREADS=1 ST_VER=std3 ST_FBDZIP="${ST_FBDZIP:-}"
lw(){ while [ "$(sysctl -n vm.loadavg | awk '{print int($3)}')" -ge 14 ]; do sleep 30; done; }
ST_FH_ONLY=S2,S6 ./rslot.sh Rscript 81fh_flu.R $T >> $L/run_fh_main.log 2>&1
ST_FH_ONLY=S7 ST_MV_STOP=1 Rscript 81fh_flu.R $T >> $L/run_fh_main.log 2>&1
lw; ST_FH_ONLY=S4,S5,S8,S9 ./rslot.sh Rscript 81fh_flu.R $T > $L/run_fh_sens.log 2>&1 &
sleep 20; lw; ST_FH_ONLY=S10,S11,S12,S13,S14,S16,S18 ./rslot.sh Rscript 81fh_flu.R $T > $L/run_fh_sens2.log 2>&1 &
sleep 20; lw; ST_FH_ONLY=S21,S22,S23,S24 ./rslot.sh Rscript 81fh_flu.R $T > $L/run_fh_sens3.log 2>&1 &
for R in 1-32 33-64 65-96 97-127; do sleep 20; lw; ST_FH_ONLY=S7 ST_MV_RANGE=$R ./rslot.sh Rscript 81fh_flu.R $T > $L/run_fh_s7_$R.log 2>&1 & done
wait
ST_FH_ONLY=S7,S17,S44 ./rslot.sh Rscript 81fh_flu.R $T >> $L/run_fh_main.log 2>&1
echo "chain done $(date)" >> $L/run_fh_main.log; echo ALLDONE
