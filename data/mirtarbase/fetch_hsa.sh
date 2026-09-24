#!/bin/bash
# Robust sequential fetch of hsa_MTI.csv (server lacks range resume; 429s on rapid requests).
cd ~/mega27/18/data/mirtarbase
U='https://awi.cuhk.edu.cn/miRTarBase/downloads/files/10.0/hsa_MTI.csv'
for i in 1 2 3 4 5 6; do
  rm -f hsa_MTI.part hdr.txt
  curl -s -L -D hdr.txt -o hsa_MTI.part "$U"; rc=$?
  cl=$(grep -i '^content-length' hdr.txt | tail -1 | tr -dc '0-9')
  sz=$(stat -c%s hsa_MTI.part 2>/dev/null)
  echo "attempt $i rc=$rc size=$sz content_length=$cl $(date +%T)" >> dl.log
  if [ "$rc" = "0" ] && { [ -z "$cl" ] || [ "$sz" = "$cl" ]; }; then
    mv hsa_MTI.part hsa_MTI.csv
    cd ~/mega27/18
    setsid nohup python3 scripts/mirtarbase_validate.py > results/mirtarbase_validate.log 2>&1 &
    echo "analysis launched pid $! $(date +%T)" >> data/mirtarbase/dl.log
    exit 0
  fi
  sleep 90
done
echo "hsa fetch FAILED after 6 attempts" >> dl.log
