#!/bin/bash
MIN_NUM=0
MAX_NUM=999
for ((NUM=$MIN_NUM; NUM < $MAX_NUM; NUM++)); do
if [ ${NUM} -lt 10 ]; then
    JOBNUM=000${NUM}
elif [ ${NUM} -gt 9 ] && [ ${NUM} -lt 100 ]; then
    JOBNUM=00${NUM}
elif [ ${NUM} -gt 99 ] && [ ${NUM} -lt 1000 ]; then
    JOBNUM=0${NUM}
else
    JOBNUM=${NUM}
fi

echo "Start ${JOBNUM}"
mkdir ${JOBNUM}
echo "Finish ${JOBNUM}"

if [ $(( $NUM % 20 )) == 0 ]; then
        echo 'SLEEP'
       wait
    fi
done
