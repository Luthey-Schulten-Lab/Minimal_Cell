#!/bin/bash
COUNTER=0

# While loop
while [ ${COUNTER} -lt 3 ]; do
COUNTER=$[${COUNTER}+1]
echo $COUNTER
done

echo Done!
