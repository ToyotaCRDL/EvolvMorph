#!/bin/bash


touch $$
while [ ! -f USPEX_IS_DONE ]&&[ ! -f ABORT ]; do
  date >> log
  date >> err
  USPEX -r -o 1>> log 2>>err
sleep 60
done

if [ -f '$$' ]; then
  rm -f '$$'
fi
