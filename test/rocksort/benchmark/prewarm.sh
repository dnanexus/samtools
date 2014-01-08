#!/bin/bash -e

dd if=/dev/zero of=/mnt/00FILLER bs=16M &
dd_pid=$!
while ps -p $dd_pid --no-heading ; do
        sleep 60s
        kill -USR1 $dd_pid || true
done
rm /mnt/00FILLER

touch /mnt/.PREWARMED
