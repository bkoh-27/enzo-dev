#!/bin/bash
PASS=0
FAIL=0
check_log() {
    local label=$1
    local logfile=$2
    local step1 step2
    step1=$(grep '\[BHSEED\]' "$logfile" | sed -n '1p')
    step2=$(grep '\[BHSEED\]' "$logfile" | sed -n '2p')
    if echo "$step1" | grep -q 'created=1'; then echo "$label step 1 created=1: PASS"; PASS=$((PASS+1)); else echo "$label step 1 created=1: FAIL — got: $step1"; FAIL=$((FAIL+1)); fi
    if echo "$step1" | grep -q 'ncand_global=2'; then echo "$label step 1 ncand_global=2: PASS"; PASS=$((PASS+1)); else echo "$label step 1 ncand_global=2: FAIL — got: $step1"; FAIL=$((FAIL+1)); fi
    if echo "$step2" | grep -q 'dist_blocked=2'; then echo "$label step 2 dist_blocked=2: PASS"; PASS=$((PASS+1)); else echo "$label step 2 dist_blocked=2: FAIL — got: $step2"; FAIL=$((FAIL+1)); fi
    if echo "$step2" | grep -q 'created=0'; then echo "$label step 2 created=0: PASS"; PASS=$((PASS+1)); else echo "$label step 2 created=0: FAIL — got: $step2"; FAIL=$((FAIL+1)); fi
}
echo "=============================="
echo "  TS3_wrap Assertion Results"
echo "=============================="
check_log "np1" "ts3wrap_np1.log"
check_log "np4" "ts3wrap_np4.log"
echo "------------------------------"
echo "Total: $PASS passed, $FAIL failed"
if [ $FAIL -eq 0 ]; then echo "TS3_WRAP: PASS"; exit 0; else echo "TS3_WRAP: FAIL"; exit 1; fi
