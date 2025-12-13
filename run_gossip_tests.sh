#!/bin/bash
# Run gossip simulation tests with varying gossip ratios and propagation models
# Usage: ./scratch/run_gossip_tests.sh

set -e  # Exit on error

SIM_TIME=360
GOSSIP_RATIOS=(0.1 0.2 0.3 0.4 0.5 0.6)

echo "=============================================="
echo "Gossip Simulation Test Suite"
echo "SimTime: ${SIM_TIME}s"
echo "Gossip Ratios: ${GOSSIP_RATIOS[*]}"
echo "=============================================="

# Build first
echo ""
echo "Building simulation..."
./ns3 build scratch/simulationWifi

# Track total runs
TOTAL_RUNS=$(( ${#GOSSIP_RATIOS[@]} * 3 ))  # 3 propagation configs
CURRENT_RUN=0

# Function to run a single test
run_test() {
    local exp=$1
    local gr=$2
    local name=$3
    
    CURRENT_RUN=$((CURRENT_RUN + 1))
    echo ""
    echo "----------------------------------------------"
    echo "Run ${CURRENT_RUN}/${TOTAL_RUNS}: ${name}"
    echo "  PathLossExp=${exp}, GossipRatio=${gr}"
    echo "----------------------------------------------"
    
    ./ns3 run "scratch/simulationWifi \
        --simTime=${SIM_TIME} \
        --gossipRatio=${gr} \
        --pathLossExp=${exp} \
        --runName=${name}"
    
    echo "Completed: ${name}"
}

# ===== Test 1: Default YansWifiChannelHelper (LogDistance exp=3) =====
echo ""
echo "=============================================="
echo "Testing: Default YansWifiChannelHelper"
echo "=============================================="
for gr in "${GOSSIP_RATIOS[@]}"; do
    gr_str=$(echo $gr | tr '.' '_')
    run_test 0 $gr "default_gr${gr_str}_t${SIM_TIME}s"
done

# ===== Test 2: Log-distance exponent 3.3 =====
echo ""
echo "=============================================="
echo "Testing: LogDistance Exponent 3.3 (Light Forest)"
echo "=============================================="
for gr in "${GOSSIP_RATIOS[@]}"; do
    gr_str=$(echo $gr | tr '.' '_')
    run_test 3.3 $gr "exp3_3_gr${gr_str}_t${SIM_TIME}s"
done

# ===== Test 3: Log-distance exponent 3.4 =====
echo ""
echo "=============================================="
echo "Testing: LogDistance Exponent 3.4 (Dense Forest)"
echo "=============================================="
for gr in "${GOSSIP_RATIOS[@]}"; do
    gr_str=$(echo $gr | tr '.' '_')
    run_test 3.4 $gr "exp3_4_gr${gr_str}_t${SIM_TIME}s"
done

echo ""
echo "=============================================="
echo "All ${TOTAL_RUNS} tests completed!"
echo "Results saved in simWifi/ directory:"
ls -d simWifi/*/ 2>/dev/null | head -20
echo "=============================================="
