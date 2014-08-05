#!/bin/bash

#joblist=("fastsimcoal_sim.sh" \
#"submit_ms.sh" \
#"macs_retain0.sh" "macs_retain500.sh" "macs_retain1000.sh" \
#"macs_retain3000.sh" "macs_retain5000.sh" "macs_retain7000.sh" "macs_retain10000.sh" \
#"macs_retain30000.sh" "macs_retain50000.sh" "macs_retain70000.sh" "macs_retain100000.sh" \
#"submit_scrm0.sh" "submit_scrm500.sh" "submit_scrm1000.sh" \
#"submit_scrm3000.sh" "submit_scrm5000.sh" "submit_scrm7000.sh" "submit_scrm10000.sh" \
#"submit_scrm30000.sh" "submit_scrm50000.sh" "submit_scrm70000.sh" "submit_scrm100000.sh" \
#)

#joblist=("submit_scrm0.sh" "submit_scrm500.sh" "submit_scrm1000.sh" \ 
#"submit_scrm3000.sh" "submit_scrm5000.sh" "submit_scrm7000.sh" "submit_scrm10000.sh" \
#"submit_scrm30000.sh" "submit_scrm50000.sh" "submit_scrm70000.sh" "submit_scrm100000.sh" \
#)

joblist=("submit_scrm300000.sh" "submit_scrm500000.sh" "submit_scrm700000.sh" "submit_scrm1000000.sh" \
"macs_retain300000.sh" "macs_retain500000.sh" "macs_retain700000.sh" "macs_retain1000000.sh" \
)

for jobi in $(seq 0 1 23)
    do
    echo ${joblist[${jobi}]}
    qsub ${joblist[${jobi}]}
    done
