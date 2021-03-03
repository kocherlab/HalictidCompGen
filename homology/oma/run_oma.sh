#!/bin/bash

sbatch --array=1-100 -N1 --qos=1wk --time=6-23:00<<EOF
#!/bin/sh 
export NR_PROCESSES=100
/Genomics/kocherlab/berubin/local/src/OMA.2.3.0/oma/OMA/bin/OMA
EOF
