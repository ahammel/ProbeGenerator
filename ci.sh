#!/bin/bash

ROOT_DIR="$HOME/.ci"
PG_WORKING_REPO="$HOME/repos/probing-pipeline/trunk/ProbeGenerator"

_python=/gsc/software/linux-x86_64/python-3.2.2/bin/python3

FUSION_PROBES="$HOME/scratch/probes/actionable_fusions/statement.txt"
SNP_PROBES="$HOME/scratch/probes/gene_snp_probes/actionable_mutation_statements.txt.csv"

EMAIL_ADDRESS="ahammel@bcgsc.ca"

send_email () {
    (
        echo "To: ahammel@bcgsc.ca"
        echo "Subject: probe-generator CI TestResults"
        echo "Content-Type: text/html; charset=utf-8"
        echo
        echo "$1"
        echo
    ) | /usr/sbin/sendmail -t
}

run_pg () {

    cd "$ROOT_DIR"/ProbeGenerator
    ssh xhost08 "export PYTHONPATH=$ROOT_DIR/usr/lib/python3.2/site-packages:$PYTHONPATH &&
                 $_python -m probe_generator \
                 -a /home/ahammel/scratch/probes/refseq_genes.txt \
                 -a /home/ahammel/scratch/probes/ucsc_genes.txt   \
                 -g /genesis/scratch/transabyss/trans-ABySS/annotations/hg19/201304/genome.fa \
                 -s $1"
}

mkdir -p "$ROOT_DIR"/usr

cd "$PG_WORKING_REPO"
revision=$(git rev-parse HEAD)

cd "$ROOT_DIR"
rm -rf ProbeGenerator
rm -f new.fa
rm -f truth.fa
git clone "$PG_WORKING_REPO"
mkdir -p usr
cat "$FUSION_PROBES" "$SNP_PROBES" > "$ROOT_DIR"/probes.txt

cd "$ROOT_DIR"/ProbeGenerator

git checkout "$revision"
$_python setup.py install --prefix="$ROOT_DIR"/usr
run_pg "$ROOT_DIR"/probes.txt > "$ROOT_DIR"/new.fa 2> "$ROOT_DIR/new_errors.txt"

git checkout master
$_python setup.py install --prefix="$ROOT_DIR"/usr
run_pg "$ROOT_DIR"/probes.txt > "$ROOT_DIR"/old.fa 2> "$ROOT_DIR/truth_errors.txt"

difference=$(diff "$ROOT_DIR"/old.fa "$ROOT_DIR"/new.fa)

if [ -z "$difference" ]; then
    send_email "<font size='64'>&#x2713;</font>"
else
    send_email "<font size='64'>&#x2620;</font><br><br><pre>$difference</pre>"
fi
