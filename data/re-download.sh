awk 'NR>1 {print $1}' filereport_read_run_PRJNA722257.tsv | while read SRR_ID; do
    echo "Downloading and splitting $SRR_ID..."
    fasterq-dump --include-technical --split-files "$SRR_ID"
done
