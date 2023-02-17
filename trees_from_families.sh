#!/bin/bash

for fasta in *.fasta
do
	if [ $(grep -c ">" ${fasta}) -gt 1 ]
		then
			echo "Processing ${fasta}..." >> trees_from_families_log.txt

			PHY=${fasta//fasta/phy}
			TREE=${fasta//fasta/tree}

			echo "Creating phylip file.." >> trees_from_families_log.txt
			clustalw ${fasta} -output=PHYLIP -outfile=${PHY}
			
			echo "Creating tree..." >> trees_from_families_log.txt
			fasttree -gtr -gamma -cat 4 -nt ${PHY} > ${TREE}
	fi
done

echo "Done." >> trees_from_families_log.txt
