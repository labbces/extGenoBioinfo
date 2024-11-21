#!/bin/bash

cd ~/
rm -rf dia4 dia5 dia3 dia6
mkdir -p ~/dia6
cd ~/dia6
wget https://labbces.cena.usp.br/shared/CEN5789/dia6/NRRLY27205.asm.reads.sorted.bam
wget https://labbces.cena.usp.br/shared/CEN5789/dia6/NRRLY27205.asm.reads.sorted.bam.bai
wget https://labbces.cena.usp.br/shared/CEN5789/dia6/NRRLY27205.asm.bp.hap1.p_ctg.softmasked.fa
wget https://labbces.cena.usp.br/shared/CEN5789/dia6/GALBA/galba.gtf
wget https://labbces.cena.usp.br/shared/CEN5789/dia6/GALBA_EGGNOG.emapper.decorated.gff