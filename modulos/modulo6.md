# Módulo 6 - Anotação de Genomas: Predição de Genes

Certifique-se de que todos os exercícios desta seção sejam executados dentro da pasta `~/dia6`. Se essa pasta não existe, crie-a antes de começar os exercícios. Isso garantirá que você esteja organizando seus arquivos da maneira apropriada.

### Realize a avaliação de contaminações.

Antes de prosseguirmos com a anotação estrutural, é crucial aprimorar nossa avaliação das montagens. Para isso, é de extrema importância realizar uma análise detalhada da possível contaminação por meio da ferramenta [BlobToolKit](https://github.com/blobtoolkit/blobtoolkit). Essa etapa adicional de avaliação nos permitirá garantir a qualidade e confiabilidade das montagens, fornecendo uma visão abrangente da integridade dos dados genômicos. Por razões de tempo, optaremos por pular esta etapa hoje, mas é fundamental ressaltar que em projetos reais, a execução desta análise é essencial. Hoje só iremos estudas as figuras do artigo [No evidence for extensive horizontal gene transfer in the genome of the tardigrade Hypsibius dujardini](https://www.pnas.org/doi/full/10.1073/pnas.1600338113). Discuta com seu instrutor.

### Mascare repetições

Antes de realizar a predição de genes, é imperativo mascarar os scaffolds/contigs da montagem por meio de uma biblioteca de repetições adaptada ao genoma em questão. Entre as estratégias frequentemente empregadas, destacam-se o uso do [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler)/[RepeatMasker](https://www.repeatmasker.org/), [EDTA](https://github.com/oushujun/EDTA) e [Earl Grey](https://github.com/TobyBaril/EarlGrey). No entanto, é importante observar que esse processo costuma demandar considerável capacidade computacional. Por isso, hoje optaremos por omiti-lo.

Se você já possui uma boa biblioteca das repetições presentes em seu genoma, pode utilizar o [NGSEP](https://github.com/NGSEP/NGSEPcore) como um mascarador rápido, nos modos TransposonsFinder e GenomeAssemblyMask. Para isso, usaremos a biblioteca do [Dfam](https://www.dfam.org/) que temos disponível em formato [FASTA](https://labbces.cena.usp.br/shared/CEN5789/dia6/Dfam_curatedonly.fasta). Observe que, primeiro, é necessário acessar o repositório do NGSEP e fazer o download do aplicativo através do link de Releases. Baixe tanto o software quanto a biblioteca do Dfam na pasta "dia6". Caso a pasta não exista, crie-a.

```
conda activate redotable
mkdir -p ~/dia6
cd ~/dia6
wget https://labbces.cena.usp.br/shared/CEN5789/dia6/Dfam_curatedonly.fasta
wget https://github.com/NGSEP/NGSEPcore/releases/download/v5.0.0/NGSEPcore_5.0.0.jar
java -jar NGSEPcore_5.0.0.jar TransposonsFinder -i NRRLY27205.asm.bp.hap1.p_ctg.g100kbp.fasta -o NRRLY27205.asm.bp.hap1.p_ctg.g100kbp.repeats -d Dfam_curatedonly.fasta -t 4
#Gerando uma versão soft-masked do genoma, com as repetições em letras minúsculas
java -jar NGSEPcore_5.0.0.jar GenomeAssemblyMask -i NRRLY27205.asm.bp.hap1.p_ctg.g100kbp.fasta -o NRRLY27205.asm.bp.hap1.p_ctg.g100kbp.softmasked.fa -d NRRLY27205.asm.bp.hap1.p_ctg.g100kbp.repeat
#Gerando uma versão hard-masked do genoma, substituindo as bases das repetições pela letra "N"
java -jar NGSEPcore_5.0.0.jar GenomeAssemblyMask -i NRRLY27205.asm.bp.hap1.p_ctg.g100kbp.fasta -o NRRLY27205.asm.bp.hap1.p_ctg.g100kbp.hardmasked.fa -d NRRLY27205.asm.bp.hap1.p_ctg.g100kbp.repeats -h 
conda deactivate
```

Quantas bases foram mascaradas? Podemos usar o programa compseq do EMBOSS para verificar isso:

```
conda activate emboss
compseq -word 1 -outfile stdout NRRLY27205.asm.bp.hap1.p_ctg.hardmasked.fa 
conda deactivate
```

Você acha que essa forma de mascaramento foi apropriada? Dica: Não, não foi suficiente.

### Obter evidência extrínseca - Proteínas de espécies próximas.

Vamos anotar apenas a montagem que apresentar as melhores métricas de completude e continuidade.

Faça o download de todas as proteínas do mesmo gênero no NCBI. Para fazer isso, acesse o banco de dados de taxonomia do NCBI em seu navegador e procure pelo nome do gênero _Kazachstania_ na lista de nomes. Clique no nome do gênero, conforme mostrado na figura:

![GetProteinsFromTaxonomy_1](images/GetProteinsFromTaxonomy_1.png)

Na página em que você chegou, identifique a tabela que mostra a figura, onde todos os registros associados a esse gênero em outros bancos de dados do NCBI estão listados. Clique na seção "Proteins", onde deve haver um número aproximado de 53.000 proteínas.

![GetProteinsFromTaxonomy_2](images/GetProteinsFromTaxonomy_2.png)

Finalmente descarregue um arquivo com essas proteinas em formato fasta.

![GetProteinsFromTaxonomy_2](images/GetProteinsFromTaxonomy_2.png)

### Anotando o Genoma - GALBA

Vamos a anotar o genoma usando [GALBA](https://github.com/Gaius-Augustus/GALBA). Outras alternativas incluem o uso do [BRAKER](https://github.com/Gaius-Augustus/BRAKER) ou [EASEL](https://gitlab.com/PlantGenomicsLab/easel), porém, essas ferramentas exigem dados de RNA-Seq e são mais intensivas em termos computacionais. Em cenários reais, a recomendação é empregar várias estratégias e gerar um conjunto de genes previstos com base nos melhores resultados das ferramentas, utilizando, por exemplo [EVidenceModeler](https://github.com/EVidenceModeler/EVidenceModeler). Sempre é importante utilizar evidência extrínseca, e na maioria dos casos, o RNA-Seq é a fonte de dados que oferece os melhores resultados.

Vamos usar um container do singularity para rodar mais facilmente o GALBA, para que isso funcione linque os arquivos de montagem do genoma e as proteínas para o diretório HOME.

```
conda activate singularitycew
#singularity build galba.sif docker://katharinahoff/galba-notebook:latest
wget https://labbces.cena.usp.br/shared/CEN5789/dia4/galba.sif
singularity shell -B $PWD:$PWD galba.sif
cp -r $AUGUSTUS_CONFIG_PATH/ /home/cen5789/dia5/augustus
export AUGUSTUS_CONFIG_PATH=/home/cen5789/dia5/augustus
galba.pl --threads=10 --species=KazachstaniaBulderi --genome=NRRLY27205.asm.bp.hap1.p_ctg.g100kbp.softmasked.fa --prot_seq=sequence.fasta
exit
conda deactivate
```

Você pode encontrar o resultado da previsão de genes na pasta GALBA. A saída do programa GALBA é gerada principalmente em três arquivos distintos:

- galba.gtf
- galba.codingseq
- galba.aa

A previsão de genes para GALBA pode levar um tempo considerável. Portanto, estou disponibilizando os principais arquivos de resultados neste [arquivo](modulos/files/GALBA_annotation.tar.gz).

Por favor, analise esses arquivos para compreender o conteúdo presente. Discuta com seus colegas e seu professor para obter uma compreensão completa.

Após o primeiro passo de anotação estrutural do genoma, é fundamental avaliar a qualidade da anotação. Uma maneira de realizar essa avaliação é examinando a completude. Nesse contexto, esperamos que o nível de completude da anotação seja pelo menos tão bom quanto a análise da completude do espaço gênico durante a avaliação do genoma. Lembre-se de que, anteriormente, avaliamos o genoma com o software `compleasm`. No entanto, a versão atual desse programa é destinada exclusivamente à avaliação de genomas. Para avaliar a previsão de genes, utilizaremos o [BUSCO](https://busco.ezlab.org/). Para isso, faremos uso de uma imagem do BUSCO executando com o Singularity, o que simplifica significativamente a instalação do software.

```
conda activate singularitycew
#singularity build busco.sif docker://ezlabgva/busco:v5.5.0_cv1
wget https://labbces.cena.usp.br/shared/CEN5789/dia6/busco.sif
singularity shell -B $PWD:$PWD busco.sif
busco -i GALBA/galba.aa -o GALBA_BUSCO -m protein -l saccharomycetes_odb10 --cpu 10
exit
conda deactivate
```

Os resultados do BUSCO estão disponíveis na pasta GALBA_BUSCO. Por favor, examine os vários arquivos e discuta-os com seus colegas e seu professor. Compare os resultados do BUSCO das proteínas previstas com os resultados do `compleasm` para o genoma montado. Quantos genes foram preditos?

Com a conclusão da anotação estrutural, estamos prontos para iniciar a anotação funcional. Para isso, faremos uso do banco de dados [EGGNOG](http://eggnog5.embl.de/) e da ferramenta de software [Eggnog-Mapper](http://eggnog-mapper.embl.de/). Procederemos à anotação das proteínas preditas, ou seja, dos produtos gênicos, utilizando os dados relacionados aos fungos do EGGNOG.". Antes de avançarmos nos exercícios, é essencial realizar a conversão do arquivo de anotação gerado por GALBA, que está no formato GFF, para o formato GTF. Faremos essa conversão utilizando a ferramenta `gffread`. Isso é um passo importante antes de prosseguir.

```
conda activate eggnogmapper
export EGGNOG_DATA_DIR=/home/cen5789/dia6
gffread --keep-genes -o GALBA/galba.gff3 GALBA/galba.gtf
download_eggnog_data.py -P  -y
emapper.py  -m diamond --cpu 10 --itype proteins -i GALBA/galba.aa -o GALBA_EGGNOG --decorate_gff GALBA/galba.gtf --target_orthologs all --tax_scope 4751
conda deactivate
```

O que significa o número 4751? Consulte o banco de dados de [taxonomia do NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=4751).

É importante observar que o EGGNOG Mapper pode exigir recursos computacionais significativos, especialmente em termos de tempo de execução. Os arquivos de resultados estão disponíveis [aqui](modulos/files/GALBA_EGGNOG_FUNGI.tar.gz).

Agora, procederemos à visualização da montagem, juntamente com as leituras mapeadas nela e a anotação estrutural do genoma, usando o [Integrative Genomics Viewer (IGV)](https://igv.org/). Primeiro vamos mapear as leituras no genoma usando o `minimap2` e o `samtools`. Favor fazer uma cópia do seu arquivo de leituras na pasta de trabalho 'dia6'.

```
conda activate jupiterplot
minimap2 -H -x map-hifi -a -t 10 NRRLY27205.asm.bp.hap1.p_ctg.softmasked.fa SRR25033384.filt.fastq.gz | samtools view -b --fast --threads 6 |samtools sort --threads 6 -o NRRLY27205.asm.reads.sorted.bam
samtools index NRRLY27205.asm.reads.sorted.bam
conda deactivate
```

Agora podemos inicializar o IGV, que tem ambiente gráfico. Primeiramente, carregaremos a montagem no menu __"File"__ -> __"Load From File"__, utilizando o arquivo `NRRLY27205.asm.bp.hap1.p_ctg.softmasked.fa`. Em seguida, utilizando o mesmo menu, procederemos ao carregamento do arquivo com as leituras mapeadas, denominado `NRRLY27205.asm.reads.sorted.bam`. Por fim, carregaremos o arquivo contendo a anotação do genoma a partir de `GALBA/galba.gtf`.

```
conda activate igv
igv
```

Por exemplo, localize o contig h1tg000003l nas posições que vão de 551,742 até 578,200. Quais são os significados das regiões coloridas nas leituras? E qual é a interpretação das regiões roxas?" O que representa a região marcada com o número 495? É possível que haja um gene anotado nesta região? Aqui pode consultar a documentação do [IGV](https://igv.org/doc/desktop/#UserGuide/tracks/alignments/viewing_alignments_basics/). Tente colorir os alinhamentos com base na fita de origem da leitura.

![IGV screenshot](images/igv_snapshot.png)

Combine os contigs com mais de 100 kbp dos dois haplótipos (hap1 e hap2) em um único arquivo e, posteriormente, mapeie as leituras contra este arquivo. Você consegue identificar alguma diferença nos resultados?