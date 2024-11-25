# Módulo 8 - Transcriptômica: Análise de Expressão Diferencial de Genes, Perfis de Genes

Vamos analisar a resposta dos genes de _Arabidopsis thaliana_ considerando dois fatores: 1) o genótipo e 2) o estresse ambiental. O fator genótipo é composto por três níveis: 1.a) selvagem, 1.b) mutante no gene ros1-3, 1.c) duplo mutante nos genes dml2 e dml3, e 1.d) triplo mutante nos genes ros1, dml2 e dml3. Enquanto o fator estresse ambiental possui também quatro níveis: 2.a) sem tratamento, 2.b) tratamento com ácido abscísico, 2.c) tratamento com cloreto de sódio e 2.d) tratamento com seca. [Mais detalhes](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-4477/). Os transcritos em cada condição foram sequenciados usando a tecnologia Illumina, resultando na geração de leituras em pares (paired-end), e os dados estão disponiveis no [SRA](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=282863).

Crie uma pasta com o nome `dia8` dentro do diretório $HOME. Todos os exercícios de hoje devem ser realizados dentro dessa pasta.

```
cd
mkdir dia8
cd  dia8
```

### Descarregando os dados de repositórios públicos

Anteriormente, utilizamos funcionalidades do SRA Toolkit para baixar dados do SRA. Agora, vamos empregar outra ferramenta que pode agilizar ainda mais o processo em comparação com o fasterq-dump que usamos anteriormente. Utilizaremos a ferramenta [`ffq`](https://github.com/pachterlab/ffq), a diferença, em relação ao fasterq-dump, faz o download dos dados diretamente no formato fastq. 

Vamos analisar dados de RNASeq de mutantes de _Arabidopsis thaliana_ em diversas condições ambientais. Os mutantes ros1-3, dml2 e dml3 são alelos mutantes que resultam na perda de função dos genes ROS1, DML2 e DML3, que fazem parte da família DEMETER (DME). Essa família está envolvida no processo de desmetilação do DNA.

O ROS1 atua como um repressor do silenciamento da expressão gênica. Sua função envolve a remoção da metilação do DNA do promotor do gene alvo. Ele interage fisicamente com RPA2/ROR1. Nos mutantes _ros1_ (como _ros1-3_), é observado um aumento na metilação nos promotores de vários genes. Dentre os loci afetados por _ros1_, alguns (RD29A e At1g76930) sofrem alterações na metilação de citosina em todos os contextos de sequência (CpG, CpNpG ou CpNpN), embora muitos outros sejam afetados principalmente em contextos não-CpG. O mutante _ros1_ demonstra maior suscetibilidade a patógenos biotróficos e possui redução em sua capacidade de resposta aos genes de defesa dependentes do ácido salicílico. O _ros1-3_ é um alelo de T-DNA, onde o T-DNA está inserido no gene ROS1, gerando a perda de funcão. ([Penterman, et al., 2007](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1847597/), [TAIR - AT2G36490](https://arabidopsis.org/servlets/TairObject?id=32321&type=locus)).

Os genes DML2 ([AT3G10010](https://arabidopsis.org/servlets/TairObject?id=40110&type=locus)) e DML3 ([AT4G34060](https://arabidopsis.org/servlets/TairObject?id=127923&type=locus)) são homólogos de ROS1. Mutacões nos genes DML resultam em hipermetilação do DNA em locais específicos. Dos locais desmetilados pelas enzimas DML, mais de 80% estão próximos ou se sobrepõem a genes. A desmetilação gênica pelas enzimas DML ocorre principalmente nas extremidades 5' e 3', um padrão oposto à distribuição geral da metilação de DNA do tipo selvagem ([Penterman, et al., 2007](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1847597/)).


A tabela a seguir apresenta uma lista dos números de acesso do SRA para cada amostra, juntamente com uma descrição das condições experimentais às quais foram submetidas.

| Number | Identifiers | Genotype | Environmental Stress | 
| --- | --- | --- | --- |
| 01 | DRR016125 | wild type | None |
| 02 | DRR016126 | wild type | abscisic acid |
| 03 | DRR016127 | wild type | sodium chloride |
| 04 | DRR016128 | wild type | drought |
| 05 | DRR016129 | ros1-3 mutant | none |
| 06 | DRR016130 | ros1-3 mutant | abscisic acid |
| 07 | DRR016131 | ros1-3 mutant | sodium chloride |
| 08 | DRR016132 | ros1-3 mutant | drought |
| 09 | DRR016133 | dml2, dml3 double mutant | none |
| 10 | DRR016134 | dml2, dml3 double mutant | abscisic acid |
| 11 | DRR016135 | dml2, dml3 double mutant | sodium chloride |
| 12 | DRR016136 | dml2, dml3 double mutant | drought |
| 13 | DRR016137 | ros1, dml2, dml3 triple mutant | none |
| 14 | DRR016138 | ros1, dml2, dml3 triple mutant | abscisic acid |
| 15 | DRR016139 | ros1, dml2, dml3 triple mutant | sodium chloride |
| 16 | DRR016140 | ros1, dml2, dml3 triple mutant | drought |

Vamos descarregar os links de acceso dos arquivos em formato `fastq.gz` para a mostra DRR016125. Para isso, é necessário ativar o ambiente Conda denominado "transcriptomics," no qual estão instalados todos os softwares que utilizaremos nas próximas semanas.

```
conda activate transcriptomics
ffq --ftp DRR016125
```

Isso deve gerar uma saída semelhante a esta:

```
[2023-11-07 00:19:53,630]    INFO Parsing run DRR016125
[
    {
        "accession": "DRR016125",
        "filename": "DRR016125_1.fastq.gz",
        "filetype": "fastq",
        "filesize": 967794560,
        "filenumber": 1,
        "md5": "36beaacaff3a9903d5c57a2c73525ae8",
        "urltype": "ftp",
        "url": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR016125/DRR016125_1.fastq.gz"
    },
    {
        "accession": "DRR016125",
        "filename": "DRR016125_2.fastq.gz",
        "filetype": "fastq",
        "filesize": 1001146319,
        "filenumber": 2,
        "md5": "d3f0203f534e35fa370b2a6c5e238944",
        "urltype": "ftp",
        "url": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR016125/DRR016125_2.fastq.gz"
    }
]
```

Observe que existem duas linhas que começam com "url":, que são os endereços dos arquivos fastq.gz na internet, localizados no servidor do SRA do Instituto Europeu de Bioinformática. Você pode usar esses endereços com o programa "wget" ou "curl" para fazer o download dos arquivos para o seu computador. Isso precisaria ser repetido para cada uma das amostras deste experimento.

Como esses arquivos são pesados, o professor já os baixou em um servidor do [CENA](https://labbces.cena.usp.br/CEN5789/transcriptomics/RAWREADS). Os arquivos contendo as leituras de todas as amostras requerem aproximadamente 30GB de armazenamento, um espaço que talvez não esteja disponível nos computadores que estamos usando. Portanto, cada aluno fará o download de apenas um par de arquivos (R1 e R2) de uma única amostra e trabalhará apenas com eles. Dessa forma, iremos paralelizar nosso trabalho, executando os processos de verificação de qualidade, limpeza e quantificação. Após a quantificação, compartilharemos os resultados de modo que todos os alunos tenham acesso às quantificações de todos os genes em todas as amostras. Siga as instruções do professor. 

Faça o download dos seus dados na pasta "RAWREADS" dentro da pasta "dia8". Se essa pasta ainda não existir, crie-a. Lembre que estamos usando o ambiente conda chamado `transcriptomics`. Lembre-se também de substituir o identificador da SUA amostra.

```
mkdir -p ~/dia8/RAWREADS
cd ~/dia8/RAWREADS
ID=DRR016125
curl -O https://labbces.cena.usp.br/shared/extGenoBioinfo/dia8/RAWREADS/${ID}_[1-2].fastq.gz

```

### Pre-processando os dados de RNASeq

Vamos a conferir a qualidade do sequenciamento usando o programa `fastqc`

```
cd ~/dia8/
mkdir FastQC_pre
fastqc --threads 2 --nogroup  --outdir FastQC_pre RAWREADS/${ID}_[1-2].fastq.gz
```

Visualize os resultados e tome as decisões necessárias para realizar a limpeza das leituras. Lembre-se de que as bibliotecas dessas amostras foram criadas usando a tecnologia [TruSeq](https://www.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/datasheet_truseq_sample_prep_kits.pdf), que pesca mRNA poliadenilados, e o cDNA foi gerado com iniciadores aleatórios (_random primers_).

Agora, vamos proceder com a limpeza usando o programa `bbduk` da suíte "bbmap". A qualidade das leituras, em termos da distribuição do escore Phred, está muito boa, e acredito que não seja necessário realizar o processo de _quality trimming_. No entanto, foi detectada a presença residual de adaptadores no extremo 3' em algumas amostras, os quais precisam ser removidos. Além disso, é aconselhável quantificar a quantidade de leituras que têm origem no rRNA, pois isso pode fornecer uma indicação da qualidade das amostras ainda nesta etapa de limpeza.

Primeiro removemos adaptadores:

```
cd ~/dia8
bbduk.sh in=RAWREADS/${ID}_1.fastq.gz in2=RAWREADS/${ID}_2.fastq.gz out=CLEANREADS/${ID}_cleana_1.fastq.gz out2=CLEANREADS/${ID}_cleana_2.fastq.gz ref=adapters refstats=CLEANREADS/${ID}_cleana_adapters_refstats ktrim=r threads=10
```

Agora, vamos filtrar (excluir) as leituras que correspondem ao rRNA, utilizando como entrada as leituras nas quais os adaptadores foram removidos no passo anterior. Mas antes de prosseguirmos, é necessário fazer o download do banco de dados contendo as sequências de rRNA. Este banco é derivado do [SILVA NR](https://www.arb-silva.de/), e as sequências foram agrupadas com 90% de identidade.

```
cd ~/dia8
mkdir -p ~/dia7/References
cd ~/dia8/References
wget https://labbces.cena.usp.br/shared/extGenoBioinfo/dia8/References/rRNA.tar.gz
tar xvzf rRNA.tar.gz
rm -rf rRNA.tar.gz
cd ..
bbduk.sh in=CLEANREADS/${ID}_cleana_1.fastq.gz in2=CLEANREADS/${ID}_cleana_2.fastq.gz out=CLEANREADS/${ID}_cleanf_1.fastq.gz out2=CLEANREADS/${ID}_cleanf_2.fastq.gz ref=References/rRNA_LSU_SILVA_Archaea.nr90.fasta,References/rRNA_LSU_SILVA_Bacteria.nr90.fasta,References/rRNA_LSU_SILVA_Eukarya.nr90.fasta,References/rRNA_SSU_SILVA_Archaea.nr90.fasta,References/rRNA_SSU_SILVA_Eukarya.nr90.fasta,References/rRNA_SSU_SILVA_Bacteria.nr90.fasta ktrim=f threads=10  minlength=85 refstats=CLEANREADS/${ID}_cleanf_rRNA_refstats
```

Confira o arquivo `*_cleanf_rRNA_refstats` dentro da pasta `CLEANEADS`. Uma proporção elevada de leituras de rRNA pode indicar problemas com a amostra. Se teve algum problema realizando a limpeza das leituras, pode descarregar os arquivos já limpos [aqui](https://labbces.cena.usp.br/shared/extGenoBioinfo/dia8/CLEANREADS/).

### Estimando os Níveis de Expressão de Transcritos e Genes

Com as leituras limpas em mãos, podemos começar a planejar a estimativa dos níveis de expressão dos transcritos e/ou genes. Para esta tarefa, utilizaremos o programa [Salmon](https://salmon.readthedocs.io/en/latest/), que emprega a estratégia de "quasi-mapping", conhecida por sua alta precisão e rapidez. É importante notar que o Salmon realiza a comparação em relação a uma referência que consiste nas sequências dos transcritos de interesse e pode lidar com a estimativa de valores de expressão para sequências muito semelhantes.

Vamos baixar a versão mais recente do Salmon e instalá-la na pasta `~/dia8/`:

```
cd ~/dia8
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
tar xvzf salmon-1.10.0_linux_x86_64.tar.gz
rm salmon-1.10.0_linux_x86_64.tar.gz
export PATH=/home/$USER/dia8/salmon-latest_linux_x86_64/bin/:$PATH
salmon
```

A última linha do bloco anterior deveria ter exibido a ajuda do programa, algo similar a:

```
salmon v1.10.0

Usage:  salmon -h|--help or
        salmon -v|--version or
        salmon -c|--cite or
        salmon [--no-version-check] <COMMAND> [-h | options]

Commands:
     index      : create a salmon index
     quant      : quantify a sample
     alevin     : single cell analysis
     swim       : perform super-secret operation
     quantmerge : merge multiple quantifications into a single file

```

A referência que usaremos é composta por todos os transcritos (cDNAs) anotados no genoma de _Arabidopsis thaliana_, os quais podem ser baixados do [TAIR](https://www.arabidopsis.org), também pode encontrar o arquivo [aqui](https://labbces.cena.usp.br//CEN5789/transcriptomics/References/TAIR10_cdna_20101214_updated.gz), descarreguelo dentro da sua pasta `~/dia8/` e dentro de uma subpasta chamada `References`. Esta referência precisa ser complementada com sequências decoy, ou seja, sequências que não deveriam estar presentes para a quantificação. Neste caso, usaremos o genoma completo como decoy. Recomendo a leitura [deste articulo](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8) para entender a importância do uso de decoy em análises de RNASeq.

```
mkdir -p ~/dia8/References
cd ~/dia8/References
#Os transcritos
wget https://labbces.cena.usp.br/shared/extGenoBioinfo/dia8/References/TAIR10_cdna_20101214_updated.gz
#O genoma
wget https://labbces.cena.usp.br/shared/extGenoBioinfo/dia8/References/TAIR10_genome.fasta.gz
```

Vamos gerar um arquivo de texto contendo os identificadores das sequências que serão usadas como decoy:

```
grep "^>" <(gunzip -c TAIR10_genome.fasta.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
```

Em seguida, podemos construir o índice para que o Salmon possa realizar a quantificação dos transcritos:

```
cat TAIR10_cdna_20101214_updated.gz TAIR10_genome.fasta.gz > gentrome.fa.gz
salmon index -t gentrome.fa.gz -d decoys.txt -p 10 -i salmon_index
cd ../
```

Finalmente, podemos iniciar o processo de quantificação das amostras propriamente dito:

```
mkdir -p ~/dia8/quantification
cd ~/dia8/
salmon quant -i References/salmon_index -l A -1 CLEANREADS/${ID}_cleanf_1.fastq.gz -2 CLEANREADS/${ID}_cleanf_2.fastq.gz --validateMappings -o ~/dia8/quantification/${ID} --threads 10 --seqBias --gcBias
```

Isso criará a pasta `quantification/${ID}`. Dentro dela, por favor, verifique o arquivo `logs/salmon_quant.log` e identifique o tipo de biblioteca detectado pelo Salmon e a taxa de mapeamento. Vamos tabular os resultados [neste arquivo](https://docs.google.com/spreadsheets/d/1EcOg9kpVz2dfL6DBEYUMTtU_jfslNgPpBeegi6N_kWc/edit?usp=sharing).

A quantificacao está disponivel no arquivo `~/dia8/quantification/${ID}/quant.sf`. Identifique as informações contidas no arquivo.

É necessário realizar esse processo para cada uma das 15 amostras restantes. Hoje, seu professor já executou esse procedimento, e você pode baixar a quantificação completa através do seguinte link:

```
cd ~/dia8/
wget https://labbces.cena.usp.br/shared/extGenoBioinfo/dia8/quantification/quantification.tar.gz
```

### Vamos nos preparar para a detecção de GENES diferencialmente expressos.

Lembre-se de que a quantificação realizada pelo Salmon foi de transcritos, e a partir de um único gene, vários transcritos podem ser produzidos. Hoje, vamos realizar a detecção de GENES diferencialmente expressos nas condições experimentais. Por isso, precisamos indicar quais transcritos são originados a partir do mesmo gene, a fim de transformar as quantificações de cada transcrito em quantificações de cada gene. Basicamente, iremos somar as quantificações de todos os transcritos do mesmo gene.

Você pode estar interessado em detectar os transcritos diferencialmente expressos, o que é outro tipo de análise. Isso requer a execução do Salmon com diferentes parâmetros e o uso de um software estatístico diferente.

Então, agora precisamos obter um mapa que relacione o identificador do transcrito com o identificador do gene que o origina. Felizmente, com os identificadores dos transcritos de _A. thaliana_, isso é muito fácil. Esses identificadores têm o formato:

AT**X**G**YYYYY**.**Z**

Onde **X** é o cromossomo e pode ter os valores 1, 2, 3, 4, 5, M ou C. **YYYYY** é um número serial. Essa cadeia de texto **ATXGYYYYY** representa o identificador do gene.

```
zcat TAIR10_cdna_20101214_updated.gz |grep  ">" | cut -f 1 -d' '|sed 's/>//'| sed -r 's/((AT.G[0-9]*)\.[0-9]*)/\1\t\2/' > tx2gene.txt
```

### Carregando os dados em R

Vamos usar o RStudio, para isso, no seu terminal, execute o comando:

```
rstudio
```

Uma janela como a seguinte deveria aparecer na sua tela:

![Screendhot Rstudio](images/rstudio.png)

A grosso modo, podemos dividir a área de trabalho do RStudio em quatro regiões:

- Na área superior esquerda, você encontrará abas que hospedam seus scripts e projetos. Isso permite que você organize e acesse facilmente seu código-fonte e trabalhe em diferentes projetos de análise de dados.
- Na área inferior esquerda, estão localizadas as abas Console, Terminal e Background Jobs, que permitem interagir com o ambiente do R. O Console é onde você pode executar comandos e ver a saída imediata. O Terminal permite interações mais avançadas com o sistema operacional, enquanto as Tarefas em Segundo Plano rastreiam e gerenciam tarefas em execução.
- A área superior direita exibe as abas Environment, History, Connections e Tutorial. São abas que fornecem informações e recursos adicionais. O Ambiente mostra os objetos em sua sessão R atual. O Histórico registra comandos anteriores. Conexões permitem a integração com outros serviços e tutoriais fornecem orientações úteis.
- A área inferior direita inclui as abas Files, Plots, Packages, Help, Viewer e Presentation. Arquivos gerencia seus arquivos e diretórios. Gráficos exibe gráficos gerados durante sua análise. Pacotes ajuda a gerenciar pacotes R. Ajuda fornece informações detalhadas sobre funções e pacotes. O Visualizador permite visualizar conjuntos de dados e gráficos interativamente. Apresentação ajuda a criar apresentações R Markdown.

Essas abas tornam o RStudio uma ferramenta versátil e eficaz para o desenvolvimento e análise de projetos em R.

No exemplo abaixo, criamos duas variáveis no console: greeting <- "Hola, Mundo!" e meuVector <- c(1, 2, 3, 4). Observe como elas são exibidas na aba Ambiente:

![Screendhot Rstudio2](images/rstudio2.png)

Para a análise de dados de transcriptômica, vamos criar um novo script no qual armazenaremos todos os passos da análise. Para começar a gravar um script, clique em Arquivo - Novo Arquivo - Script R. Isso abrirá um editor de texto no canto superior esquerdo da interface do RStudio (acima da guia Console).

Primeiramente, precisamos carregar algumas bibliotecas do R que estendem as funcionalidades básicas da linguagem e nos permitem lidar com dados de expressão gênica.

```R
library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(reshape2)
```
Recomendo que se acostume a limpar todos os objetos que estão presentes em seu ambiente. Possivelmente, neste momento, você tem apenas os objetos greeting e meuVector, mas se decidir reexecutar o script, é conveniente começar do zero.

```R
rm(list=ls())
```
Agora, vamos definir o seu diretório de trabalho. Ele deve ser o diretório que contém todos os resultados do Salmon. Lembre-se de que você deve ter 16 pastas com os resultados do Salmon.

```R
wd<-"~/dia8/quantification/"
setwd(wd)
```

Agora, antes de importar os dados de quantificação para o R, é crucial possuir uma descrição do plano experimental. Para carregar essa descrição, utilizamos a função `read.delim` para importar um arquivo de texto chamado [targets.txt](files/targets.txt) com colunas separadas por TAB. Os dados são armazenados em um objeto denominado `targets`, o qual é um `data.frame`.

```R
targets<-read.delim("targets.txt",header=T)
targets
rownames(targets)<-targets$SampleName
targets
targets$Genotype<-as.factor(targets$Genotype)
targets$EnvironmentalStress<-as.factor(targets$EnvironmentalStress)
targets$Number<-NULL
targets
```
Você pode inspecionar o objeto `targets` agora, basta digitar o nome do objeto em seu Console. Como alternativa, você pode clicar no nome do objeto na aba "Environment" localizada na área superior direita.

```R
> targets
#          Number SampleName     Genotype EnvironmentalStress
#DRR016125      1  DRR016125           WT                None
#DRR016126      2  DRR016126           WT                 ABA
#DRR016127      3  DRR016127           WT                NaCl
#DRR016128      4  DRR016128           WT             Drought
#DRR016129      5  DRR016129        ros13                None
#DRR016130      6  DRR016130        ros13                 ABA
#DRR016131      7  DRR016131        ros13                NaCl
#DRR016132      8  DRR016132        ros13             Drought
#DRR016133      9  DRR016133     dml2dml3                None
#DRR016134     10  DRR016134     dml2dml3                 ABA
#DRR016135     11  DRR016135     dml2dml3                NaCl
#DRR016136     12  DRR016136     dml2dml3             Drought
#DRR016137     13  DRR016137 ros1dml2dml3                None
#DRR016138     14  DRR016138 ros1dml2dml3                 ABA
#DRR016139     15  DRR016139 ros1dml2dml3                NaCl
#DRR016140     16  DRR016140 ros1dml2dml3             Drought
```
Precisamos carregar o arquivo que contém a correspondência entre identificadores de genes e transcritos que você criou anteriormente.

```R
tx2gene<-read.delim("tx2gene.txt",header=F)
head(tx2gene)
```
Vamos ver um exemplo de um gene que possui 3 transcritos.

```R
tx2gene[which(tx2gene$V2=='AT1G51370'),]
#               V1        V2
#1     AT1G51370.2 AT1G51370
#8384  AT1G51370.1 AT1G51370
#10296 AT1G51370.3 AT1G51370
```
Para nos prepararmos para carregar as quantificações do Salmon, precisamos criar um objeto que contenha os nomes das amostras e a localização no disco dos arquivos `quant.sf` produzidos pelo Salmon.

Os nomes das amostras estão no objeto "target" na coluna `SampleName`. Além disso, sabemos que no nosso diretório de trabalho temos uma pasta para cada amostra, com o mesmo nome da coluna `SampleName` e o arquivo `quant.sf` está dentro dessa pasta.

```R
myFiles<-paste(wd, targets$SampleName,"/quant.sf",sep="")
myFiles
```
Inspecione o objeto `myFiles` e verifique se os caminho dos arquivos estão corretos. Agora vamos a adicionar o nome da amostra como o nome de cada caminho.

```R
names(myFiles)<-targets$SampleName
myFiles
```

E verificamos se os arquivos de fato existem no disco. O resultado dessa operação deve ser TRUE na sua tela. Caso contrário, significa que algum arquivo não está presente ou que o caminho está errado.

```R
all(file.exists(myFiles))
```

Agora, vamos importar os dados de quantificação do Salmon para o R, utilizando o pacote tximport. É importante observar que, neste caso, o tximport irá resumir os dados de expressão no nível dos genes.

```R
txi.salmon<-tximport(files = myFiles, type = 'salmon', tx2gene = tx2gene, txIn = TRUE, txOut = FALSE)
```

Vamos conferir o objeto `txi.salmon`; ele contém várias informações diferentes, como as abundâncias dos genes em _TPM_ (_Transcripts Per Million_), seus comprimentos e as contagens (número de leituras/fragmentos por gene por amostra). Vamos dar uma olhada nas primeiras linhas do objeto.

```R
head(txi.salmon$counts)
```

Esta matriz é o ponto de entrada no DESeq2, onde cada linha representa um gene (_g_), e cada coluna uma amostra (_i_). As células $` K_{ij} `$ indicam o número de fragmentos sequenciados que foram observados para o gene na amostra.

Apenas para fins de comparação entre a quantificação no nível de gene e no nível de transcrito, carregaremos os dados no nível de transcritos.

```R
txi.salmon.tx<-tximport(files = myFiles, type = 'salmon', tx2gene = tx2gene, txIn = TRUE, txOut = TRUE)
```

Vamos examinar os dados para um gene e, em seguida, remover o objeto `txi.salmon.tx`.

```R
head(txi.salmon$counts['AT1G51370',])
head(txi.salmon.tx$counts[c('AT1G51370.1','AT1G51370.2','AT1G51370.3'),])
rm(txi.salmon.tx)
```

Observe que o valor de expressão do gene é a soma dos valores de expressão de seus transcritos. Isso ocorre porque o Salmon, ao realizar a quantificação, fornece estimativas de abundância para cada transcrito individualmente. Ao importar esses dados no R usando o pacote tximport, com os argumentos `txIn = TRUE, txOut = FALSE`, a função automaticamente realiza a agregação dos valores de expressão dos transcritos para calcular a expressão a nível de gene.

Essa abordagem é útil porque muitos experimentos de RNA-Seq fornecem dados de expressão em nível de transcrito, mas, muitas vezes, estamos mais interessados na expressão a nível de gene para análises mais globais. A soma dos valores de expressão dos transcritos para um gene específico oferece uma medida geral da expressão desse gene no contexto do experimento.

Portanto, ao usar Salmon para a quantificação e tximport para a importação e agregação dos dados no R, obtemos uma representação mais consolidada da expressão gênica em termos de valores agregados a nível de gene.

Agora, vamos avaliar se todas as amostras foram sequenciadas com a mesma intensidade (esforço de sequenciamento). Embora seja ideal que isso aconteça, na prática, devido a diversos fatores, é difícil alcançar essa uniformidade. Portanto, é necessário calcular fatores de normalização que levem em consideração a profundidade de sequenciamento de cada amostra.

```R
seqEffort<-as.data.frame(colSums(txi.salmon$counts))
seqEffort$Samples<-rownames(seqEffort)
colnames(seqEffort)<-c('NumberFragments', 'Samples')
seqEffort
```

Vamos gerar uma figura com esses dados utilizando o pacote [ggplot2](https://ggplot2.tidyverse.org/) (observe que já carregamos o pacote no início do script). O ggplot2 é uma poderosa biblioteca de visualização em R, amplamente utilizada para criar gráficos estatísticos de alta qualidade. Desenvolvido por [Hadley Wickham](https://en.wikipedia.org/wiki/Ggplot2), o ggplot2 segue uma abordagem de "[gramática de gráficos](https://link.springer.com/book/10.1007/0-387-28695-0)", o que significa que você constrói visualizações adicionando camadas de elementos gráficos.

Para criar um gráfico com ggplot2, você especifica os dados, mapeia variáveis às propriedades estéticas (como cor, forma ou tamanho), e adiciona camadas geométricas para representar os dados visualmente (pontos, linhas, barras, etc.). Além disso, é possível personalizar a aparência do gráfico ajustando temas, escalas e legendas.

Aqui, utilizamos a função `ggplot` com o data.frame `seqEffort`, definindo os nomes das amostras no eixo X e a quantidade de fragmentos sequenciados no eixo Y. Em seguida, especificamos que desejamos criar um gráfico de colunas com esses dados.

```R
p<- ggplot(seqEffort, aes(x=Samples,y=NumberFragments))+ 
  geom_col()
p
```

Esse gráfico inicial não é esteticamente muito agradável. Vamos aprimorá-lo removendo o fundo cinza e posicionando os nomes das amostras em um ângulo.

```R
p + theme_bw() +
  theme(axis.text.x = element_text(angle=60, size = 12, hjust = 1))
```
Isso já está muito melhor. Vamos agora alterar os nomes dos eixos para deixá-los em português.

```R
p + theme_bw() +
  theme(axis.text.x = element_text(angle=60, size = 12, hjust = 1)) +
  xlab("Amostras") +
  ylab("Número de fragmentos sequenciados")
```

Agora, por favor, discuta com seus colegas o efeito da profundidade de sequenciamento na estimação dos valores de expressão. Pode ser interessante examinar os dados de alguns genes, por exemplo:

```R
txi.salmon$counts['AT1G51370',]
```

### DESeq2

Existem diversos pacotes disponíveis para a identificação de Genes Diferencialmente Expressos (DEGs, por suas sigla em inglês) em R, sendo que [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) ([Love et al., 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)) e [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) ([Robinson et al., 2010](https://pubmed.ncbi.nlm.nih.gov/19910308/)) são dois dos mais amplamente utilizados e reconhecidos na comunidade científica. No entanto, é importante mencionar que há uma variedade de outras ferramentas e pacotes que também desempenham um papel significativo na análise de expressão gênica diferencial no ambiente R. A escolha do pacote pode depender de diversos fatores, como a natureza dos dados, as premissas do experimento e as preferências metodológicas do pesquisador.

Nesta seção, vamos utilizar o pacote DESeq2 e, inicialmente, apresentar o modelo estatístico empregado por este pacote. O DESeq2 utiliza um modelo linear generalizado baseado na distribuição binomial negativa das contagens de cada gene em cada condição ([Negative Binomial](https://en.wikipedia.org/wiki/Negative_binomial_distribution) [Generalized Linear Model](https://en.wikipedia.org/wiki/Generalized_linear_model)). Esse modelo considera a variabilidade intrínseca aos dados de RNA-Seq.


#### Modelo estatístico no DESeq2

No modelo de DESEq2, as contagens _y_ do gene _g_ na amostra _i_, são amostradas de uma distribuição binomial negativa com média _μ_ e parâmetro de dispersão _α_.

```math
y_{gi} ~ NB(μ_{gi}, α_g)
```

Uma estimativa precisa do parâmetro de dispersão $\alpha_{g}$ é crucial para uma correta inferência estatística. No entanto, não é simples devido ao baixo número de replicatas normalmente utilizadas em experimentos de RNA-Seq, o que resulta em estimativas de dispersão muito variáveis. Uma maneira de lidar com isso é compartilhar informações de vários genes. Inicialmente, obtemos estimativas de dispersão gene a gene usando máxima verossimilhança, com base nos dados de cada gene individualmente. Isso nos permite estimar a relação entre a dispersão e o valor médio das contagens, ou seja, obter uma estimativa precisa para o valor da dispersão esperada. No entanto, essa abordagem não captura as variações individuais de cada gene em relação a essa tendência. Em seguida, ajustamos as estimativas de dispersão gene a gene em direção aos valores previstos pela tendência para obter os valores finais de dispersão, em um processo chamado de "shrinkage". Dessa forma, o procedimento captura a variação específica de cada gene tanto quanto os dados permitem estimar. Ao mesmo tempo, a curva de tendência facilita a estimativa e inferência em regiões com menos conteúdo de informação. Esse método equilibra a precisão das estimativas individuais com a estabilidade proporcionada pela tendência global.

A abordagem de compartilhamento de informação é particularmente útil em experimentos com um pequeno número de replicatas, pois permite que genes semelhantes forneçam informações uns aos outros, melhorando assim a estabilidade das estimativas de dispersão. Essa estratégia é implementada em métodos como o DESeq2 para aprimorar a precisão das inferências em experimentos de RNA-Seq.
O modelo implementado no DESeq tem algumas suposições (todos os modelos têm) importantes e só deveria ser usado se você tiver certeza de que seus dados satisfazem essas suposições. Essas suposições são:

- As observações (contagens) são assumidas como independentes entre si.
- O parametro de dispersão _α_ é constante para todos os genes.
- O média das contagens para um gene numa amostra ($` μ_{gi} `$) está diretamente relacionada a verdadeira abundancia desse gene ($` q_{gi} `$), ajustada (_scaled_) por um fator especifico (_size factor_) para a amostra  $` s_i `$.

```math
μ_{gi} =: q_{gi} * s_i
```

O parametro $` s_i `$, pode incorporar a profundidade da amostra, a composicao da amostra, etc., i.e., todos os fatores que precisem ser normalizados.  O DESeq2 estima esse parametro para cada amostra, e o utiliza para normalizar as contagens de cada gene. O DESeq2 utiliza o modelo de regressão log-linear para estimar os parametros $` s_i `$ e $` q_{gi} `$. 

Geralmente, os _size factors_ incluem principalmente a profundidade do sequenciamento e assumem que uma minoria de genes está sendo afetada pelas condições experimentais. Se essa suposição não for correta, os _size factors_ podem ser estimados de outras maneiras. Por exemplo, quando há informações prévias, pode-se utilizar um conjunto de genes que não deveriam mudar em sua expressão em relação às condições experimentais. O DESeq2 possui uma função específica para isso chamada controlGenes.

Essa abordagem, utilizando genes de controle, proporciona uma maneira mais personalizada de estimar os _size factors_, levando em consideração genes específicos que se espera que não sejam influenciados pelas condições experimentais. Isso pode ser particularmente útil em experimentos nos quais a maioria dos genes não atende à suposição padrão de não serem afetados pelas condições experimentais

```math
log_2(q_{gi}) = \sum_r x_{ri} \beta_{rg}
```

Onde, $` x_{ri} `$  é uma matriz com o planejamento experimental e $` \beta_{rg} `$ é o coeficiente de regressão para o gene _g_ na amostra _i_, é está relacionao a mudança de abundânçia (fold change) do gene _g_ na amostra _i_ em relação a uma amostra de referência. 

O que representa essa matriz de planejamento experimental? Vamos explorar um experimento simples com um único fator e dois níveis desse fator, ou seja, controle e tratamento, cada nível com duas replicatas biológicas. Na matriz apresentada, as linhas representam as unidades experimentais, enquanto cada coluna corresponde a um nível do fator em estudo. Os valores nas células são 0 ou 1, indicando se uma unidade experimental foi (1) ou não (0) atribuída a um dos níveis do fator.

```math
x_{ri} = \begin{bmatrix} 1 & 0 \\\ 1 & 0 \\\ 0 & 1 \\\ 0 & 1 \end{bmatrix}
```

Os parametros $` \beta `$ estão no vector:

```math
\beta_g = \begin{bmatrix} \beta_{g0} \\\ \beta_{g1} \end{bmatrix}
```

Depois de ajustar o modelo com os dados disponiveis, o seja calculados os coeficientes $` \beta_{gi} `$, a mudança na expressão do gene pode ser expressa em escala de $` \log_2 `$ como 

```math
Log2FC(g) = (\beta_{g1} - \beta_{g0})
```

Ao utilizar esse modelo, o DESeq2 possibilita a identificação de genes diferencialmente expressos com maior precisão, controlando eficazmente a taxa de erro e fornecendo resultados estatisticamente significativos. Isso faz dele uma ferramenta poderosa na análise de expressão gênica diferencial, especialmente quando lidamos com dados complexos e experimentos de RNA-Seq.

#### Identificacão de genes diferencialmente expressos com DESseq2

Agora sim, vamos fornecer os dados ao DESeq2. Este pacote possui rotinas que podem importar os dados diretamente do objeto criado pelo ``tximport``. Durante este processo, vamos especificar quais são os níveis de referência para os dois fatores controlados no planejamento experimental: Genotype e EnvironmentalStress.

```R
dds <- DESeqDataSetFromTximport(txi.salmon, colData = targets, design = ~EnvironmentalStress)
dds$Genotype <-relevel(dds$Genotype, ref='WT')
dds$EnvironmentalStress <-relevel(dds$EnvironmentalStress, ref='None')
```

Vale a pena discutir um pouco sobre o tipo de objeto criado pelo DESeq com a função DESeqDataSetFromTximport. Este é um objeto do tipo `RangedSummarizedExperiment`, que pos sua vez é uma extensão do `SummarizedExperiment`, como apresentado de forma resumida na figura seguinte.

[![SummarizedExperiment object](images/SummarizedExperiment.png)](https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
)

O objeto `SummarizedExperiment` é uma classe que serve como um recipiente eficiente para armazenar dados experimentais. Ele é utilizado para representar conjuntos de dados tabulares associados a uma matriz de dados principal, onde as linhas geralmente representam recursos biológicos (genes, por exemplo) e as colunas representam amostras. O `SummarizedExperiment` é projetado para armazenar informações adicionais sobre os dados, como metadados sobre as amostras e características biológicas, proporcionando uma estrutura organizada e completa.

A classe `RangedSummarizedExperiment` inclui informações adicionais sobre a localização física dos recursos biológicos em um genoma. Isso é particularmente útil para dados genômicos, como dados de sequenciamento de RNA (RNA-Seq), onde é importante entender a posição dos genes em um genoma.

Os dados das contagens por gene e amostra podem ser acessados com a função `assay`, os nomes dos genes e das amostras com a funcão `dimnames`:

```R
assay(dds)
dimnames(dds)
```

A informação das amostras pode ser acessada com `colData`, e a dos genes com `rowData`:

```R
colData(dds)
rowData(dds)
```

Geralmente, é conveniente excluir genes com níveis muito baixos de expressão. Neste caso, vamos remover os genes nos quais a soma das contagens seja menor que 10.

Essa prática é comumente adotada para filtrar genes com expressão extremamente baixa, os quais podem contribuir pouco para a análise global e introduzir ruído desnecessário nos resultados. Eliminar genes com somas de contagens muito baixas pode melhorar a sensibilidade da análise, focando nos genes mais relevantes e robustos para o experimento em questão.

É importante ressaltar que a escolha do limiar de expressão a ser utilizado pode depender do contexto específico do estudo e da natureza dos dados. Portanto, a definição desse limiar deve ser cuidadosamente considerada em cada análise para garantir resultados mais precisos e interpretáveis.

Vamos a estudar primeiro a distribucao do número de fragmentos por gene por amostra, para isso vamos reformar nosso objeto `assay(dds)` numa forma que o ggplot reconheca facilmente e criaremos entao um boxplot

```R
counts_melt<-melt(assay(dds))
colnames(counts_melt)<-c("Gene","Sample","Count")
counts_melt
ggplot(counts_melt,aes(x=Sample,y=Count)) + 
  geom_boxplot()+
  scale_y_log10()+
  theme_bw() +
  theme(axis.text.x = element_text(angle=60, size = 12, hjust = 1)) +
  xlab("Amostras") +
  ylab("Número de fragmentos sequenciados")
```

Agora, vamos remover os genes com soma de contagens menor que 10. O resultado da funcão `table` nos mostra quantos genes serão removidos.

```R
keep <- rowSums(counts(dds)) >= 10
table(keep)
dim(dds)
dds <- dds[keep,]
dim(dds)
```

Chegou o momento crucial de realizar a análise diferencial de expressão (DESeq), a qual estima os dados para normalizacão, realiza o ajuste dos dados ao modelo e estima os parâmetros essenciais, tais como $` \beta `$ (coeficientes de regressão) e $` \log_2FC `$ (logaritmo na base 2 do fold change).

```R
dds <- DESeq(dds)
```

Por favor, reveja os resultados da função `rowData`:

```R
rowData(dds)
```

Vamos agora transformar os dados de expressão usando o método de estabilização de variância (vst), a fim de conduzir análises exploratórias. 

A transformação de estabilização de variância (vst) é um método comumente aplicado em análises de dados de expressão gênica, especialmente em dados provenientes de sequenciamento de RNA (RNA-seq). Seu objetivo principal é reduzir a heteroscedasticidade, ou seja, a variabilidade não constante ao longo da faixa de intensidades de expressão gênica.

Ao aplicar a transformação vst, busca-se tornar a distribuição das variações de expressão mais homogênea, o que pode ser crucial para garantir a validade estatística de análises subsequentes. Isso é particularmente relevante em estudos de expressão gênica, nos quais as variações nas contagens de RNA podem ser desproporcionalmente grandes em regiões de baixa expressão, tornando difícil a comparação e interpretação correta.

Portanto, a transformação vst desempenha um papel fundamental na preparação dos dados, melhorando a capacidade de detecção de padrões biologicamente relevantes e contribuindo para resultados mais robustos em análises exploratórias e inferenciais.

```R
vsd <- vst(dds, blind=FALSE)
```

Vamos conduzir uma análise de componentes principais com o objetivo de compreender o comportamento dos dados. Para isso, utilizaremos a função plotPCA do pacote DESeq2, que seleciona automaticamente os 500 genes mais variáveis e os utiliza como base para a análise.

A análise de componentes principais (PCA) é uma técnica exploratória poderosa que nos permite visualizar a estrutura e a variabilidade nos dados. Ao selecionar os genes mais variáveis, podemos capturar as principais fontes de variação e examinar como as amostras se diferenciam em relação a esses genes.

A interpretação dos resultados da PCA nos ajudará a identificar padrões, agrupamentos ou possíveis outliers nos dados. Esse passo é fundamental para a compreensão inicial da estrutura dos dados antes de prosseguirmos com análises mais específicas de expressão gênica diferencial. Esperamos observar uma tendência de agrupamento nas amostras de acordo com o planejamento experimental, ou seja, que a separação nos dois primeiros componentes principais seja congruente com o fator experimental em questão.

```R
plotPCA(vsd, intgroup=c("Genotype"),ntop=500) + theme_bw()
plotPCA(vsd, intgroup=c("EnvironmentalStress"),ntop=500) + theme_bw()
plotPCA(vsd, intgroup=c("EnvironmentalStress", "Genotype"),ntop=500) + theme_bw()
```

Desta análise, torna-se evidente que o principal fator que influencia os níveis de expressão é a condição ambiental. O genótipo não exerce uma influência significativa na expressão dos genes; na nossa matriz de design, ele atua mais como um fator de lotes ou bloqueio. Observe que esses dois primeiros componentes explicam 82% da variância observada nos dados.

![PCA](images/transcriptomicsPCA_DESeq2.png)

Outra análise exploratória importante é a análise de correlação entre as amostras. Novamente, aqui esperamos que amostras sob as mesmas condições apresentem valores de correlação mais elevados. Note que, nesta análise, ao contrário do `plotPCA`, não estamos utilizando os 500 genes mais variáveis, mas sim todos os genes.

```R
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$EnvironmentalStress, vsd$Genotype, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$EnvironmentalStress, vsd$Genotype, sep="-")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
```

Agora, vamos realizar os testes estatísticos para selecionar os genes que apresentam mudanças significativas em expressão em resposta à condição ambiental. Para isso, utilizaremos a função results do DESeq2. Essa função calcula os valores de _p-valor_ e _fold change_ para cada gene, utilizando o modelo estatístico que descrevemos anteriormente. O resultado é um objeto do tipo `data.frame` com os valores de _p-valor_ e _fold change_ para cada gene.

No DESeq2, é possível conduzir explicitamente o teste de hipótese para avaliar se o $` \log_2FC `$ é significativamente diferente, maior ou menor que um limiar escolhido pelo pesquisador. Neste caso, usaremos 2 como limiar, indicando que a expressão do gene muda em pelo menos 4 vezes nas condições comparadas.

```R
res_SALT_vs_Control <- results(dds, lfcThreshold=1, altHypothesis="greaterAbs", contrast = c('EnvironmentalStress','NaCl','None'), alpha=0.05)
res_SALT_vs_Control
res_ABA_vs_Control <- results(dds, lfcThreshold=1, altHypothesis="greaterAbs", contrast = c('EnvironmentalStress','ABA','None'), alpha=0.05)
res_ABA_vs_Control
res_Drought_vs_Control <- results(dds, lfcThreshold=1, altHypothesis="greaterAbs", contrast = c('EnvironmentalStress','Drought','None'), alpha=0.05)
res_Drought_vs_Control
```

Vamos agora explorar esses resultados por meio de dois tipos comuns de gráficos: o MA plot e o Volcano plot.

O _MA plot_ é uma representação gráfica que exibe a relação entre a média de abundância (M) e a diferença (A) entre as condições experimentais. É uma ferramenta útil para visualizar a variabilidade dos dados e identificar genes diferencialmente expressos.

O _Volcano plot_, por outro lado, destaca os genes que são estatisticamente significativos em termos de p-valor e fold change. Ele oferece uma visão mais abrangente das diferenças de expressão, permitindo a identificação rápida dos genes mais relevantes.

Vamos explorar essas representações gráficas para obter insights adicionais sobre os genes que respondem significativamente às condições ambientais. Primeiro os _MA plots_

```R
drawLines <- function() abline(h=c(-1,1),col="dodgerblue",lwd=2)
plotMA(res_SALT_vs_Control, main = "SALT_vs_Control"); drawLines()
plotMA(res_ABA_vs_Control, main = "ABA_vs_Control"); drawLines()
plotMA(res_Drought_vs_Control, main = "Drought_vs_Control"); drawLines()
```

Agora, os _Volcano plots_. O DESeq2 não possui uma função pronta para criar esse tipo de gráfico, mas podemos facilmente criá-lo utilizando as tabelas de resultados e o ggplot. Em seguida, faremos a comparação _NaCl vs Control_. Deixo para você realizar as outras duas comparações.

```R
df_res_SALT_vs_Control<-as.data.frame(res_SALT_vs_Control)
df_res_SALT_vs_Control$diffExpressed <- "NO"
df_res_SALT_vs_Control$diffExpressed[df_res_SALT_vs_Control$padj < 0.05 & df_res_SALT_vs_Control$log2FoldChange >0] <- "UP"
df_res_SALT_vs_Control$diffExpressed[df_res_SALT_vs_Control$padj < 0.05 & df_res_SALT_vs_Control$log2FoldChange <0] <- "DOWN"
table(df_res_SALT_vs_Control$diffExpressed)
ggplot(df_res_SALT_vs_Control, aes(x=log2FoldChange, y=-log10(padj), col=diffExpressed)) + 
  theme_bw() +
  geom_point() +
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed')+
  scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), 
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  ggtitle('SALT_vs_Control DEGs')
```
Nesse último gráfico, o eixo y está muito condensado. Vamos aplicar uma transformação em log10 para obter uma melhor visualização:

```R
ggplot(df_res_SALT_vs_Control, aes(x=log2FoldChange, y=-log10(padj), col=diffExpressed)) + 
  theme_bw() +
  geom_point() +
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed')+
  scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), 
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  scale_y_log10()+
  ggtitle('SALT_vs_Control DEGs')
```

Agora, para os outros contrastes, gere o grafico usando como template o codigo anterior:

```R
df_res_ABA_vs_Control<-as.data.frame(res_ABA_vs_Control)
df_res_ABA_vs_Control$diffExpressed <- "NO"
df_res_ABA_vs_Control$diffExpressed[df_res_ABA_vs_Control$padj < 0.05 & df_res_ABA_vs_Control$log2FoldChange >0] <- "UP"
df_res_ABA_vs_Control$diffExpressed[df_res_ABA_vs_Control$padj < 0.05 & df_res_ABA_vs_Control$log2FoldChange <0] <- "DOWN"
table(df_res_ABA_vs_Control$diffExpressed)

df_res_Drought_vs_Control<-as.data.frame(res_Drought_vs_Control)
df_res_Drought_vs_Control$diffExpressed <- "NO"
df_res_Drought_vs_Control$diffExpressed[df_res_Drought_vs_Control$padj < 0.05 & df_res_Drought_vs_Control$log2FoldChange >0] <- "UP"
df_res_Drought_vs_Control$diffExpressed[df_res_Drought_vs_Control$padj < 0.05 & df_res_Drought_vs_Control$log2FoldChange <0] <- "DOWN"
table(df_res_Drought_vs_Control$diffExpressed)
```

Agora, vamos criar um mapa de calor exclusivamente com os genes diferencialmente expressos e os valores de expressão, gerados pela função vst.

Esta visualização nos permitirá explorar padrões de expressão gênica em amostras específicas, concentrando-se nos genes que mostraram alterações significativas em resposta às condições ambientais. A transformação VST (variance stabilizing transformation) aplicada aos dados de expressão ajuda a estabilizar a variação e fornece uma representação mais precisa das diferenças de expressão entre as condições experimentais.

Ao observar o mapa de calor, será possível identificar claramente como a expressão desses genes específicos varia em diferentes contextos ambientais, contribuindo para uma compreensão mais profunda dos padrões de regulação gênica no experimento.

```R
pheatmap(assay(vsd)[row.names(assay(vsd)) %in%
                      row.names(df_res_SALT_vs_Control[which(df_res_SALT_vs_Control$diffExpressed %in% 
                                                               c('UP','DOWN')),]),
                    row.names(targets[which(targets$EnvironmentalStress %in% c("NaCl","None")),])],
         scale='row', 
         annotation_col = targets,
         main = "SALT_vs_Control DEGs")

pheatmap(assay(vsd)[row.names(assay(vsd)) %in%
                      row.names(df_res_ABA_vs_Control[which(df_res_ABA_vs_Control$diffExpressed %in% 
                                                               c('UP','DOWN')),]),
                    row.names(targets[which(targets$EnvironmentalStress %in% c("ABA","None")),])],
         scale='row', 
         annotation_col = targets,
         main = "ABA_vs_Control DEGs")

pheatmap(assay(vsd)[row.names(assay(vsd)) %in%
                      row.names(df_res_Drought_vs_Control[which(df_res_Drought_vs_Control$diffExpressed %in% 
                                                               c('UP','DOWN')),]),
                    row.names(targets[which(targets$EnvironmentalStress %in% c("Drought","None")),])],
         scale='row', 
         annotation_col = targets,
         main = "Drought_vs_Control DEGs")
```

Consulte a ajuda do pacote pheatmap para não mostrar os nomes dos genes no gráfico. Você pode usar o comando ?pheatmap.

Vamos salvar os dados com os identificadores dos genes diferencialmente expressos para cada comparação.

```R
write.table(df_res_SALT_vs_Control[which(df_res_SALT_vs_Control$diffExpressed %in% c('UP','DOWN')),],
            file = "SALT_vs_Control_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(df_res_ABA_vs_Control[which(df_res_ABA_vs_Control$diffExpressed %in% c('UP','DOWN')),],
            file = "ABA_vs_Control_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(df_res_Drought_vs_Control[which(df_res_Drought_vs_Control$diffExpressed %in% c('UP','DOWN')),],
            file = "Drought_vs_Control_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

```

#### Diagramas de Venn dos genes differencialmente expressos

Quando temos vários contrastes, é comum querer saber quais genes são comuns entre as diferentes comparações. Para isso, podemos utilizar os diagramas de Venn. Vamos empregar o pacote [ggVennDiagram](https://cran.r-project.org/web/packages/ggVennDiagram/index.html) para criar esses diagramas utilizando o ggplot2.

Primeiramente, vamos comparar todos os genes diferencialmente expressos entre os diversos contrastes.

```R
library(ggVennDiagram)
diffExpGenes<-list(Drought_vs_Control=rownames(df_res_Drought_vs_Control[which(df_res_Drought_vs_Control$diffExpressed != 'NO'),]),
              ABA_vs_Control    =rownames(df_res_ABA_vs_Control[which(df_res_ABA_vs_Control$diffExpressed != 'NO'),]),
              SALT_vs_Control   =rownames(df_res_SALT_vs_Control[which(df_res_SALT_vs_Control$diffExpressed != 'NO'),]))

ggVennDiagram(diffExpGenes) + scale_fill_gradient(low="blue",high = "red")

process_region_data(Venn(diffExpGenes))
```

Agora, faça as comparações apenas para os genes UP e DOWN por separado.

#### Identificando funções sobre-representadas nos grupos de genes diferencialmente expressos.

Outra pergunta comum é se existe alguma função ou funções que aparecem mais do que esperado por acaso nos conjuntos de genes diferencialmente expressos.

Existem várias estratégias para avaliar esse aspecto, e um pacote que implementa diversas abordagens é o topGO. Neste contexto, utilizaremos o [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html) para identificar termos GO que surgem como sobre-representados. Especificamente, aplicaremos a prova exata de Fisher. No entanto, é importante observar que o pacote topGO oferece outras provas que podem ser de interesse para abordar diversas questões biológicas. Essa diversidade de métodos estatísticos proporciona uma análise abrangente e adaptável às diferentes nuances biológicas que podem surgir durante a avaliação dos termos GO sobre-representados nos conjuntos de genes diferencialmente expressos.

[A prova exata de Fisher](https://en.wikipedia.org/wiki/Fisher%27s_exact_test) é um teste estatístico que avalia a associação entre duas variáveis categóricas, neste caso, a associação entre a presença de um gene em um conjunto diferencialmente expresso e sua associação com termos GO. Essa abordagem estatística nos permite determinar se a ocorrência desses termos nos conjuntos de genes diferencialmente expressos é significativamente maior do que o esperado ao acaso.

Primeiro, é necessário baixar as associações entre os identificadores dos genes de _A. thaliana_ e os termos GO. Você pode realizar o download dessas associações a partir [daqui](https://arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/gene_association.tair.gz). No seu terminal, descomprima o arquivo utilizando o comando `gunzip` e inspecione o conteúdo do arquivo com o comando `less`. Desse arquivo, precisamos apenas da coluna 2 (GeneID) e da coluna 5 (GOID). No seu terminal, utilizando bash, vamos filtrar o conteúdo para manter apenas os campos que nos interessam:

```bash
grep AGI_LocusCode gene_association.tair |cut -f 2,5|sort -u > gene_association.tair.simplified
```

Agora, vamos importar esses dados no R e e reformatar essa tabela numa lista de R, formato que o topGO reconhece.

```R
library(topGO)
GTOGO<-read.delim2("gene_association.tair.simplified", header = FALSE)
geneID2GO<- by(GTOGO$V2,
               GTOGO$V1,
               function(x) as.character(x))
```

Para cada análise de sobre-representação, é necessário criar um vetor com 2 valores: o código do gene para todos os genes e o segundo valor, indicando se o gene está ou não no conjunto de genes diferencialmente expressos. Vamos gerar esses vetores para cada contraste.

```R
allGenes<-as.factor(rownames(assay(dds)))
geneListDrought_vs_Control<-factor(as.integer(allGenes %in% diffExpGenes$Drought_vs_Control))
names(geneListDrought_vs_Control)<-allGenes
```

Agora, temos os dados necessários para realizar a primeira análise de sobre-representação de termos GO. Vamos criar um objeto topGO com essas informações. E vamos visualizar um resumo desse objeto.

```R
Drought_vs_Control_topGO<-new("topGOdata",
                              description = "Drought_vs_Control_ALL", ontology = "BP",
                              allGenes = geneListDrought_vs_Control,
                              nodeSize = 4,
                              annot =  annFUN.gene2GO, gene2GO = geneID2GO)

Drought_vs_Control_topGO
```

Agora, finalmente, podemos realizar o teste estatístico para detectar os termos GO que aparecem mais do que esperado por acaso. Nesse contexto, o topGO nos auxiliará na identificação dos termos GO que estão sobre-representados nos conjuntos de genes diferencialmente expressos.

```R
resultFisherDrought_vs_Control_topGO <- runTest(Drought_vs_Control_topGO, algorithm = "classic", statistic = "fisher")
```

O topGO não calcula ou corrige automaticamente os p-values para testes múltiplos. Para realizar essa correção, vamos primeiro criar um dataframe com os resultados e, em seguida, utilizando o método `p.adjust`, corrigiremos os p-values usando a abordagem de _false discovery rate_.

```R
resultFisherDrought_vs_Control_table<-cbind(termStat(Drought_vs_Control_topGO),
                                            score(resultFisherDrought_vs_Control_topGO))

colnames(resultFisherDrought_vs_Control_table)<-c('Annotated','Significant','Expected','pvalue')
resultFisherDrought_vs_Control_table$padj<-p.adjust(resultFisherDrought_vs_Control_table$pvalue, method = 'fdr')
```

Para verificar quantos termos GO foram detectados como sobre-representados, com um FDR < 0.05, você pode realizar a seguinte análise no R:

```R
dim(resultFisherDrought_vs_Control_table[which(resultFisherDrought_vs_Control_table$padj < 0.05),])
```

Podemos visualizar os termos GO sobre-representados na topologia da ontologia GO. Para isso, vamos utilizar a função `showSigOfNodes` do topGO.

```R
showSigOfNodes(Drought_vs_Control_topGO, score(resultFisherDrought_vs_Control_topGO), firstSigNodes = 10, useInfo = 'all')
```

Realize a análise para os outros dois contrastes.

```R
sesionInfo()
```
####  Identificando grupos de genes com perfis de expressão semelhantes.

Primeiro, vamos selecionar todos os identificadores dos genes diferencialmente expressos nos três contrastes. Em seguida, iremos recuperar os valores de expressão desses genes em todas as amostras. Para representar a expressão, utilizaremos o método VST. Repare que é preciso reduzir as replicatas a uma única coluna, já que as replicatas da mesma amostra não contribuem com informação independente do perfil de expressão do gene nessa condição.

```R
allDEGs<-unique(c(rownames(df_res_Drought_vs_Control[which(df_res_Drought_vs_Control$diffExpressed != 'NO'),]),
                  rownames(df_res_ABA_vs_Control[which(df_res_ABA_vs_Control$diffExpressed != 'NO'),]),
                  rownames(df_res_SALT_vs_Control[which(df_res_SALT_vs_Control$diffExpressed != 'NO'),])))	

allDEGsExp<-assay(vsd)[allDEGs,]

allDEGsExpMeans<-as.data.frame(matrix(data=NA,ncol=4,nrow=nrow(allDEGsExp)))
rownames(allDEGsExpMeans)<-rownames(allDEGsExp)
colnames(allDEGsExpMeans)<-c('ControlMean','ABAMean','SaltMean','DroughtMean')
head(allDEGsExpMeans)
allDEGsExpMeans$ControlMean<-rowMeans(allDEGsExp[,rownames(targets[which(targets$EnvironmentalStress == 'None'),])])
allDEGsExpMeans$ABAMean<-rowMeans(allDEGsExp[,rownames(targets[which(targets$EnvironmentalStress == 'ABA'),])])
allDEGsExpMeans$SaltMean<-rowMeans(allDEGsExp[,rownames(targets[which(targets$EnvironmentalStress == 'NaCl'),])])
allDEGsExpMeans$DroughtMean<-rowMeans(allDEGsExp[,rownames(targets[which(targets$EnvironmentalStress == 'Drought'),])])

head(allDEGsExpMeans)
dim(allDEGsExpMeans)         

```

Vamos a transformar os valores de expressao usando o método [Z-score](https://en.wikipedia.org/wiki/Standard_score), para focar mais na forma do perfil de expressão do que nos valores absolutos de expressão. Para isso vamos criar uma funcao nova em R.

Para identificar os grupos de genes com padrões (perfil) de expressão semelhante, utilizaremos o pacote [Mclust](https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html). O Mclust é uma ferramenta de análise de misturas que implementa a modelagem de mistura gaussiana para agrupar dados em subpopulações ou clusters com base na distribuição normal multivariada, e geralmente identifica clusters robustos e coesos. Essa análise pode levar alguns minutos.

```R
library(mclust)
clusters<-Mclust(allDEGsExpMeans,G=2:90)
max(clusters$classification)
table(clusters$classification)
```

Agora vamos a visualizar os perfis de expressão dos genes em cada grupo. Para isso, vamos criar um objeto do tipo `data.frame` com os valores de expressão e os grupos de cada gene.

```R
ClusterSel<-4
ClusterSelExp<-allDEGsExpMeansZscore[which(rownames(allDEGsExpMeansZscore)
                                           %in% 
                                             names(clusters$classification[clusters$classification == ClusterSel])),]


head(ClusterSelExp)
conditionMeans<-colMeans(ClusterSelExp)
ClusterSelExp<-as.data.frame(ClusterSelExp)
ClusterSelExp$Gene<-rownames(ClusterSelExp)
head(ClusterSelExp)
ClusterSelExpMelt<-melt(ClusterSelExp,id.vars = 'Gene')
head(ClusterSelExpMelt)
tittle=paste("Group ",ClusterSel, sep='')
ggplot(ClusterSelExpMelt, aes(y=value,x=variable,group=Gene))+
  geom_line(colour='gray')+
  theme_bw()+
  ylab('Relative expression value')+
  xlab('Condition')+
  ggtitle(tittle)
```

Agora, vamos criar um mapa de calor com os valores de expressão dos genes em cada grupo. Para isso, vamos utilizar a função `pheatmap` do pacote pheatmap.

```R
pheatmap(assay(vsd)[which(rownames(assay(vsd)) %in% ClusterSelExp$Gen),], 
         scale = 'row', 
         annotation_col = targets,
         main= tittle)
```

Ou só as médias:

```R
pheatmap(allDEGsExpMeans[which(rownames(allDEGsExpMeans) %in% ClusterSelExp$Gene),], 
         scale = 'row',
         cluster_cols = FALSE,
         main= tittle)
```