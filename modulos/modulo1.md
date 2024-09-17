# Módulo 1 - Bancos de dados biológicos 

## Testando nosso ambiente

Antes de iniciar a sessão prática, vamos conferir se o software que iremos usar está instalado corretamente.

Execute o seguinte comando:

```bash
apt list libcurl4-openssl-dev zlib1g-dev libbz2-dev build-essential libudunits2-dev libgdal-dev gdal-bin libfreetype6-dev libfontconfig1-dev software-properties-common dirmngr libblas-dev liblapack-dev libatlas-base-dev gfortran zlib1g-dev libcurl4-openssl-dev libxml2-dev git r-base | grep -P -w 'amd64|all'
```

Estamos testando a instalação de 19 pacotes a partir do repositório do sistema operacional. Esses pacotes são de propósito geral, não especificamente de bioinformática, mas são necessários (são dependências) para a instalação de alguns pacotes de bioinformática. Na tela, deve aparecer a lista dos 19 pacotes com a palavra '[installed]'. Se algum pacote não aparecer, consulte o seu instrutor para verificar como ele pode ser instalado.

A maioria dos pacotes de bioinformática foi instalada com o sistema de gerenciamento [Conda](https://conda.org/). Para testá-los, vamos entrar em cada um dos ambientes Conda instalados e verificar se o principal comando do pacote é encontrado. Se aparecer na tela 'comando não encontrado', consulte o seu instrutor para instalar o pacote.

```bash
conda activate emboss
which wossname
which seqret
conda deactivate
```

```bash
conda activate blast
which makeblastdb
conda deactivate
```

```bash
conda activate flye
which flye
conda deactivate
```

```bash
conda activate hifiasm
which hifiasm
conda deactivate
```
```bash
conda activate bandage
which Bandage
conda deactivate
```
```bash
conda activate quast
which quast
conda deactivate
```
```bash
conda activate fastqc
which fastqc
conda deactivate
```
```bash
conda activate bbmap
which bbduk.sh
conda deactivate
```

```bash
conda activate sratoolkit
which fastq-dump
conda deactivate
```
```bash
conda activate tidk
which tidk
conda deactivate
```
```bash
conda activate compleasm
which compleasm
conda deactivate
```
```bash
conda activate merqury
which merqury.sh
which meryl
conda deactivate
```
```bash
conda activate jupiterplot
which circos
which samtools
which minimap2
conda deactivate
```
```bash
conda activate igv
which igv
which java
conda deactivate
```
```bash
conda activate eggnogmapper
which emapper.py
conda deactivate
```

```bash
conda activate transcriptomics
which ffq
which fastqc
which bbkuk.sh
which multiqc
conda deactivate
```

```bash
conda activate genomescope2
which genomescope2
conda deactivate
```

```bash
conda activate redotable
which redotable
which java
conda deactivate
```

Em algumas das sessões práticas, vamos usar R com algumas bibliotecas específicas. Vamos testar se essas bibliotecas estão disponíveis.

Inicie uma sessão de R e execute os seguintes comandos, um a um:

```R
library(ggplot2)
library(tximport)
library(DESeq2)
library(topGO)
library(Rgraphviz)
library(pheatmap)
library(mclust)
library(reshape2)
library(readr)
library(ggVennDiagram)
```

## Buscas em bancos de dados web

O [NCBI](https://www.ncbi.nlm.nih.gov/), uma agência dos EUA, reúne alguns dos bancos de dados mais utilizados nas biociências. Vamos explorar esses recursos para desenvolver habilidades avançadas de busca.

Vamos procurar todos os artigos publicados pelo autor "Weisshaar, B" em 2009. Para isso, utilize os modificadores de busca DP (Data de Publicação) e AU (Autor) para restringir a busca aos campos apropriados. Você pode encontrar mais informações sobre os campos de busca disponíveis [neste link](https://pubmed.ncbi.nlm.nih.gov/help/#using-search-field-tags). Junto com seu instrutor, identifique os diferentes componentes das entradas dos dados recuperados.

Explore as opções de busca avançada.

## Operações básicas em Bioinformática

### Ferramentas do Unix úteis na bioinformática.

Após adquirir alguma familiaridade com os fundamentos do [sistema operacional Linux](unix.md), vamos explorar como alguns de seus comandos mais básicos podem ser extremamente úteis na área de bioinformática. Você entenderá por que o Linux é o sistema operacional de escolha na bioinformática.

Para realizar esses exercícios, você precisa usar o arquivo [file1.tar.gz](files/file1.tar.gz). Após baixá-lo, o arquivo deve estar na sua pasta "Downloads". Você deve descompactá-lo em seu diretório HOME.

```
cd
mv ~/Downloads/file1.tar.gz ~/
tar xvzf file1.tar.gz
```

#### Algumas operações básicas com arquivos.

Usando alguns comandos do UNIX, podemos obter informações sobre arquivos e o conteúdo deles de forma rápida e eficiente, muitas vezes sem a necessidade de abrir o arquivo, que pode ser muito grande, para obter essas informações.

No subdiretório "~/dia1/", encontre o arquivo "TAIR10_pep_20101214_updated.fasta.gz", que corresponde à base de dados de sequências de proteínas previstas no genoma da planta modelo _Arabidopsis thaliana_. Para saber quantas linhas este arquivo possui, descomprioma o arquivo e conte o número de linhas com os comandos:

```
cd 
cd ~/dia1
gunzip TAIR10_pep_20101214_updated.fasta.gz
wc -l TAIR10_pep_20101214_updated.fasta
```
Pode conhecer o tamnho do arquivo com o comando _ls_:

```
ls -l -h TAIR10_pep_20101214_updated.fasta
```

O que faz a opção -h no comando 'ls'? Consulte a página de manual do _ls_ para saber.

Na maioria das vezes, é importante visualizar o conteúdo do arquivo, seja no início ou no final.  No entanto, devido ao grande tamanho dos arquivos com os quais normalmente se trabalha, não é conveniente abrir o arquivo com nenhum editor de texto, pois isso pode reduzir o tempo de resposta do computador. Podemos visualizar as primeiras ou ultimas linhas de um arquivo de texto com os comandos  _head_ e _tail_ respetivamente

```
head TAIR10_pep_20101214_updated.fasta
tail TAIR10_pep_20101214_updated.fasta
```
Esses comandos mostram as primeirais ou ultimas 10 linhas do arquivo. O que você pode fazer para mostrar um número maior de linhas? Consulte a página de manual do  _head_ 

Repare na saída do comando `head`. Está mostrando um registro de sequência no formato 'fasta'. Este formato é o mais simples para armazenar sequências, tanto de ácidos nucleicos quanto de proteínas. Sua estrutura é muito simples. Cada registro começa com uma linha que tem no seu início o sinal _>_ seguido de uma cadeia de caracteres de comprimento arbitrário que funciona como o identificador da sequência. Em seguida, nas linhas subsequentes, aparece a sequência em si, em quantas linhas forem necessárias.

Pode usar o comando _grep_ para localizar todas as linhas que tem um padrão de texto específico, ou seja, uma cadeia de texto específica. Vamos identificar todas a linas que comecam com o sinal _>_

```
grep ">" TAIR10_pep_20101214_updated.fasta
```

São muitas linhas, vamos usar o _pipe_ para examinar só as primeiras 4 linas que tem o padrão de texto procurado, assim:

```
grep ">" TAIR10_pep_20101214_updated.fasta | head -n 4
```

#### Formatos de sequências

Existem diferentes formatos para sequências, geralmente em texto simples. Isso significa que elas podem ser visualizadas e editadas com qualquer editor de texto, como vi ou pico. Alguns desses formatos são mais comuns do que outros, e muitos programas de bioinformática aceitam vários dos formatos mais comuns ([Leonard et al., 2007](https://pubmed.ncbi.nlm.nih.gov/18428774/)).

Todos os formatos de sequências têm uma característica (campo) em comum: um identificador para cada sequência, para que esta possa ser reconhecida de forma única.

##### Fasta

O formato mais simples é conhecido como Fasta. Nele, uma entrada, que é uma sequência, é dividida em duas partes: a linha de identificação, que deve começar com o símbolo ">" seguido imediatamente pelo identificador da sequência, que pode ser qualquer cadeia de caracteres sem espaços. As linhas imediatamente após o identificador correspondem à própria sequência.

O formato Fasta é o formato de sequências mais amplamente utilizado em aplicações de bioinformática.

```
>gi|110742030|dbj|BAE98952.1| putative NAC domain protein [Arabidopsis thaliana]
MEDQVGFGFRPNDEELVGHYLRNKIEGNTSRDVEVAISEVNICSYDPWNLRFQSKYKSRDAMWYFFSRRE
NNKGNRQSRTTVSGKWKLTGESVEVKDQWGFCSEGFRGKIGHKRVLAFLDGRYPDKTKSDWVIHEFHYDL
LPEHQRTYVICRLEYKGDDADILSAYAIDPTPAFVPNMTSSAGSVVNQSRQRNSGSYNTYSEYDSANHGQ
QFNENSNIMQQQPLQGSFNPLLEYDFANHGGQWLSDYIDLQQQVPYLAPYENESEMIWKHVIEENFEFLV
DERTSMQQHYSDHRPKKPVSGVLPDDSSDTETGSMIFEDTSSSTDSVGSSDEPGHTRIDDIPSLNIIEPL
HNYKAQEQPKQQSKEKVISSQKSECEWKMAEDSIKIPPSTNTVKQSWIVLENAQWNYLKNMIIGVLLFIS
VISWIILVG
```

##### GenBank

O formato GenBank é utilizado pelo 'National Center for Biotechnology Information' ([NCBI](https://www.ncbi.nlm.nih.gov/)), o maior repositório de sequências de ácidos nucleicos e proteínas do mundo. O NCBI, juntamente com o [EMBL-EBI](https://www.ebi.ac.uk/) e o [DDBJ](https://www.ddbj.nig.ac.jp/), mantém conjuntamente o 'The International Nucleotide Sequence Database' ([Mizrachi, 2008](https://pubmed.ncbi.nlm.nih.gov/27896718/)).

Uma entrada neste formato é composta por duas partes. A primeira parte abrange as posições de 1 a 10 e geralmente contém o nome do campo, como LOCUS, DEFINITION, ACCESSION ou SOURCE. A segunda parte de cada entrada contém informações correspondentes ao campo em questão. Cada entrada é finalizada com o símbolo '//'. Você pode encontrar mais informações sobre esse tipo de arquivo [aqui](http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html).

```
LOCUS       BAE98952                 429 aa            linear   PLN 27-JUL-2006
DEFINITION  putative NAC domain protein [Arabidopsis thaliana].
ACCESSION   BAE98952
VERSION     BAE98952.1
DBSOURCE    accession AK226863.1
KEYWORDS    .
SOURCE      Arabidopsis thaliana (thale cress)
  ORGANISM  Arabidopsis thaliana
            Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
            Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae;
            Pentapetalae; rosids; malvids; Brassicales; Brassicaceae;
            Camelineae; Arabidopsis.
REFERENCE   1
  AUTHORS   Totoki,Y., Seki,M., Ishida,J., Nakajima,M., Enju,A., Morosawa,T.,
            Kamiya,A., Narusaka,M., Shin-i,T., Nakagawa,M., Sakamoto,N.,
            Oishi,K., Kohara,Y., Kobayashi,M., Toyoda,A., Sakaki,Y.,
            Sakurai,T., Iida,K., Akiyama,K., Satou,M., Toyoda,T., Konagaya,A.,
            Carninci,P., Kawai,J., Hayashizaki,Y. and Shinozaki,K.
  TITLE     Large-scale analysis of RIKEN Arabidopsis full-length (RAFL) cDNAs
  JOURNAL   Unpublished
REFERENCE   2  (residues 1 to 429)
  AUTHORS   Totoki,Y., Seki,M., Ishida,J., Nakajima,M., Enju,A., Morosawa,T.,
            Kamiya,A., Narusaka,M., Shin-i,T., Nakagawa,M., Sakamoto,N.,
            Oishi,K., Kohara,Y., Kobayashi,M., Toyoda,A., Sakaki,Y.,
            Sakurai,T., Iida,K., Akiyama,K., Satou,M., Toyoda,T., Konagaya,A.,
            Carninci,P., Kawai,J., Hayashizaki,Y. and Shinozaki,K.
  TITLE     Direct Submission
  JOURNAL   Submitted (26-JUL-2006) Motoaki Seki, RIKEN Plant Science Center;
            1-7-22 Suehiro-cho, Tsurumi-ku, Yokohama, Kanagawa 230-0045, Japan
            (E-mail:mseki@psc.riken.jp, URL:http://rarge.gsc.riken.jp/,
            Tel:81-45-503-9625, Fax:81-45-503-9586)
COMMENT     An Arabidopsis full-length cDNA library was constructed essentially
            as reported previously (Seki et al. (1998) Plant J. 15:707-720;
            Seki et al. (2002) Science 296:141-145).
            This clone is in a modified pBluescript vector.
            Please visit our web site (http://rarge.gsc.riken.jp/) for further
            details.
FEATURES             Location/Qualifiers
     source          1..429
                     /organism="Arabidopsis thaliana"
                     /db_xref="taxon:3702"
                     /chromosome="1"
                     /clone="RAFL08-19-M04"
                     /ecotype="Columbia"
                     /note="common name: thale cress"
     Protein         1..429
                     /product="putative NAC domain protein"
     Region          5..137
                     /region_name="NAM"
                     /note="No apical meristem (NAM) protein; pfam02365"
                     /db_xref="CDD:426740"
     CDS             1..429
                     /gene="At1g01010"
                     /coded_by="AK226863.1:89..1378"
ORIGIN      
        1 medqvgfgfr pndeelvghy lrnkiegnts rdvevaisev nicsydpwnl rfqskyksrd
       61 amwyffsrre nnkgnrqsrt tvsgkwkltg esvevkdqwg fcsegfrgki ghkrvlafld
      121 grypdktksd wvihefhydl lpehqrtyvi crleykgdda dilsayaidp tpafvpnmts
      181 sagsvvnqsr qrnsgsynty seydsanhgq qfnensnimq qqplqgsfnp lleydfanhg
      241 gqwlsdyidl qqqvpylapy enesemiwkh vieenfeflv dertsmqqhy sdhrpkkpvs
      301 gvlpddssdt etgsmifedt ssstdsvgss depghtridd ipslniiepl hnykaqeqpk
      361 qqskekviss qksecewkma edsikippst ntvkqswivl enaqwnylkn miigvllfis
      421 viswiilvg
//
```

#### Algumas operações básicas com sequências no formato Fasta

Durante o restante desta seção e na maior parte do que resta do curso, usaremos apenas sequências no formato Fasta. Por favor, verifique se as sequências de _A. thaliana_ no arquivo 'TAIR10_pep_20101214_updated.fasta' estão neste formato. Você pode usar o comando `head nome_do_arquivo` ou o comando `less nome_do_arquivo` para fazer isso.

Você já teve que contar o número de sequências num arquivo ou alterar o identificador das sequências no formato Fasta? Se fossem apenas algumas sequências, isso poderia ser feito facilmente em qualquer editor de texto. No entanto, quando se trata de milhares de sequências, a opção de usar um editor de texto deixa de ser viável. Felizmente, alguns comandos do Unix nos permitem realizar essas tarefas simples de forma rápida.

Como observado anteriormente, o comando `grep` pode nos ajudar a contar o número de sequências em um arquivo Fasta. O modificador `-c` conta o número de linhas que contêm um padrão específico em um arquivo, e podemos aproveitar o fato de que em um arquivo Fasta o símbolo _>_ aparece uma única vez para cada sequência. Para contar o número de sequencias armazenadas no arquivo pode usar o comando:

```
grep -c ">" TAIR10_pep_20101214_updated.fasta
```

Em outras situações, é importante modificar o identificador de cada sequência para incluir, por exemplo, uma abreviatura que represente o nome da espécie à qual a sequência pertence. Novamente, o Unix nos permite fazer essa alteração rapidamente usando o comando `sed`. Vamos adicionar a partícula 'ATH' a cada um dos identificadores do arquivo, aproveitando o fato de que à esquerda de cada identificador temos o símbolo '>'. Observe que os resultados desta operação estão sendo armazenados em um novo arquivo:

```
sed 's/>/>ATH_/' TAIR10_pep_20101214_updated.fasta > TAIR10_pep_20101214_updated.mod.fasta
```
As linhas com os identificadores neste arquivo são muito extensas, para muitos programas isto não é desejável. Vamos eliminar tudo que aparecer depois do primeiro ' |'. Para isso vamos usar [regular expressions](https://www.gnu.org/software/sed/manual/html_node/Regular-Expressions.html).

```
sed -r 's/ | .*$//' TAIR10_pep_20101214_updated.mod.fasta > TAIR10_pep_20101214_updated.mod2.fasta
```

### The European Molecular Biology Open Software Souite - EMBOSS

Alguns dos exercícios desta seção são baseados no [Tutorial do Usuário do EMBOSS](http://emboss.open-bio.org/html/use/ch04.html).

[EMBOSS](http://emboss.open-bio.org/) é um pacote de software de análise gratuito, desenvolvido em C, voltado para as necessidades da comunidade de biologia molecular e bioinformática. O software lida automaticamente com dados em vários formatos e permite a obtenção de sequências diretamente da web. Além disso, o pacote inclui bibliotecas extensas, oferecendo uma plataforma para que outros cientistas possam desenvolver e compartilhar softwares, seguindo o espírito de código aberto. EMBOSS também integra diversos pacotes e ferramentas de análise de sequências em uma solução unificada.

Todos os programas do EMBOSS são executados a partir da linha de comando do Unix. A utilidade do EMBOSS, chamada _wossname_, gera uma lista de todas as aplicações do EMBOSS.

Primeiro, você deverá ativar o ambiente Conda que tem o EMBOSS instalado:

```bash
conda activate emboss
```

Repare que seu prompt mudou; agora aparece _(emboss)_, indicando que este é o nome do ambiente atualmente ativo.

Digite _wossname_ no prompt do Linux e pressione a tecla Enter. Em seguida, digite a palavra _protein_ para listar todos os programas do EMBOSS que lidam de alguma forma com proteínas. Os programas do EMBOSS começam com uma breve descrição em uma linha e, depois, solicitam informações. Neste caso, você verá:

```bash
wossname
Find programs by keywords in their short description
Text to search for, or blank to list all programs: protein
SEARCH FOR 'PROTEIN'
antigenic      Find antigenic sites in proteins
backtranambig  Back-translate a protein sequence to ambiguous nucleotide sequence
backtranseq    Back-translate a protein sequence to a nucleotide sequence
charge         Draw a protein charge plot
checktrans     Report STOP codons and ORF statistics of a protein
compseq        Calculate the composition of unique words in sequences
emowse         Search protein sequences by digest fragment molecular weight
epestfind      Find PEST motifs as potential proteolytic cleavage sites
freak          Generate residue/base frequency table or plot
fuzzpro        Search for patterns in protein sequences
fuzztran       Search for patterns in protein sequences (translated)
garnier        Predict protein secondary structure using GOR method
helixturnhelix Identify nucleic acid-binding motifs in protein sequences
hmoment        Calculate and plot hydrophobic moment for protein sequence(s)
iep            Calculate the isoelectric point of proteins
makeprotseq    Create random protein sequences
maskambigprot  Mask all ambiguity characters in protein sequences with X
msbar          Mutate a sequence
mwcontam       Find weights common to multiple molecular weights files
mwfilter       Filter noisy data from molecular weights file
octanol        Draw a White-Wimley protein hydropathy plot
oddcomp        Identify proteins with specified sequence word composition
patmatdb       Search protein sequences with a sequence motif
patmatmotifs   Scan a protein sequence with motifs from the PROSITE database
pepcoil        Predict coiled coil regions in protein sequences
pepdigest      Report on protein proteolytic enzyme or reagent cleavage sites
pepinfo        Plot amino acid properties of a protein sequence in parallel
pepnet         Draw a helical net for a protein sequence
pepstats       Calculate statistics of protein properties
pepwheel       Draw a helical wheel diagram for a protein sequence
pepwindow      Draw a hydropathy plot for a protein sequence
pepwindowall   Draw Kyte-Doolittle hydropathy plot for a protein alignment
preg           Regular expression search of protein sequence(s)
profit         Scan one or more sequences with a simple frequency matrix
prophecy       Create frequency matrix or profile from a multiple alignment
prophet        Scan one or more sequences with a Gribskov or Henikoff profile
pscan          Scan protein sequence(s) with fingerprints from the PRINTS database
psiphi         Calculates phi and psi torsion angles from protein coordinates
showpep        Display protein sequences with features in pretty format
shuffleseq     Shuffle a set of sequences maintaining composition
sigcleave      Report on signal cleavage sites in a protein sequence
tcode          Identify protein-coding regions using Fickett TESTCODE statistic
tmap           Predict and plot transmembrane segments in protein sequences
tranalign      Generate an alignment of nucleic coding regions from aligned proteins
wordcount      Count and extract unique words in molecular sequence(s)
```

#### Recuperando Sequências de Bancos de Dados

Os programas do EMBOSS podem ler sequências de vários bancos de dados, desde que a sequência seja referenciada no formato banco_de_dados:identificador. Este formato é um exemplo de Uniform Sequence Address [(Endereço Uniforme de Sequência), ou USA](http://emboss.open-bio.org/html/use/ch06s06.html).

_seqret_ lê uma sequência e a escreve de volta. É provavelmente o programa mais utilizado do EMBOSS. Para extrair a sequência da proteína correspondente ao identificador **AT3G48060.1** do arquivo *TAIR10_pep_20101214_updated.fasta*, você pode usar o seguinte comando com o seqret do EMBOSS:

```bash
seqret -sequence TAIR10_pep_20101214_updated.fasta:AT3G48060.1 -outseq AT3G48060.1.fasta
```
Este comando faz o seguinte:

- sequence: indica o arquivo de entrada e o identificador específico da sequência de interesse.
- outseq: define o arquivo de saída, onde a sequência extraída será salva (neste caso, AT3G48060.1.fasta).

Após executar o comando, a sequência da proteína será extraída e salva no arquivo especificado.

Usando ferramentas do Linux, identifique em qual linha do arquivo original está presente o identificador **AT3G48060.1**.

Com o seqret, podemos reformatar sequências, por exemplo, convertendo uma sequência do formato FASTA para o formato GenBank. Você pode executar o seguinte comando:

```bash
seqret AT3G48060.1.fasta genbank::AT3G48060.1.gb
```

No entanto, lembre-se de que o seqret não cria informações adicionais. Compare a entrada no formato FASTA com o arquivo gerado em formato GenBank, e observe que o arquivo GenBank resultante conterá apenas as informações presentes no arquivo FASTA original.

Programas do EMBOSS também podem lidar com múltiplas sequências. Vamos usar o programa infoseq para obter informações sobre as proteínas de Arabidopsis thaliana que estão presentes na mitocôndria.

```bash
infoseq TAIR10_pep_20101214_updated.fasta:ATMG*
```

Repare que podemos limitar quais campos são exibidos no resultado final. Vamos manter só o identificador da sequência

```bash
infoseq -only -name -length TAIR10_pep_20101214_updated.fasta:ATMG*
```
De forma semelhante, você pode recuperar várias sequências a partir de um arquivo multifasta. Outra alternativa, quando você tem uma lista de identificadores para serem recuperados, é:

```bash
seqret @selected.ids selected_seqs.fasta
```

Como você poderia verificar se há o mesmo número de identificadores nos arquivos _selected.ids_ e *selected_seqs.fasta*?


