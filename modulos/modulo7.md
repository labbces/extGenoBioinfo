# Módulo 7 - Pangenômica e Genômica Comparativa

## Calculando Average Nucleotide Identity (ANI) entre genomas

O fastANI é uma ferramenta utilizada para análise de similaridade entre genomas bacterianos. Ele calcula o ANI, um índice de similaridade que determina se duas amostras genômicas pertencem à mesma espécie com base em um corte de 95% de similaridade. O ANI é útil para estudos de diversidade e filogenia bacteriana, sendo um método preferível ao uso do gene 16S rRNA por oferecer uma resolução mais alta ([Jain et al., 2018](https://www.nature.com/articles/s41467-018-07641-9)). 

Aqui vamos usar o software [fastANI](https://github.com/ParBLiSS/FastANI), que implementa uma análise livre de alinhamento, precisa e rápida para comparar genomas de Bacillus spp.

Primeiramente, vamos instalar o fastANI, dentro de uma nova pasta chamada dia7 no seu $HOME:

```bash
mkdir ~/dia7
cd ~/dia7
wget https://github.com/ParBLiSS/FastANI/releases/download/v1.34/fastANI-linux64-v1.34.zip
unzip fastANI-linux64-v1.34.zip
rm fastANI-linux64-v1.34.zip
./fastANI -h
```

O comando anterior deve exibir a mensagem de ajuda do fastANI. Agora, vamos baixar os genomas de Bacillus spp. que serão comparados. Você pode facilmente descarregar os genomas de Bacillus spp. do [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/), mas para hoje, vamos usar os genomas que estão disponíveis no nosso servidor, vamos descarregar eles no seu computador e conferir que a descarga foi feita corretamente:

```bash
wget https://labbces.cena.usp.br/shared/extGenoBioinfo/dia7/bacillus_MAGs.tar.gz
wget https://labbces.cena.usp.br/shared/extGenoBioinfo/dia7/bacillus_MAGs.tar.gz.md5
md5sum -c bacillus_MAGs.tar.gz.md5
wget https://labbces.cena.usp.br/shared/extGenoBioinfo/dia7/bacillus_TypeMaterial.tar.gz
wget https://labbces.cena.usp.br/shared/extGenoBioinfo/dia7/bacillus_TypeMaterial.tar.gz.md5
md5sum -c bacillus_TypeMaterial.tar.gz.md5
```

Na sua tele deve aparecer a mensagem "bacillus_MAGs.tar.gz: OK" e "bacillus_TypeMaterial.tar.gz: OK". Agora, vamos descompactar os arquivos baixados:

```bash
tar -xzf bacillus_MAGs.tar.gz
tar -xzf bacillus_TypeMaterial.tar.gz
```
Isso vai criar as pastas bacillus_MAGs e bacillus_TypeMaterial, que contêm os genomas de Bacillus spp. que vamos comparar. Na pasta TypeMaterial, há apenas genomas de cepas de Bacillus spp. usadas para descrever essas espécies. Já na pasta MAGS, estão os genomas de Bacillus spp. montados a partir de dados metagenômicos. Os arquivos com a extensão fofn contêm uma lista dos caminhos onde estão os arquivos dos genomas, que podemos usar com o fastANI para realizar uma comparação muitos contra muitos. Agora, vamos calcular o ANI entre os genomas de Bacillus spp. usando o fastANI.

```bash
./fastANI --ql bacillus_MAGs.fofn --rl bacillus_TypeMaterial.fofn -o Bacillus_spp.fastANI_MAGs_vs_TypeMaterial.tsv -t 20
```	

Peça ajuda novamente ao comando fastANI para entender o que cada opção faz. O comando anterior calculará o ANI entre os genomas de Bacillus spp. presentes nas pastas MAGs e TypeMaterial. O resultado será salvo no arquivo Bacillus_spp.fastANI_MAGs_vs_TypeMaterial.tsv. O parâmetro -t 20 indica que o fastANI usará 20 threads para realizar o cálculo; ajuste esse parâmetro de acordo com a quantidade de núcleos disponíveis no seu computador.

Após a execução do comando, você pode visualizar o resultado com o comando less:

```bash
less Bacillus_spp.fastANI_MAGs_vs_TypeMaterial.tsv
```
Na saida teremos um arquivo de texto separado por tabulações com os resultados do ANI entre os genomas de Bacillus spp. Na primeira coluna temos o nome do genoma de Bacillus spp. da pasta MAGs (query), na segunda coluna temos o nome do genoma de Bacillus spp. da pasta TypeMaterial (reference), na terceira coluna temos o ANI entre esses genomas, na quarta coluna a contagem de fragmentos mapeados bidirecionalmente e na quinta coluna o total de fragmentos de consulta. A fração de alinhamento (em relação ao genoma de consulta (query)) é simplesmente a razão entre os mapeamentos e o total de fragmentos.

Vamos a reformatar um pouco esse aquivo para que seja mais facil de ser processado por outros pacotes:

```bash
sed -E 's#bacillus_MAGs/data/GC[A-Z]_[0-9.]+/##g' Bacillus_spp.fastANI_MAGs_vs_TypeMaterial.tsv|sed -E 's#bacillus_TypeMaterial/data/GC[A-Z]_[0-9.]+/##'| sed -E 's/_genomic.fna//g'> Bacillus_spp.fastANI_MAGs_vs_TypeMaterial.formated.tsv
```

Visualize o novo arquivo criado e explique as mudanças que foram realizadas.

Vamos tentar visualizar esses dados de uma forma mais clara. Para isso, vamos usar o R, com o rstudio. Pode executar o comando rstudio no seu terminal e abrir um novo script. Vamos usar o seguinte script para visualizar os dados:

```R
library(ggplot2)
rm(list=ls())
threshold <- 90
setwd("~/dia7")
dataANI<-read.delim("Bacillus_spp.fastANI_MAGs_vs_TypeMaterial.formated.tsv", header = FALSE, dec='.')
dataANI$V6<-dataANI$V4*100/dataANI$V5
colnames(dataANI)<-c('query','reference','ANI','AlignedFragmentes','TotalFragments', 'FractionAligned')
head(dataANI)
ggplot(dataANI, aes(x=query,y=reference))+
  geom_tile(aes(fill = ANI), colour = "white") +
  scale_fill_gradient(low = "red", high = "steelblue")

library(reshape2)
df_wide <- dcast(dataANI, query ~ reference, value.var = "ANI")
df_wide[is.na(df_wide)] <- 0
rownames(df_wide)<-df_wide$query
rownamesANI<-df_wide$query
df_wide$query<-NULL
head(df_wide)
df_wide <- as.matrix(sapply(df_wide, as.numeric))
rownames(df_wide)<-rownamesANI
display_matrix <- ifelse(df_wide > threshold, round(as.numeric(df_wide), 2), "")
display_matrix <- matrix(display_matrix, nrow = nrow(df_wide), ncol = ncol(df_wide))

# Ensure row and column names are retained in display_matrix
rownames(display_matrix) <- rownames(df_wide)
colnames(display_matrix) <- colnames(df_wide)
dim(display_matrix)
head(df_wide)

color_gradient <- colorRampPalette(c("green", "white", "red"))(100)
library(pheatmap)
pheatmap(df_wide, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = display_matrix,
         fontsize_row = 4,
         fontsize_number = 4, 
         fontsize_col = 4,
         na_col = "blue",
         color = color_gradient)

```

Vamos tentar adicionar os nomes das espécies ao gráfico. Podemos obter os nomes das especies a partir do arquivo multifasta do genoma que usamos para calcualr o ANI. Vamos usar o seguinte comando para obter os nomes das espécies:

```bash
cd ~/dia7/
for i in bacillus_TypeMaterial/*/*/*_genomic.fna; do ASSEMBLY=`basename $i`; ASSEMBLY=${ASSEMBLY/_genomic.fna};LINE=`grep  ">" $i|head -n1|cut -f2,3 -d' '`; echo $ASSEMBLY $LINE;done > bacillus_TypeMaterial.species.txt
```

Depois, voltando no R, vamos adicionar os nomes das espécies ao gráfico:

```R
speciesNames<-read.delim("bacillus_TypeMaterial.species.txt",header=FALSE)
head(speciesNames)
rownames(speciesNames)<-speciesNames$V1
speciesNames$V1<-NULL
colnames(speciesNames)<-'Species'

pheatmap(df_wide, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = display_matrix,
         fontsize_row = 4,
         fontsize_number = 4, 
         fontsize_col = 4,
         annotation_col = speciesNames,
         color = color_gradient)
```

Another alternative is to use the ComplexHeatMap package, which allows for more customization of the heatmap. The following code generates a heatmap using the ComplexHeatMap package:

```R
library(ComplexHeatmap)

species_filtered <- speciesNames[rownames(speciesNames) %in% colnames(df_wide), , drop = FALSE]
head(species_filtered)
extended_colors <- colorRampPalette(c("blue", "green", "red", "purple", "orange"))(length(unique(species_filtered$Species)))

# Example data and heatmap with color and texture annotations
Heatmap(df_wide,
        name = "ANI",
        row_names_gp = gpar(fontsize = 6), 
        column_names_gp = gpar(fontsize = 6),
        top_annotation = HeatmapAnnotation(
          col = list(column_variable = extended_colors),  # Extended color palette
          df = species_filtered
        ),
        column_split = factor(species_filtered$Species),
        column_title_gp = gpar(fontsize = 8, fontface = "bold"))

```

Discuta com seu professor os resultados. Você consegue identificar quais genomas de Bacillus spp. são mais semelhantes entre si? E quais são mais diferentes? O que isso pode nos dizer sobre a diversidade de Bacillus spp.?

## Trabalhando com genes ortólogos em OMA

O [Orthologous MAtrix (OMA)](https://omabrowser.org/oma/home/) é um banco de dados de genes ortólogos que permite a comparação de genomas de diferentes espécies. O OMA é uma ferramenta útil para a identificação de genes ortólogos, que são genes que divergiram a partir de um ancestral comum. A identificação de genes ortólogos é importante para a inferência de funções gênicas e para a comparação de genomas.

Vamos realizar alguns exercícios no OMA browser que foram apresentados originalmente [neste artigo](https://pmc.ncbi.nlm.nih.gov/articles/PMC7014581/). 

### Gene UniProtKB S100P_HUMAN. 

Cada gene (também conhecido como entrada) no OMA possui um identificador OMA, que consiste no código de espécie de cinco letras do UniProtKB e um número único de 5 dígitos. Por exemplo, o gene da proteína P de ligação ao cálcio S100 do ser humano possui o identificador OMA HUMAN74803. É possível buscar um gene usando um identificador, uma sequência de proteína ou uma busca em texto completo no navegador OMA:

**Busca por um identificador.** Busque um ID de gene no navegador digitando ou colando no campo de busca na página inicial ou na barra superior de qualquer página. No OMA, uma busca com o ID do UniProtKB (S100P_HUMAN), o número de acesso primário (P25815), RefSeq (NP_005971) ou EMBL (AAH06819) recupera a entrada HUMAN74803. Note que o OMA não permite buscas por números de acesso secundários do UniProtKB ou IDs Unigene. Repare que o identificar OMA (e.g., HUMAN74803 não é estavel entre vesões, especialmente para as espeices que são atualizadas frequentemente).

**Busca por uma sequência de proteína**. Copie e cole a sequência da proteína no campo de busca e escolha Protein sequence no menu suspenso. Neste exemplo, a sequência "MTELETAMGMIIDVFSRYSGSEGSTQTLTKGELKVLMEKELPGFLQSGKDKDAVDKLLKDLDANGDAQVDFSEFIVFVAAITSACHKYFEKAGLK" recupera a mesma entrada, HUMAN74803. Correspondências exatas de sequência de outras espécies também são exibidas. Espaços e números de linha na sequência são ignorados. Recomenda-se buscar por sequência de proteína para evitar ambiguidades. Caso a entrada esperada não seja recuperada, é possível usar a função Approximate Sequence Search.

**Busca por uma palavra-chave usando a busca em texto completo**. Selecione Full-text search no menu suspenso e digite sua palavra-chave para buscar. Podem ser identificadores alternativos; por exemplo, uma busca com o gene Ensembl (ENSG00000163993), transcrito (ENST00000296370) ou identificadores de proteína (ENSP00000296370), ou identificadores PubMed (PMID:15632002), todos retornam nosso gene original, HUMAN22168. Outros textos na descrição do gene também podem ser usados. Aspas (") podem ser usadas para buscar uma sequência exata de palavras. Por exemplo, uma busca em texto completo por “S100 calcium-binding protein P” também recupera a entrada HUMAN74803.

**Obtenha mais informações sobre seu gene**. Após buscar pelo gene, você será levado à página do gene, que fornece informações externas. Você também pode encontrar essas informações clicando na aba Information. A página de informações inclui o ID OMA, descrição, organismo, locus, outros IDs e referências cruzadas, arquiteturas de domínios e anotações de Gene Ontology.

### Encontrando ortólogos para um gene

#### P53_RAT. 

Após encontrar seu gene no OMA browser, o próximo passo é recuperar os ortólogos. Consulte a seção sobre tipos de ortólogos identificados pelo OMA: https://omabrowser.org/oma/homologs/. A página de informações da entrada no OMA fornece links para os ortólogos (Figura 6). A seguir, descrevemos como encontrar cada um dos diferentes tipos de ortólogos.

#### Encontrando ortólogos pareados

A proteína p53 do rato atua como um supressor de tumores em muitos tipos de tumores; ela induz a parada do crescimento ou a apoptose, dependendo das circunstâncias fisiológicas e do tipo celular. Primeiro, recupere a entrada buscando o identificador “P53_RAT.” Na aba Orthologs, é possível ver que atualmente existem 40 ortólogos pareados. Esse valor pode mudar com a inclusão de novas espécies nas atualizações do OMA. Clicar na aba Orthologs exibe a lista de todos os ortólogos pareados.

Todas as cardinalidades de relacionamento entre o p53 do rato e os ortólogos pareados são exibidas. Vamos considerar um exemplo para cada um dos diferentes tipos de relacionamento:

A entrada para o camundongo (Mus musculus) é listada como tendo uma relação de “1:1 ortólogo” com o gene p53 do rato. Isso indica que o gene p53 do rato tem apenas um ortólogo no camundongo, e esse ortólogo é o único no camundongo.

As duas entradas para o peixe-zebra (*Danio rerio*) têm uma relação 1 para muitos ortólogos (1:n) com o gene p53 do rato. Isso indica que o gene p53 do rato tem mais de um ortólogo *Danio rerio*, mas que ambos os genes ortólogos no peixe têm apenas um gene ortólogo no rato. Isso implica que o gene p53 foi duplicado em um ancestral do peixe, mas após o evento de especiação.

Curiosamente, o gene p53 humano não está listado, o que indica que o gene humano não é um ortólogo pareado do gene p53 do rato. No entanto, tanto os genes p53 humano quanto o do rato são encontrados no grupo hierárquico HOG:0430403, indicando que estão relacionados.

#### Encontrando grupos OMA

Um grupo OMA contém conjuntos de genes que são todos ortólogos entre si dentro do grupo. Isso implica que há, no máximo, uma entrada de cada espécie em um grupo.

Ao clicar na aba OMA Groups na entrada do gene p53 do rato, é retornado o grupo OMA 1388790. O grupo tem o fingerprint: "HKKGEPC". O fingerptint é uma subsequência específica de um determinado grupo na versão atual do OMA.

O grupo OMA HKKGEPC contém 88 entradas de espécies diferentes. Note que as entradas apresentam uma arquitetura de domínios semelhante à dos ortólogos pareados.

Agora realize os exercicios propostos no OMA browser: https://omabrowser.org/oma/academy/module/OMA_browser