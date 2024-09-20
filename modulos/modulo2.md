# Módulo 2 - Trabalho com dados genômicos

## Visualização de dados genômicos

### Arquivos de anotação

Os formatos [GFF (General Feature Format) e GTF (Gene Transfer Format)](https://www.ensembl.org/info/website/upload/gff.html) são amplamente utilizados para anotar genomas. Estes arquivos descrevem as localizações de diferentes elementos genômicos, como genes, exons, introns e outras características, usando uma estrutura tabular que pode ser visualizada no terminal Linux. Vmmos explorar cada formato com comandos do terminal.

#### Estrutura Básica dos Arquivos GFF e GTF

Ambos os formatos compartilham a ideia de dividir a anotação em colunas. As colunas comuns em ambos são:

- **Seqname:** Nome da sequência (cromossomo ou scaffolds) onde a característica foi encontrada.
- **Fonte (Source):** Ferramenta ou banco de dados que produziu a anotação.
- **Tipo (Feature):** Tipo da característica anotada, como "gene", "mRNA", "exon", etc.
- **Início (Start):** Primeira posição do genoma onde a característica está localizada.
- **Fim (End):** Última posição da característica.
- **Score:** Valor de confiança (ou . se não for aplicável).
- **Strand:** Cadeia de DNA (+ para sentido direto, - para sentido reverso).
- **Frame:** Informações de fase de leitura do codon (usado em CDSs).
- **Atributos (Attributes):** Informações adicionais no formato chave-valor.

Diferenças:
- O GFF (versão mais comum é GFF3) é mais simples e genérico. O campo de atributos usa formato separado por ponto-e-vírgula (;) com valores nomeados.
- O GTF, mais específico para genes, contém atributos mais detalhados e foi projetado para ser usado pelo Ensembl. Os atributos são pares chave-valor separados por espaços e terminados por ponto-e-vírgula.

Exemplo GFF:

```
scaffold_1  maker gene  1000  2000  .  +  .  ID=gene0001;Name=my_gene
scaffold_1  maker mRNA  1000  2000  .  +  .  ID=mRNA0001;Parent=gene0001;Name=my_mRNA
scaffold_1  maker exon  1000  1500  .  +  .  ID=exon0001;Parent=mRNA0001

```

Exemplo GTF:

```
scaffold_1  maker gene  1000  2000  .  +  .  gene_id "gene0001"; gene_name "my_gene";
scaffold_1  maker transcript 1000  2000  .  +  .  gene_id "gene0001"; transcript_id "mRNA0001";
scaffold_1  maker exon  1000  1500  .  +  .  gene_id "gene0001"; transcript_id "mRNA0001";

```

#### Anotação do genoma de _Spodoptera frugiperda_

No módulo anterior, você descarregou o genoma de Spodoptera frugiperda, que ficou armazenado na pasta `~/dia1`. Por favor, mova esse arquivo para a pasta `~/dia2`, entre nessa pasta e descarregue a anotação do mesmo genoma em formato GFF usando o seguinte comando:

```bash
cd ~/dia2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/101/765/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.gff.gz
gunzip GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.gff.gz
```

Visualize o conteúdo do arquivo descarregado usando o comando `less` e identifique cada uma das colunas mencionadas acima.

Quais valores podem aparecer na coluna de **Features**? Podemos usar o Linux para obter uma lista dos valores únicos que aparecem na coluna Feature do arquivo. Note que essa é a terceira coluna no arquivo:

```bash
cat GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.gff.gz|grep -v "#"| cut -f 3|sort -u
```

Você também pode contar quantas vezes cada um desses valores aparece no arquivo de anotação:

```bash
cat GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.gff.gz|grep -v "#"| cut -f 3|sort | uniq -c
```

Repare no número de features do tipo exon e CDS que aparecem no arquivo. Por que esses números são diferentes?

Vamos extrair uma regiao de interesse do arquivo GFF
### Arquivos de mapeamento de leituras
