![Curso de Extensão: Introdução à Genômica e Bioinformática](Figs/FondoCENA.jpg)

# Curso de Extensão: Introdução à Genômica e Bioinformática

## Introdução

### Objetivo 

O curso busca capacitar os alunos com conhecimentos sobre bancos de dados biológicos, métodos de sequenciamento de ácidos nucléicos, alinhamento de sequências, análise de dados de sequenciamento em larga escala, montagem e anotação de genomas, pangenômica, genômica comparativa, e análise de expressão diferencial de genes. Através de uma combinação de aulas teóricas baseadas em leituras de artigos científicos e componentes práticos usando dados reais e ferramentas bioinformáticas, os alunos desenvolverão habilidades essenciais para realizar pesquisas e análises bioinformáticas utilizando seus próprios dispositivos. O curso também enfatiza a importância do uso do sistema operacional Linux e de ferramentas de linha de comando, preparando os alunos para enfrentar os desafios práticos da área. 

### Generalidades 

O curso está dividido em módulos, cada um com um componente teórico e um componente prático. Cada módulo é realizado em dois encontros, um por semana, com duração de 3 horas cada. No primeiro encontro do módulo, abordamos o componente teórico, enquanto no segundo encontro, focamos no componente prático.

O componente teórico de cada módulo é construído a partir de leituras de livros ou artigos científicos que os alunos devem realizar previamente ao encontro com o professor. Durante o encontro, a discussão será centrada nesses artigos, guiada por uma série de perguntas que serão disponibilizadas antecipadamente.

O material para o componente prático será disponibilizado em este repositório. Os primeiros dois módulos cobrem os aspectos gerais do curso. A partir do terceiro módulo, aprofundaremos na parte metodológica da Genômica, com um enfoque computacional.

O componente prático consiste em trabalhar com dados disponíveis publicamente nos computadores. Para isso, os computadores a serem usados serão os dos mesmos alunos, seguindo o conceito de “[Bring Your Own Device](https://en.wikipedia.org/wiki/Bring_your_own_device)”, que permitirá que os alunos ainda fiquem familiarizados com os processos de instalação de software de bioinformática, além de ficar com o software nas suas máquinas e revisitar em qualquer momento os componentes práticos do curso. O único requerimento é que esteja instalado o sistema operacional Linux Ubuntu 22.04 (ou superior) com ambiente gráfico. Para conseguir desenvolver todas as aulas práticas com facilidade será necessário que o computador tenha pelo menos 16GB de RAM, 8 cores de processamento e 200GB de espaço em disco. Se os computadores não tiverem essa configuração mínima não é garantido que as sessões práticas serão executas com sucesso. 

Antes do início do curso, é imperativo que os alunos treinem por conta própria e se familiarizem com o uso do sistema operacional Linux, especialmente com a linha de comandos usando a shell Bash. Para isso, será disponibilizado material de apoio que os alunos deverão estudar de forma completamente autônoma. Esse é um pré-requisito essencial para o desenvolvimento dos temas do curso. Um segundo pré-requisito é que os alunos cheguem para o primeiro encontro com todo o software instalado, para isso será disponibilizado um [script em linguagem bash](setting_env.sh) que realizará o processo de instalação de todo o software de forma semi-automática, caso tiver problemas os alunos terão um canal de acesso a discussão pelo Discord, que será mantido durante a duração do curso para resolver dúvidas pontoais de forma assíncrona. 

## Módulos

### Módulo 0 - Linux

### Módulo 1 - Bancos de dados biológicos 

[Neste módulo](modulo1.md), os alunos explorarão a evolução e a definição da bioinformática, examinando artigos que traçam a história da área e discutem seu impacto tanto globalmente quanto no Brasil. Serão apresentados aos conceitos de bancos de dados primários e secundários, com um foco especial nos recursos do National Center for Biotechnology Information e na cobertura de sequências de proteínas pelo AlphaFold Protein Structure Database. Os alunos também investigarão a anotação de genomas procarióticos na era dos metagenomas e como esses avanços impulsionam a ciência brasileira. O componente prático incluirá o uso de ferramentas básicas de Linux para operações com arquivos de texto, busca de dados em bancos de dados do NCBI, trabalho com sequências nos formatos Fasta e GenBank, e utilização do EMBOSS para recuperar sequências de grandes bancos de dados e extrair regiões específicas de um genoma. As discussões em classe serão guiadas por perguntas sobre a definição, origem, e principais avanços da bioinformática, assim como a identificação dos principais bancos de dados na área.

### Módulo 2 - Métodos de Sequenciamento de Ácidos Nucléicos

[Neste módulo](modulo2.md), os alunos explorarão os métodos de sequenciamento de ácidos nucléicos, abordando desde a complexidade dos genomas eucarióticos até a coevolução e estrutura dos genomas procarióticos. Discutiremos a importância da metagenômica e da genômica de célula única na reconstrução de genomas de procariontes não cultivados e as limitações desses métodos. Os alunos aprenderão sobre o Projeto Earth BioGenome (EBP), seus objetivos e sua contribuição para a conservação da biodiversidade global e a compreensão da evolução biológica. Também será abordado o uso das sequências genômicas como material de tipo para descrições taxonômicas de procariontes, bem como os desafios e benefícios dessa abordagem. O módulo incluirá discussões sobre como o sequenciamento genômico pode unificar diferentes áreas da pesquisa biológica, integrando dados com outras abordagens biológicas. Analisaremos o "Nanopore Adaptive Sampling" e suas vantagens para a detecção de espécies de baixa abundância em amostras metagenômicas, comparando tecnologias de sequenciamento de terceira geração. Por fim, exploraremos as inovações químicas que permitiram o desenvolvimento do sequenciamento de próxima geração (NGS) e como essas inovações impactaram a precisão e a eficiência do NGS. No componente prático, os alunos calcularão a quantidade de dados necessária para sequenciar genomas de diferentes tamanhos e tecnologias, analisarão dados de sequenciamento e trabalharão com genomas completos, anotações e arquivos de mapeamento de leituras 
