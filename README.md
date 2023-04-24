## Introduction

For now, ChloroAnno has six functions and could be split into three categories based on annotaion stages.

Before annotaion: Isomerism, Curate

During annotation: Tidy, Correct 

After  annotaion: Check, Convert 

### Isomerism
If you assemble the plastome with IRs, the assembly result should be two equimolar isomeric sequences. Both of them are right and coexist in the plant (Palmer 1983; Walker et al. 2015). For publication, people tend to use just one of them. You can pick the configuration of which the SSC region is in the same direction as most of your data or outgroup, which makes downstream annotation/analysis easier. see [GetOrganelle wiki](https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#what-are-the-differences-of-complete1path_sequencefasta-and-complete2path_sequencefasta-output-files) and [Weiwen Wang & Robert Lanfear 2019](https://academic.oup.com/gbe/article/11/12/3372/5637229) for more details.

### Curate
This function was written for curate the raw sequence before annotation
1. The strand direction. Although none of public database specicfies a specific
   strand, it always assumes the minus strand for trnH as the default direction.
2. The start position. Since the chloroplast is the circular genome. The start position
     in fasta is defined by people. It always assumes the trnH gene should at the head to 
     avoid any gene locates on the gap.
3. Whether there were amubiguous nucleotides in sequences

### Tidy
This function is written for GeSeq. Its gff results is not in line with [Sequence Ontology definetion for GFF3](http://www.sequenceontology.org/resources/gff3.html) and thus need a tidy.

### Correct

This is the most important function in ChloroAnno. `correct` task could infer information from two annotation 
files and correct one of them using another one.

A question here is why it is necessary to combine information from different software for annotation? Could I just use one annotation software result? First of all, the annotation from a single software may not be accurate, for example, there may be errors in start or stop codons. There are many reasons for such errors, mainly including the inaccuracy of the reference and the characteristics of the software. For example, MAKER often miss too short segments (<30 bp). Based on our experience, the annotation results of PGA are often more accurate when the quality of the reference sequence is relatively high. Can we just use PGA for annotation directly? The answer is yes, but the workload is large when it comes to annotating chloroplasts in large quantities. When annotating the chloroplast genome, if only the basic or model plants (such as _Amborella trichopoda_ or _Arabidopsis thaliana_) are used as references, some genes that appear in specific lineages may not be annotated. On the other hand, if we choose related species for annotation, we may face the problem of poor quality of reference, such as confusion in gene names. If the reference genome mistakenly writes atpI as aptI or abbreviates trnK-UUU as trnK, similar errors will also appear in the annotation results. Therefore, when using PGA for annotation, it is necessary to correct each reference chloroplast genome, which requires a lot of workload. Therefore, our strategy is to use the basic or model plants as references and PGA for annotation; use the reference dataset in CPGAVAS2 and CPGAVAS2 for annotation. Finally, combine the results of both to generate an annotation file. **It should be noted that the annotation of rps12 gene could only be inferred from PGA results for now due to its particularity.**

Of course, users could replace CPGAVAS2 with GeSeq. We choose CPGAVAS2 beacause it provides an docker image that enable use to run a batch of annotation tasks locally. According to our evaluation, the results of the GeSeq and CPGAVAS2 are very similar. But, GeSeq always generate some error tRNA annotaions.

~~~text
G2023021001     blatN   gene    135684  136697  69.5    -       1       ID=gene-blatn_trnI-GAU_1;gbkey=gene;gene=trnI-GAU;gene_biotype=tRNA;Note=blatN_hit_trnI-GAU_NC_035998.1%2C_position_1_-_72%2C_psl_score_69.5%2C_coverage_100.00%25%2C_match_100.00%25
G2023021001     blatN   tRNA    135684  135718  .       -       1       ID=trna-blatn_trnI-GAU_1;gbkey=tRNA;gene=trnI-GAU
G2023021001     blatN   tRNA    136661  136697  .       -       1       ID=trna-blatn_trnI-GAU_1;gbkey=tRNA;gene=trnI-GAU
G2023021001     blatN   exon    135684  135718  .       -       1       ID=blatn_trnI-GAU_1_exon_2;Parent=gene-blatn_trnI-GAU_1;gbkey=exon;gene=trnI-GAU
G2023021001     blatN   intron  135719  136660  .       -       1       ID=blatn_trnI-GAU_1_intron_1;Parent=gene-blatn_trnI-GAU_1;gbkey=intron;gene=trnI-GAU
G2023021001     blatN   exon    136661  136697  .       -       1       ID=blatn_trnI-GAU_1_exon_1;Parent=gene-blatn_trnI-GAU_1;gbkey=exon;gene=trnI-GAU

~~~

After correct, ChloroAnno will run a check automaticly. See the next section for check detail.

### Check
While annotation files are generated, it does not imply that the annotated result is correct; therefore, annotation validation is necessary. At present, automatic check achieves the following functions:
1. Starting codon check (for RNA editing phenomena in the starting codon, as it is impossible to introduce RNAseq data for verification during the assembly process, they are all marked as pseudo processing)
2. Is there a legal stop codon 
3. Is there a stop codon inside the sequence
4. Is the sequence length a multiple of 3
5. Is the CDS too short (less than 33bp)
6. Is there any overlap in genomic regions

### Convert
The format convertion is one of the eternal themes of bioinformatics. For chloroplast genome annotation, there are two major formats: gff and GenBank. However, even for the same format, there are some small variations across files from different sources. For now, ChloroAnno could recognise annotaion files from following sources:

| Source |  Name | Format|
|  :---: | :---: | :---: |
|Database|NCBI   |  gff  |
|Database|NCBI   |GenBank|
|Database|GWH    |  gff  |
|Database|GenBase|GenBank|
|Database|GenBank|  tbl* |
|Database|GenBase|  tbl* |
|Software|GeSeq  |  gff  |
|Software|PGA    |GenBank|
|Software|CPGAVAS|GenBank|

*Note:* tbl could only be used as output files.

## Manual
### meta file
1. The ChloroAnno use a flexible but unified meta file informat. Users need to choose the appropriate
header based on functions they want to use. Here is an overview of all illegal headers. For 
specific headers used in each function, please refer to the following sections.

- `inpath1`  : Input file path. **Every meta file should include this column.**
- `inpath2`  : Some functions need two input files. This is for the second input file.
- `refpath`  : Some functions need a reference file.
- `organism` : The organism name.
- `prefix`   : The prefix character used for naming gene. e.g. 'AAA' in 'AAA001'
- `informat1`: The input file format (e.g. genbank, fasta). ChloroAnno could determine the file type based on file
suffix (e.g. '.fasta'). So, this column is unnecessary.
- `informat2`: he input file format for inpath2

Since ChloroAnno use this flexible input format, **users need to include the header name in each meta file**. However, one benefit of specifying headers in meta file is that there is no fixed order for the columns in meta file.

2. All meta file should be tab delimited.

3. Although a relative path should work, it is recommended to use an absolute path for safety.

### Isomerism

Select the path based on a reference fasta. This function will determine the quartered structure of chloroplast and then
construct a new genomic sequence based on the SSC direction on reference chloroplast.

#### meta file

Required  headers: inpath1, refpath

Optionanl headers: informat1

~~~shell
inpath1 refpath
/path/to/assembly1.fasta    /path/to/ref.fasta
~~~

#### command

~~~shell
python ChloroAnno -t iso -i meta.info -o output_directory
~~~
*Note:* The `outfmt` option is not applicable to this function. All output would be written to fasta format.

### Curate

Only need input fasta files

#### meta file

Required headers: inpath1

~~~shell
inpath1
/path/to/assembly1.fasta
~~~

#### command

~~~shell
python ChloroAnno -t curate -i meta.info -o /path/to/output/directory
~~~

### Tidy

This function is used to tidy GeSeq GFF files and is not needed if GeSeq is not used to annotate your chloroplast.

#### meta file
~~~shell
inpath1
/path/to/GeSeq.gff
~~~
#### command
~~~shell
python ChloroAnno -t tidy -i meta.info --out /path/to/output/directory
~~~
*Note:* The output format is gff, no need to specify.


### Correct
#### meta file

Required  headers: inpath1, inpath2, prefix

Optionanl headers: organism

*Note:* 
1. The more reliable annotaion result (e.g. PGA) should be placed under `inpath1` and the more informative result (e.g. GeSeq, CPGAVAS2) should be placed under `inpath2`.
2. The organism name is needed if you want to output a formal GenBank file, otherwise, the default value is 'plant'.
 
~~~shell
inpath1 inpath2 prefix  organism
/path/to/PGA.gb /path/to/CPGAVAS2.gb    AAA Arabidopsis thaliana
~~~

#### command
~~~shell
python ChloroAnno -t correct -i meta.info --outfmt tbl --out tbl_dir
~~~

### Check
This function only checks the input annotation file and does not ouput any file.

#### meta file

Required  headers: inpath1

Optionanl headers: inpath2. If you want to check the gff, then you need to provide a fasta file under inpath2.

~~~shell
# Case 1: Check GenBank file
inpath1
/path/to/assembly1.gb
# Case 2: Check gff file
inpath1 inpath2
/path/to/assembly1.gff  /path/to/assembly1.fasta
~~~

#### command
~~~shell
python ChloroAnno -t check -i meta.info
~~~

### Convert
This function convert the annotation file from one format to another. The accepted file format includes:

*Note:* 
1. The rps12 gene will be ignored for now due to its particularity. I hope a future updation could fix this.
2. The tbl format only supports as an output format, not as an input format.

#### meta file

Required  headers: inpath1, inpath2, refpath

Optionanl headers: informat1, informat2

~~~shell
inpath1
/path/to/assembly1.gb
~~~

#### command
~~~shell
python ChloroAnno -t convert -i meta.info --outfmt gff -o output_directory
~~~