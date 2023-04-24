### 不使用PGA注释反向重复区的原因
trnH基因有一种情况是正好在序列头尾gap区，pre_check.py是通过SeqIO读取时产生的warning来捕捉这一信息的。明显不可能有两个基因同时跨gap区，所以捕获的warning信息只可能是一个，但是如果注释反向重复区的话，反向重复区在gff注释中也可能产生跨gap的注释，这时会捕获多个warning，无法判断真实情况。
### PGA和GeSeq转gff使用不同前缀的原因
correct.py脚本进行自动校正的时候是根据gene ID**逐个**进行删除、插入的，如果使用相同前缀的话可能出现前面插入了AAAAA001，后面又删除了AAAAA001的现象。至于为何要逐个插入，而不使用全部先记录再删除、插入的策略，是因为两个gff中的ID不是一一对应的，实际上我是通过基因名匹配的，如果批量记录、删除的话，那么多拷贝基因就很难处理（逐个处理的过程中是通过gene区域判断是多拷贝中的哪个的）
### 被子植物和其他植物的rps12基因
被子植物的rps12是分两个区域，三个CDS，通常双拷贝，且其中只有一个CDS的区域会被用两次;其他植物的是两个区域，两个CDS，其中一个区域会被用两次。目前脚本判断谁用两次是通过基因区域中的CDS数量（一个CDS的用两次。两个CDS的用一次。），但是这没法对其他植物的起作用。后期考虑为其他植物单写一个基于长度判断的method。
### CPGAVAS2无法注释rps12
目前CPGAVA2无法注释rps12，所以仍然需要PGA进行辅助注释
### GenBank文件和gff文件不一致
使用`BCBio.GFF`和`Bio.SeqIO`处理这两种文件的话，可以保证他们在Python内部均为`SeqRecord`类，便于处理。但是NCBI不知道出于什么原因，同一个序列在两种格式下，其属性是不同的。例如基因名，在gff中，是"Name"，在gb中是"gene"；基因id在gff中是"ID",在gb中是"locus_tag"。此外，还有嵌套结构的问题。上述问题造成了很多类需要分开写。对于GenBank来说，最重要的是有无"locus_tag"属性，如果没有这一属性，则gene和子feature的匹配难度直线上升。
## 开发过程中碰到的各种`SeqRecord`类的结构
计划解决方案：弄一个包括所有可能qulifiers键的SeqRecord类，然后全部转成这个类，再写不同地的写入、写出的方法。

GenBank要求CDS是CompoundLocation，gff则要求是SimpleLocation，如何统一这个矛盾？初步设想的解决方案：在CoreSeqFeature中分开存储，因为分开好合并，合并的不好分开

### GFF NCBI
~~~shell
GFF(NCBI)
    - scaffold -> features[Gene]
        - gene -> subfeatures[CDS/tRNA/rRNA]
        - 
GFF(CPGAVAS2)
GenBank(NCBI)
GenBank(PGA)
~~~
### 对临时文件夹的处理，全部移交给main

#TODO
1. 自动识别不同来源GenBank文件rps12基因功能
2. 更统一、但更灵活的输入文件形式