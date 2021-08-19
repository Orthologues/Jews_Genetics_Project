## Lund University BINP50 Master's Thesis Project in bioinformatics: 
## Identifying a fine genetic structure across Ashkenazic Jews

#### Author: Jiawei Zhao (ji8842zh-s@student.lu.se)
#### Supervisor: Prof. Eran Elhaik

## Section I: Basic Setup for the project

### Create a directory called "../main_folder" to conduct our major analysis. However, the directory should be git ignore or it would be too large to be put on Github.

```python
mkdir -p ../main_folder && cd ../main_folder && ls .
```

    /lunarc/nobackup/projects/snic2019-34-3/jzhao/MasterThesis/main_folder


#### Download PLINK 1.9

```python
!curl -sSL https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip -o plink.zip && \
unzip plink.zip
```

```python
!git clone https://github.com/DReichLab/EIG && ls -l EIG/CONVERTF/*
```

#### Download the tool to convert .geno & .snp & .ind formats into PLINK formats: https://github.com/DReichLab/EIG


```python
!git clone https://github.com/DReichLab/EIG && ls -l EIG/CONVERTF/*
```

```python
!cat EIG/CONVERTF/par.ANCESTRYMAP.EIGENSTRAT
```

    genotypename:    example.ancestrymapgeno
    snpname:         example.snp
    indivname:       example.ind
    outputformat:    EIGENSTRAT
    genotypeoutname: example.eigenstratgeno
    snpoutname:      example.snp
    indivoutname:    example.ind



```python
!cat EIG/CONVERTF/par.PED.PACKEDPED
```

    genotypename:    example.ped
    snpname:         example.pedsnp # or example.map, either works
    indivname:       example.pedind # or example.ped, either works
    outputformat:    PACKEDPED
    genotypeoutname: example.bed
    snpoutname:      example.pedsnp
    indivoutname:    example.pedind
    familynames:     NO


#### use 
#### <code>sudo apt-get install -y libgsl-dev -y libopenblas-dev -y liblapacke-dev -y libblas-dev -y libatlas-base-dev</code> 
#### to install necessary dependencies for makefile


```python
!cd EIG/src && make && make install
```

#### We might need to use the 1240k dataset at https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/ as the reference panel in our study. However, we need to reformat the 1240k dataset into .vcf format as input for BEAGLE imputation reference.

```python
cd .. #Head back to the main folder
```

```python
! curl -SL \
  https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V44/V44.3/SHARE/public.dir/v44.3_1240K_public.tar \
  -o v44.3.1240K.tar
```

      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
    100 2805M  100 2805M    0     0  24.9M      0  0:01:52  0:01:52 --:--:-- 24.2M



```python
! tar -xf v44.3.1240K.tar
```


```python
! du -sh v44.3* && rm -rf v44.3.1240K.tar
```

    3,7M	v44.3_1240K_public.anno
    2,7G	v44.3_1240K_public.geno
    372K	v44.3_1240K_public.ind
    75M	v44.3_1240K_public.snp

#### Convert .geno format into PLINK format

```python
!echo "\
genotypename:    v44.3_1240K_public.geno\n\
snpname:         v44.3_1240K_public.snp\n\
indivname:       v44.3_1240K_public.ind\n\
outputformat:    PACKEDPED\n\
genotypeoutname: v44.3_1240K_public.bed\n\
snpoutname:      v44.3_1240K_public.bim\n\
indivoutname:    v44.3_1240K_public.fam\n\
familynames:     NO" > par.1240K_MAP.PACKEDPED
```


```python
!EIG/src/convertf -p par.1240K_MAP.PACKEDPED
```

    parameter file: par.1240K_MAP.PACKEDPED
    genotypename: v44.3_1240K_public.geno
    snpname: v44.3_1240K_public.snp
    indivname: v44.3_1240K_public.ind
    outputformat: PACKEDPED
    genotypeoutname: v44.3_1240K_public.bed
    snpoutname: v44.3_1240K_public.bim
    indivoutname: v44.3_1240K_public.fam
    familynames: NO
    read 1073741824 bytes
    read 2147483648 bytes
    read 2859357147 bytes
    packed geno read OK
    end of inpack
    numvalidind:   9275  maxmiss: 9275001
    packedped output
    ##end of convertf run

#### Rename the 1240K reference panel

```python
! ls v44.3_1240K_public*|while read ref;do \
    newname=$(echo $ref|sed -r 's/_1240K_public/\.1240K/'); \
    mv $ref $newname; \
  done
```

#### Create two directories 
#### <code>panel_1240k_vcf</code> to store the entire VCF file of the v44.3 1240K reference panel
#### <code>panel_1240k_vcf_chr</code> to store VCF files by chromosome of the v44.3 1240K reference panel


```python
! mkdir -p panel_1240k_vcf panel_1240k_vcf_chr
```


```python
! plink --bfile v44.3.1240K --recode vcf --out panel_1240k_vcf/v44.3.1240K
```

    PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
    (C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to panel_1240k_vcf/v44.3.1240K.log.
    Options in effect:
      --bfile v44.3.1240K
      --out panel_1240k_vcf/v44.3.1240K
      --recode vcf
    
    64039 MB RAM detected; reserving 32019 MB for main workspace.
    1233013 variants loaded from .bim file.
    9275 people (5217 males, 3926 females, 132 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to panel_1240k_vcf/v44.3.1240K.nosex .
    9275 phenotype values loaded from .fam.
    Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
    phenotypes to be ignored, use the --allow-no-sex flag.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 9275 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Warning: 1103 het. haploid genotypes present (see
    panel_1240k_vcf/v44.3.1240K.hh ); many commands treat these as missing.
    Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    treat these as missing.
    Total genotyping rate is 0.632602.
    1233013 variants and 9275 people pass filters and QC.
    Among remaining phenotypes, 0 are cases and 9275 are controls.
    Warning: Underscore(s) present in sample IDs.
    --recode vcf to panel_1240k_vcf/v44.3.1240K.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.



```python
! du -sh panel_1240k_vcf/*
```

    32K	panel_1240k_vcf/v44.3.1240K.hh
    512	panel_1240k_vcf/v44.3.1240K.log
    512	panel_1240k_vcf/v44.3.1240K.nosex
    43G	panel_1240k_vcf/v44.3.1240K.vcf


#### Install bcftools and tabix by conda

```python
! conda install -c bioconda/label/cf201901 bcftools -c bioconda/label/cf201901 tabix && bgzip panel_1240k_vcf/v44.3.1240K.vcf
```


```python
! du -sh panel_1240k_vcf/*
```

    32K	panel_1240k_vcf/v44.3.1240K.hh
    512	panel_1240k_vcf/v44.3.1240K.log
    512	panel_1240k_vcf/v44.3.1240K.nosex
    3.0G	panel_1240k_vcf/v44.3.1240K.vcf.gz



```python
! tabix panel_1240k_vcf/v44.3.1240K.vcf.gz
```

#### Split the genotype file in VCF by chromosome

```python
! for chr in {1..22}; \
    do bcftools view -r $chr panel_1240k_vcf/v44.3.1240K.vcf.gz \
                     -O z -o panel_1240k_vcf_chr/v44.3.1240K_chr${chr}.vcf.gz; \
  done; \
```

## Section II: Data collection

### Collect genotype data of modern AJ individuals

### We obtained, processed and merged six datasets of autosomal Jewish genotype data. The merged dataset consists of 1413 Jewish individuals before Quality Control.

#### 1. 471 AJs from Bray et al. 2010 (https://pubmed.ncbi.nlm.nih.gov/20798349/) and 7 AJs from Lazaridis et al. 2016 (https://pubmed.ncbi.nlm.nih.gov/27459054/) at ~144K SNPs
#### 2. 20 AJs and 101 other Jews from Behar et al. 2010 (https://pubmed.ncbi.nlm.nih.gov/20531471/) at ~570K SNPs
#### 3. 1 AJ and 83 other Jews from Behar et al. 2013 (https://pubmed.ncbi.nlm.nih.gov/25079123/) at ~300K SNPs
#### 4. 343 AJs (221 mixed-origin and 122 single-origin) from Gladstein et al. 2019 (https://pubmed.ncbi.nlm.nih.gov/30840069/) at ~710K SNPs
#### 5. 147 AJs and 240 other Jews from Kopelman et al. 2020 (https://pubmed.ncbi.nlm.nih.gov/31919450/) at ~470K SNPs

#### The annotation file for the pre-QC 1413 Jewish individuals is hosted at 

#### The original merged dataset is hosted at 

#### The merged dataset post extraction of Genochip-SNPs is hosted at 

#### The merged dataset post imputation (for IBD inference) is hosted at

### Collect genotype data of modern non-Jewish individuals
