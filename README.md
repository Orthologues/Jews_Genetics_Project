# Lund University BINP50 Master's Thesis Project in bioinformatics: 
# Identifying a fine genetic structure across Ashkenazic Jews
## April 2021 - August 2021
### Author: Jiawei Zhao (ji8842zh-s@student.lu.se)
### Supervisor: Prof. Eran Elhaik

<hr />
<br />

[Section I: Basic Setup for the project](#sec1) <br />
[Section II: Download and extract genotype data of Jews for IBD and GPS analyses](#sec2) <br />
[Section III: Extract genotype data of non-Jews from regions of interest in Afro-Eurasia for subsequent IBD analyses](#sec3) <br />
[Section IV: IBD inference by <code>phasedibd</code> and Extract of IBD segments](#sec4) <br />
+ [Run phasedibd](#run-phasedibd)
+ [Removing replicate individuals and close relatives(minimally half-cousins) from Jewish individuals](#jew-qc)

<hr />

<a name="sec1"></a>
# Section I: Basic Setup for the project

### Create a directory called "../main_folder" to conduct our major analysis. However, the directory should be git ignore or it would be too large to be put on Github.

```bash
%%bash
mkdir -p ../main_folder && cd ../main_folder
```


#### Download PLINK 1.9

```bash
%%bash
curl -sSL https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip -o plink.zip &&
unzip plink.zip
```

#### Download the tool to convert .geno & .snp & .ind formats into PLINK formats: https://github.com/DReichLab/EIG


```bash
%%bash
git clone https://github.com/DReichLab/EIG && ls -l EIG/CONVERTF/*
```

```bash
%%bash
cat EIG/CONVERTF/par.ANCESTRYMAP.EIGENSTRAT
```

    genotypename:    example.ancestrymapgeno
    snpname:         example.snp
    indivname:       example.ind
    outputformat:    EIGENSTRAT
    genotypeoutname: example.eigenstratgeno
    snpoutname:      example.snp
    indivoutname:    example.ind



```bash
%%bash
cat EIG/CONVERTF/par.PED.PACKEDPED
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


```bash
%%bash
cd EIG/src && make && make install && cd ../..
```

<a name="sec2"></a>
# Section II: Download and extract genotype data of Jews for IBD and GPS analyses

### Collect genotype data of modern AJ individuals

### We obtained, processed and merged four datasets of autosomal Jewish genotype data. The merged dataset consists of 935 Jewish individuals before Quality Control.

#### 1. 20 AJs and 101 other Jews from Behar et al. 2010 (https://pubmed.ncbi.nlm.nih.gov/20531471/) at ~570K SNPs
#### 2. 1 AJ and 83 other Jews from Behar et al. 2013 (https://pubmed.ncbi.nlm.nih.gov/25079123/) at ~300K SNPs
#### 3. 343 AJs (221 mixed-origin and 122 single-origin) from Gladstein et al. 2019 (https://pubmed.ncbi.nlm.nih.gov/30840069/) at ~710K SNPs
#### 4. 147 AJs and 240 other Jews from Kopelman et al. 2020 (https://pubmed.ncbi.nlm.nih.gov/31919450/) at ~470K SNPs

#### The annotation file for the 935 pre-QC Jewish individuals is hosted at <a>https://github.com/Orthologues/Jews-Genetics-Project/blob/main/pre_QC_935_Jews_annotation.csv</a>

#### The original merged dataset is hosted at <a>https://github.com/Orthologues/Jews-Genetics-Project/blob/main/pre_QC_935_Jews_original_gt.tar.gz</a>

#### The text file for autosomal genochip SNPs is hosted at <a>https://github.com/Orthologues/Jews_Genetics_Project/blob/main/genochip_autosomal_snps.txt</a>

#### The merged dataset post imputation (for IBD inference) is hosted at <a>https://github.com/Orthologues/Jews_Genetics_Project/blob/main/HRC_eu_ref_imputed_935Jews_for_IBD.tar.gz</a>

### Clone the github repository of this project and unzip its genotype files

```bash
%%bash
git clone git@github.com:Orthologues/Jews_Genetics_Project &&
tar -xzvf Jews_Genetics_Project/pre_QC_935Jews_original_gt.tar.gz &&
tar -xzvf Jews_Genetics_Project/HRC_eu_ref_imputed_935Jews_for_IBD.tar.gz
```

### Extract autosomal genochip SNPs from the original genotype dataset of Jews for ADMIXTURE and GPS analyses
```bash
%%bash
plink --bfile Jews_Genetics_Project/pre_QC_935Jews_original_gt/pre_QC_935Jews_original_gt --extract Jews_Genetics_Project/genochip_autosomal_snps.txt \
--allow-no-sex --make-bed --out pre_QC_935Jews_genochip
```

<a name="sec3"></a>
# Section III: Extract genotype data of non-Jews from regions of interest in Afro-Eurasia for subsequent IBD analyses

### In order to perform a more comprehensive IBD analysis, we would need to extract modern populations from regions of interest available at https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data from the v44.3 1240K+HO panel

### Link to the original .tar file: <a>https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V44/V44.3/SHARE/public.dir/v44.3_HO_public.tar</a>

### Reference for searching samples: https://www.coriell.org/0/Sections/Search/

### An extra high-coverage SNP dataset consisting of non-Jewish modern individuals: 
#### Paper: https://www.nature.com/articles/nature19792 (Pagani et. al., 2016)
#### Link to data: https://evolbio.ut.ee/


```bash
%%bash
cd ../main_folder/
```


```bash
%%bash
mkdir -p 1240K_HO_pop_of_interest
```


```bash
%%bash
cd 1240K_HO_pop_of_interest/
```

```bash
%%bash
curl -sSL \
https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V44/V44.3/SHARE/public.dir/\
v44.3_HO_public.tar -o v44.3_HO_public.tar
```


```bash
%%bash
tar -xf v44.3_HO_public.tar && rm v44.3_HO_public.tar
```


```bash
%%bash
head -n 5 v44.3_HO_public.ind
```

                 MAL-005 M Malawi_Yao
                 MAL-009 M Malawi_Yao
                 MAL-011 M Malawi_Chewa
                 MAL-012 M Malawi_Chewa
                 MAL-014 M Malawi_Chewa



```bash
%%bash
echo "\
genotypename:    v44.3_HO_public.geno
snpname:         v44.3_HO_public.snp
indivname:       v44.3_HO_public.ind
outputformat:    PACKEDPED
genotypeoutname: v44.3_HO_public.bed
snpoutname:      v44.3_HO_public.bim
indivoutname:    v44.3_HO_public.fam
familynames:     NO" > par.HO_MAP.PACKEDPED
```


```bash
%%bash
../EIG/src/convertf -p par.HO_MAP.PACKEDPED
```

    parameter file: par.HO_MAP.PACKEDPED
    genotypename: v44.3_HO_public.geno
    snpname: v44.3_HO_public.snp
    indivname: v44.3_HO_public.ind
    outputformat: PACKEDPED
    genotypeoutname: v44.3_HO_public.bed
    snpoutname: v44.3_HO_public.bim
    indivoutname: v44.3_HO_public.fam
    familynames: NO
    read 1073741824 bytes
    read 1971990900 bytes
    packed geno read OK
    end of inpack
    numvalidind:  13197  maxmiss: 13197001
    packedped output
    ##end of convertf run



```bash
%%bash
head -n 5 v44.3_HO_public.fam
```

         1        MAL-005 0 0 1 1
         2        MAL-009 0 0 1 1
         3        MAL-011 0 0 1 1
         4        MAL-012 0 0 1 1
         5        MAL-014 0 0 1 1



```bash
%%bash
grep Armenian v44.3_HO_public.ind
```

                  ARK-59 M Armenian_Hemsheni
                  ARK-72 M Armenian_Hemsheni
                  ARK-76 M Armenian_Hemsheni
                  ARK-78 M Armenian_Hemsheni
                  ARK-82 M Armenian_Hemsheni
                  ARK-84 M Armenian_Hemsheni
                  ARK-88 M Armenian_Hemsheni
                  ARK-89 M Armenian_Hemsheni
                  ARM012 M Armenian.WGA
                  ARM013 M Armenian.WGA
                  ARM014 M Armenian.WGA
              armenia176 M   Armenian
              armenia191 M   Armenian
               armenia86 M   Armenian
              armenia279 M   Armenian
               armenia91 M   Armenian
              armenia293 M   Armenian
              armenia102 M   Armenian
              armenia106 M   Armenian
              armenia139 M   Armenian
              armenia162 M   Armenian
         S_Armenian-1.DG M Armenian.DG
         S_Armenian-2.DG M Armenian.DG


### Suitable populations in summary: 
#### Turkic : Turkish, Chuvash, Tatar, Bashkir, Uzbek, Turkmen
#### Mid-eastern and Iranic: Iranian, Iranian(Bandari, Gulf region), Bedouin, Druze, Palestinian, Syrian, Assyrian, Kurd, Lebanese (Christian, Muslim, unknown), Jordanian, Saudi, Pathan
#### Caucasus: Adygei, Armenian, Hemsheni, Chechen, Georgian, Azeri,  Ingushian, Avar, Darginian
#### Eastern-Europe: Albanian, Bulgarian, Croatian, Greek, Hungarian, Polish(Only 1+plus 4 from Pagani et. al., 2016), Czech, Ukrainian, Belarusian, Moldavian, Mordovian(Ugric people in Russia) Romanian, Serbian(0), Slovakian(0), Russian, Latvian(3 from Pagani et. al., 2016), Lithuanian, Estonian, German(3 from Pagani et. al., 2016), Swedes(2 from Pagani et. al., 2016), Romani(3 from Pagani et. al., 2016)
#### Western-Europe: Dutch(0), German(0), Basque, French, French_south, Northern Italian, Southern Italian, Sardinian, Tuscan, Spanish, Northern Spanish, Portuguese(0), English, Scottish, Icelandic, Maltese
#### Northern Europe: Norwegian, Finnish
#### North Africa: Tunisian, Algerian, Egyptian, Moroccan, Libyan


## Mid-eastern and Iranic

```bash
%%bash
grep Kurd v44.3_HO_public.ind|grep -v [Jjew]
```

                   KRD-5 F       Kurd
                   KRD-8 F       Kurd
                  KRD-14 F       Kurd
                  KRD-31 M       Kurd
                  KRD-40 M       Kurd
                  KRD-43 F       Kurd
                  KRD-45 F       Kurd
                  KRD-49 M       Kurd
                 KRD_010 M   Kurd.WGA
                 KRD_011 M   Kurd.WGA



```bash
%%bash
grep Kurd v44.3_HO_public.ind|grep -v [Jjew]|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Kurd.id && wc -l Kurd.id
```

    10 Kurd.id



```bash
%%bash
grep Saudi v44.3_HO_public.ind|grep -v [Jjew]
```

                 SaudiA5 M      Saudi
                 SaudiA6 M      Saudi
                 SaudiA7 M      Saudi
               saudi1434 M      Saudi
               saudi1424 M      Saudi
                 SaudiA9 M      Saudi
                 SaudiA1 M      Saudi
               saudi1403 M      Saudi



```bash
%%bash
grep Saudi v44.3_HO_public.ind|grep -v [Jjew]|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Saudi.id && wc -l Saudi.id
```

    8 Saudi.id



```bash
%%bash
grep Jordanian v44.3_HO_public.ind|grep -v [Jj]ew
```

               Jordan543 M  Jordanian
               Jordan214 M  Jordanian
               Jordan445 M  Jordanian
                Jordan62 M  Jordanian
               Jordan603 M  Jordanian
               Jordan307 M  Jordanian
               Jordan444 M Ignore_Jordanian
               Jordan503 M  Jordanian
               Jordan646 M  Jordanian
               Jordan384 M  Jordanian
        S_Jordanian-2.DG M Jordanian.DG
        S_Jordanian-3.DG M Jordanian.DG
        S_Jordanian-1.DG M Jordanian.DG



```bash
%%bash
grep Jordanian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Jordanian.id && wc -l Jordanian.id
```

    9 Jordanian.id



```bash
%%bash
grep Lebanese v44.3_HO_public.ind|grep -v [Jj]ew
```

                Lebanon5 M   Lebanese
                Lebanon6 M   Lebanese
                Lebanon7 M   Lebanese
                Lebanon8 M   Lebanese
                Lebanon1 M   Lebanese
                Lebanon2 M   Lebanese
                Lebanon3 M   Lebanese
                Lebanon4 M   Lebanese
          Lebanese1AQ127 M Lebanese_Christian
          Lebanese1AQ170 F Lebanese_Christian
          Lebanese2AQ121 F Lebanese_Muslim
          Lebanese2AQ127 M Lebanese_Muslim
          Lebanese4AQ115 M Lebanese_Christian
          Lebanese6AQ115 M Lebanese_Christian
          Lebanese6AQ170 M Lebanese_Christian
           Lebanese6AS15 M Lebanese_Muslim
          Lebanese7AQ150 M Lebanese_Muslim
           Lebanese7AR20 M Lebanese_Muslim
           Lebanese7AR23 M Lebanese_Muslim
           Lebanese8AS15 M Lebanese_Christian
         Lebanese10AQ127 M Lebanese_Muslim
          Lebanese10AR37 M Lebanese_Christian
          Lebanese11AS14 F Lebanese_Muslim
          Lebanese15AR37 M Lebanese_Christian
          Lebanese20AR21 F Lebanese_Muslim
          Lebanese22BA23 F Lebanese_Christian
          Lebanese24AR27 M Lebanese_Muslim
          Lebanese30AR21 M Lebanese_Muslim



```bash
%%bash
grep Lebanese_Christian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Lebanese_Christian.id && wc -l Lebanese_Christian.id
```

    9 Lebanese_Christian.id



```bash
%%bash
grep Lebanese_Muslim v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Lebanese_Muslim.id && wc -l Lebanese_Muslim.id
```

    11 Lebanese_Muslim.id



```bash
%%bash
grep Lebanese v44.3_HO_public.ind|grep -v Muslim|grep -v Christian|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Lebanese_unknown.id && wc -l Lebanese_unknown.id
```

    8 Lebanese_unknown.id



```bash
%%bash
grep Iranian v44.3_HO_public.ind|grep -v [Jj]ew
```

                  iran19 M    Iranian
                  iran14 M    Iranian
                   iran2 M    Iranian
                  iran20 M    Iranian
                   iran3 M    Iranian
                   iran7 M Ignore_Iranian
                  iran11 M    Iranian
                  iran16 M    Iranian
                  iran17 M    Iranian
                   PV001 M Iranian_Bandari
                   PV002 M Iranian_Bandari
                   PV003 M Iranian_Bandari
                   PV004 M Iranian_Bandari
                   PV005 M Iranian_Bandari
                   PV006 M Iranian_Bandari
                   PV007 M Iranian_Bandari
                   PV008 M Iranian_Bandari
                   PV019 M    Iranian
                   PV020 M    Iranian
                   PV021 M    Iranian
                   PV022 M    Iranian
                   PV023 M    Iranian
                   PV024 M    Iranian
                   PV025 M    Iranian
                   PV026 M    Iranian
                   PV027 M    Iranian
                   PV028 M    Iranian
                   PV009 M    Iranian
                   PV010 M    Iranian
                   PV011 M    Iranian
                   PV012 M    Iranian
                   PV013 M    Iranian
                   PV014 M    Iranian
                   PV015 M    Iranian
                   PV016 M    Iranian
                   PV017 M    Iranian
                   PV018 M    Iranian
                   PV029 M    Iranian
                   PV030 M    Iranian
                   PV031 M    Iranian
                   PV032 M    Iranian
                   PV033 M    Iranian
                   PV034 M    Iranian
                   PV035 M    Iranian
                   PV036 M    Iranian
                   PV037 M    Iranian
                   PV038 M    Iranian
          S_Iranian-1.DG M Iranian.DG
          S_Iranian-2.DG M Iranian.DG


#### "Bandari" shall be regarded as Southern Iranians. It status would be considered in our analysis later


```bash
%%bash
grep Iranian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|grep -v Bandari|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Iranian.id && wc -l Iranian.id
```

    38 Iranian.id



```bash
%%bash
grep Iranian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|grep Bandari|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Iranian_gulf.id && wc -l Iranian_gulf.id
```

    8 Iranian_gulf.id



```bash
%%bash
grep Palestinian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Palestinian.id && wc -l Palestinian.id
```

    38 Palestinian.id



```bash
%%bash
grep Bedouin v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Bedouin.id && wc -l Bedouin.id
```

    44 Bedouin.id



```bash
%%bash
grep Druze v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|grep -v Africa|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Druze.id && wc -l Druze.id
```

    39 Druze.id



```bash
%%bash
grep Syrian v44.3_HO_public.ind
```

                syria464 M     Syrian
                syria361 M     Syrian
                syria485 M     Syrian
                syria520 M     Syrian
                  syria4 M     Syrian
                syria461 M     Syrian
                  syria6 M     Syrian
                  syria7 M     Syrian



```bash
%%bash
grep Syrian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Syrian.id && wc -l Syrian.id
```

    8 Syrian.id



```bash
%%bash
grep Assyrian v44.3_HO_public.ind
```

                 ASR_001 M Assyrian.WGA
                 ASR_002 M Assyrian.WGA
                 ASR_003 M Assyrian.WGA
                 ASR_004 M Assyrian.WGA_outlier
                 ASR_005 M Assyrian.WGA
             Assyrian151 F   Assyrian
             Assyrian152 M   Assyrian
             Assyrian153 M   Assyrian
             Assyrian155 M   Assyrian
             Assyrian159 F   Assyrian
             Assyrian160 M   Assyrian
             Assyrian161 M   Assyrian
             Assyrian162 M   Assyrian
             Assyrian163 F   Assyrian
             Assyrian164 M   Assyrian
             Assyrian165 F   Assyrian



```bash
%%bash
grep Assyrian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Assyrian.id && wc -l Assyrian.id
```

    15 Assyrian.id



```bash
%%bash
grep Pathan v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Pathan.id && wc -l Pathan.id
```

    17 Pathan.id


## North Africa


```bash
%%bash
grep Egyptian v44.3_HO_public.ind
```

                  Egypt7 M   Egyptian
                  Egypt8 M Ignore_Egyptian
                  Egypt9 M   Egyptian
                 Egypt10 M   Egyptian
                  Egypt3 M   Egyptian
                 Egypt11 M   Egyptian
                  Egypt1 M   Egyptian
                 Egypt12 M   Egyptian
              Egypt8BA65 M Ignore_Egyptian
            Egypt17AQ176 M   Egyptian
            Egypt19AQ172 M Ignore_Egyptian
            Egypt20AQ172 M Ignore_Egyptian
            Egypt15AQ172 M   Egyptian
             Egypt3AQ172 M   Egyptian
             Egypt1AQ172 M   Egyptian
             Egypt5AQ172 M   Egyptian
             Egypt8AT113 M   Egyptian
            Egypt12AQ172 M   Egyptian
             Egypt9AQ177 M   Egyptian
             Egypt9AQ172 M   Egyptian
             Egypt7AQ172 M   Egyptian
             Egypt22TD21 M   Egyptian



```bash
%%bash
grep Egyptian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Egyptian.id && wc -l Egyptian.id
```

    18 Egyptian.id



```bash
%%bash
grep Libyan v44.3_HO_public.ind
```

           LibyanJew1438 F Jew_Libyan
           LibyanJew1601 M Jew_Libyan
           LibyanJew1104 F Jew_Libyan
           LibyanJew1462 F Jew_Libyan
           LibyanJew1263 F Jew_Libyan
           LibyanJew1611 M Jew_Libyan
           LibyanJew1605 M Jew_Libyan
           LibyanJew1639 M Jew_Libyan
           LibyanJew1659 M Jew_Libyan
                    LIB7 M     Libyan
                   LIB13 M     Libyan
                   LIB18 M     Libyan
                   LIB27 M     Libyan
                   LIB30 F     Libyan



```bash
%%bash
grep Libyan v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Libyan.id && wc -l Libyan.id
```

    5 Libyan.id



```bash
%%bash
grep Tunisian v44.3_HO_public.ind
```

          Tunisian200000 F   Tunisian
         TunisianJew1421 F Jew_Tunisian
            Tunisian20C1 M   Tunisian
            Tunisian20F4 M   Tunisian
         TunisianJew1763 M Jew_Tunisian
         TunisianJew1507 M Jew_Tunisian
            Tunisian20D4 F   Tunisian
         TunisianJew1544 F Jew_Tunisian
            Tunisian20A5 F   Tunisian
            Tunisian20C4 F   Tunisian
         TunisianJew1170 M Jew_Tunisian
         TunisianJew1531 F Jew_Tunisian
         TunisianJew1511 F Jew_Tunisian
            Tunisian20D1 F   Tunisian
            Tunisian20B4 M   Tunisian



```bash
%%bash
grep Tunisian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Tunisian.id && wc -l Tunisian.id
```

    8 Tunisian.id



```bash
%%bash
grep Algerian v44.3_HO_public.ind
```

           Algerian43A22 F   Algerian
           Algerian43A21 F   Algerian
           Algerian43A34 M   Algerian
           Algerian43A13 M   Algerian
           Algerian43A24 F   Algerian
           Algerian43A32 F   Algerian
           Algerian43A23 F   Algerian



```bash
%%bash
grep Algerian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Algerian.id && wc -l Algerian.id
```

    7 Algerian.id



```bash
%%bash
grep Moroccan v44.3_HO_public.ind
```

         MoroccanJew5134 M Ignore_Jew_Moroccan
         MoroccanJew5126 M Jew_Moroccan
         MoroccanJew4634 M Jew_Moroccan
         MoroccanJew4789 F Jew_Moroccan
         MoroccanJew5168 M Jew_Moroccan
         MoroccanJew4692 F Jew_Moroccan
         MoroccanJew4683 F Jew_Moroccan
                    MBE3 M Ignore_Moroccan
                   MBE11 M Ignore_Moroccan
                   MBE13 M Ignore_Moroccan
                   MBE16 M Ignore_Moroccan
                   MBE19 M Ignore_Moroccan
                   MBE22 M Ignore_Moroccan
                   MBE23 M Ignore_Moroccan
                    MCA7 M   Moroccan
                    MCA8 M   Moroccan
                    MCA9 F   Moroccan
                   MCA14 M   Moroccan
                   MCA16 M   Moroccan
                   MCA19 M   Moroccan
                   MCA24 F   Moroccan
                   MCA37 M   Moroccan
                   MCA38 M   Moroccan
                   MCA39 F   Moroccan



```bash
%%bash
grep Moroccan v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Moroccan.id && wc -l Moroccan.id
```

    10 Moroccan.id


## Northern Europe


```bash
%%bash
grep Finnish v44.3_HO_public.ind
```

                 HG00171 F    Finnish
                 HG00174 F    Finnish
                 HG00190 M    Finnish
                 HG00266 F    Finnish
                 HG00183 M    Finnish
                 HG00173 F    Finnish
                 HG00181 M    Finnish
                 HG00182 M    Finnish
          S_Finnish-3.DG M Finnish.DG
          S_Finnish-1.DG F Finnish.DG
          S_Finnish-2.DG M Finnish.DG



```bash
%%bash
grep Finnish v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Finnish.id && wc -l Finnish.id
```

    8 Finnish.id



```bash
%%bash
grep Norwegian v44.3_HO_public.ind
```

                  NOR119 M  Norwegian
                  NOR124 M  Norwegian
                  NOR106 M  Norwegian
                  NOR101 M  Norwegian
                  NOR146 M  Norwegian
                  NOR108 M  Norwegian
                  NOR126 M  Norwegian
                  NOR107 M  Norwegian
                  NOR109 M  Norwegian
                  NOR148 M  Norwegian
                  NOR111 M  Norwegian
        S_Norwegian-1.DG U Norwegian.DG



```bash
%%bash
grep Norwegian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Norwegian.id && wc -l Norwegian.id
```

    11 Norwegian.id


## Western European


```bash
%%bash
grep Maltese v44.3_HO_public.ind
```

              Malta4AM91 M    Maltese
              Malta8AM91 M    Maltese
              Malta7AM91 M    Maltese
             Malta17AM91 M    Maltese
              Malta2AM91 M    Maltese
             Malta16AM91 M    Maltese
             Malta15AM91 M    Maltese
             Malta12AM91 M    Maltese



```bash
%%bash
grep Maltese v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Maltese.id && wc -l Maltese.id
```

    8 Maltese.id



```bash
%%bash
grep Spanish v44.3_HO_public.ind|head -n 20 && grep Spanish v44.3_HO_public.ind|wc -l
```

                 HG01500 M    Spanish
                 HG01501 F    Spanish
                 HG01503 M    Spanish
                 HG01504 F    Spanish
                 HG01506 M    Spanish
                 HG01507 F    Spanish
                 HG01509 M    Spanish
                 HG01510 F    Spanish
                 HG01512 M    Spanish
                 HG01513 F    Spanish
                 HG01515 M Spanish_North
                 HG01516 F Spanish_North
                 HG01518 M Spanish_North
                 HG01524 M    Spanish
                 HG01527 M    Spanish
                 HG01528 F    Spanish
                 HG01530 M    Spanish
                 HG01536 M    Spanish
                 HG01537 F    Spanish
                 HG01605 F    Spanish
    181



```bash
%%bash
grep Spanish v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|grep -v North|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Spanish.id && wc -l Spanish.id
```

    173 Spanish.id



```bash
%%bash
grep Spanish v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|grep North|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Spanish_north.id && wc -l Spanish_north.id
```

    5 Spanish_north.id


### Therefore, we would like to seperate "Spanish" and "Spanish_North" (the latter closer to Basque?)


```bash
%%bash
grep Portuguese v44.3_HO_public.ind ## No Portuguese individuals
```


```bash
%%bash
grep Basque v44.3_HO_public.ind
```

               HGDP01357 M     Basque
               HGDP01358 M     Basque
               HGDP01359 M     Basque
               HGDP01360 M     Basque
               HGDP01362 M     Basque
               HGDP01363 F     Basque
               HGDP01364 M     Basque
               HGDP01365 F     Basque
               HGDP01366 F     Basque
               HGDP01367 F     Basque
               HGDP01368 U     Basque
               HGDP01369 F Ignore_Basque
               HGDP01370 M     Basque
               HGDP01371 M     Basque
               HGDP01372 M Ignore_Basque
               HGDP01373 F     Basque
               HGDP01374 M     Basque
               HGDP01375 M     Basque
               HGDP01377 M     Basque
               HGDP01378 M     Basque
               HGDP01379 M     Basque
               HGDP01380 F     Basque
                   BAS35 M     Basque
                   BAS31 F     Basque
                   BAS22 F     Basque
                   BAS25 M     Basque
                   BAS32 M     Basque
                   BAS30 F     Basque
                   BAS29 F Ignore_Basque
                   BAS28 F     Basque
                   BAS27 M     Basque
                   BAS33 F     Basque
           S_Basque-2.DG F  Basque.DG
           S_Basque-1.DG M  Basque.DG
           HGDP01368.SDG U Basque.SDG
           HGDP01361.SDG M Basque.SDG
           HGDP01358.SDG M Basque.SDG
           HGDP01359.SDG M Basque.SDG
           HGDP01357.SDG M Basque.SDG
           HGDP01363.SDG F Basque.SDG
           HGDP01360.SDG M Ignore_Basque.SDG
           HGDP01362.SDG M Basque.SDG
           HGDP01376.SDG M Basque.SDG
           HGDP01375.SDG M Basque.SDG
           HGDP01367.SDG F Basque.SDG
           HGDP01366.SDG F Basque.SDG
           HGDP01369.SDG F Ignore_Basque.SDG
           HGDP01379.SDG M Basque.SDG
           HGDP01378.SDG M Basque.SDG
           HGDP01373.SDG F Basque.SDG
           HGDP01370.SDG M Basque.SDG
           HGDP01372.SDG M Ignore_Basque.SDG
           HGDP01377.SDG M Basque.SDG
           HGDP01380.SDG F Basque.SDG
           HGDP01374.SDG M Basque.SDG
           HGDP01364.SDG M Basque.SDG
           HGDP01365.SDG F Basque.SDG



```bash
%%bash
grep Basque v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Basque.id && wc -l Basque.id
```

    29 Basque.id



```bash
%%bash
grep French v44.3_HO_public.ind
```

               HGDP00511 M     French
               HGDP00512 M     French
               HGDP00513 F     French
               HGDP00514 F     French
               HGDP00515 M     French
               HGDP00516 F     French
               HGDP00517 F     French
               HGDP00518 M     French
               HGDP00519 M     French
               HGDP00520 F Ignore_French
               HGDP00521 M Ignore_French(discovery)
               HGDP00522 M     French
               HGDP00523 F     French
               HGDP00524 F     French
               HGDP00525 M     French
               HGDP00526 F     French
               HGDP00527 F     French
               HGDP00528 M     French
               HGDP00529 F     French
               HGDP00530 M Ignore_French
               HGDP00531 F     French
               HGDP00533 M     French
               HGDP00534 F     French
               HGDP00535 F     French
               HGDP00536 F     French
               HGDP00537 F     French
               HGDP00538 M     French
               HGDP00539 F     French
         SouthFrench3326 M     French
         SouthFrench3947 M     French
         SouthFrench1323 M     French
         SouthFrench3951 M     French
         SouthFrench3068 M     French
         SouthFrench1112 M     French
         SouthFrench4018 M     French
             French23812 M     French
             French23814 M     French
             French23821 F     French
             French23830 M     French
             French23833 F     French
             French23862 M     French
             French23915 M     French
             French23919 F     French
             French23952 M Ignore_French
             French23989 F     French
             French24061 M     French
             French24075 F     French
             French24076 M     French
             French24090 M     French
             French24118 M     French
             French24120 F     French
             French24124 M     French
             French24144 F     French
             French24148 M     French
             French24178 F     French
             French24247 F     French
             French24381 F     French
             French24400 M     French
             French24408 F     French
             French24433 M     French
             French24434 F     French
             French24437 F     French
             French24690 U     French
             French24817 M     French
             French25068 F     French
           HGDP00521_WGA M Ignore_French(discovery)
           A_French-4.DG M Ignore_French(discovery).DG
           S_French-1.DG M  French.DG
           B_French-3.DG M  French.DG
           S_French-2.DG F  French.DG
                  TAP002 M French_Polynesia_200BP
                  TAP003 M French_Polynesia_400BP
                  TAP004 M French_Polynesia_200BP
           HGDP00536.SDG F French.SDG
           HGDP00539.SDG F French.SDG
           HGDP00519.SDG M French.SDG
           HGDP00512.SDG M French.SDG
           HGDP00513.SDG F French.SDG
           HGDP00527.SDG F French.SDG
           HGDP00524.SDG F French.SDG
           HGDP00511.SDG M French.SDG
           HGDP00520.SDG F Ignore_French.SDG
           HGDP00522.SDG M French.SDG
           HGDP00535.SDG F French.SDG
           HGDP00523.SDG F French.SDG
           HGDP00518.SDG M French.SDG
           HGDP00531.SDG F French.SDG
           HGDP00537.SDG F French.SDG
           HGDP00528.SDG M French.SDG
           HGDP00516.SDG F French.SDG
           HGDP00529.SDG F French.SDG
           HGDP00515.SDG M French.SDG
           HGDP00538.SDG M French.SDG
           HGDP00525.SDG M French.SDG
           HGDP00517.SDG F French.SDG
           HGDP00514.SDG F French.SDG
           HGDP00534.SDG F French.SDG
           HGDP00521.SDG M Ignore_French(discovery).SDG
           HGDP00526.SDG F French.SDG
           HGDP00533.SDG M French.SDG
           HGDP00530.SDG M Ignore_French.SDG



```bash
%%bash
grep French v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|grep -v Polynesia|\
sed -r 's/^\s+//'|grep -v [Ss]outh|cut -d ' ' -f 1 > French.id && wc -l French.id
```

    54 French.id



```bash
%%bash
grep French v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|grep -v Polynesia|\
sed -r 's/^\s+//'|grep [Ss]outh|cut -d ' ' -f 1 > French_south.id && wc -l French_south.id
```

    7 French_south.id



```bash
%%bash
grep Sardinian v44.3_HO_public.ind
```

               HGDP00665 M Ignore_Sardinian(discovery)
               HGDP00666 M  Sardinian
               HGDP00667 F  Sardinian
               HGDP00668 M  Sardinian
               HGDP00669 F  Sardinian
               HGDP00670 M  Sardinian
               HGDP00671 M  Sardinian
               HGDP00672 F  Sardinian
               HGDP00673 F  Sardinian
               HGDP00674 M  Sardinian
               HGDP01062 F  Sardinian
               HGDP01063 M  Sardinian
               HGDP01064 F  Sardinian
               HGDP01065 F  Sardinian
               HGDP01066 M  Sardinian
               HGDP01067 M  Sardinian
               HGDP01068 F  Sardinian
               HGDP01069 M  Sardinian
               HGDP01070 F  Sardinian
               HGDP01071 M  Sardinian
               HGDP01072 F  Sardinian
               HGDP01073 M  Sardinian
               HGDP01074 F  Sardinian
               HGDP01075 M  Sardinian
               HGDP01076 M  Sardinian
               HGDP01077 M  Sardinian
               HGDP01078 F  Sardinian
               HGDP01079 M  Sardinian
           HGDP00665_WGA M Ignore_Italian_Sardinian(discovery)
        A_Sardinian-4.DG M Ignore_Sardinian(discovery).DG
        B_Sardinian-3.DG M Sardinian.DG
        S_Sardinian-1.DG M Sardinian.DG
        S_Sardinian-2.DG F Sardinian.DG
           HGDP01074.SDG F Sardinian.SDG
           HGDP00671.SDG M Sardinian.SDG
           HGDP01073.SDG M Sardinian.SDG
           HGDP01065.SDG F Sardinian.SDG
           HGDP01066.SDG M Sardinian.SDG
           HGDP00670.SDG M Sardinian.SDG
           HGDP01071.SDG M Sardinian.SDG
           HGDP00668.SDG M Sardinian.SDG
           HGDP00674.SDG M Sardinian.SDG
           HGDP01068.SDG F Sardinian.SDG
           HGDP01075.SDG M Sardinian.SDG
           HGDP01070.SDG F Sardinian.SDG
           HGDP00667.SDG F Sardinian.SDG
           HGDP01063.SDG M Sardinian.SDG
           HGDP01067.SDG M Sardinian.SDG
           HGDP01069.SDG M Sardinian.SDG
           HGDP00673.SDG F Sardinian.SDG
           HGDP00669.SDG F Sardinian.SDG
           HGDP01072.SDG F Sardinian.SDG
           HGDP01077.SDG M Sardinian.SDG
           HGDP00666.SDG M Sardinian.SDG
           HGDP01064.SDG F Sardinian.SDG
           HGDP00665.SDG M Ignore_Sardinian(discovery).SDG
           HGDP01079.SDG M Sardinian.SDG
           HGDP01078.SDG F Sardinian.SDG
           HGDP01076.SDG M Sardinian.SDG



```bash
%%bash
grep Sardinian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Sardinian.id && wc -l Sardinian.id
```

    27 Sardinian.id



```bash
%%bash
grep Tuscan v44.3_HO_public.ind
```

           S_Tuscan-1.DG F Tuscan_1.DG
           S_Tuscan-2.DG M Tuscan_1.DG
           HGDP00672.SDG F Tuscan_2.SDG
           HGDP01062.SDG F Tuscan_2.SDG



```bash
%%bash
grep Tuscan v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Tuscan.id && wc -l Tuscan.id
```

    4 Tuscan.id



```bash
%%bash
grep Italian_South v44.3_HO_public.ind
```

                   BEL57 M Italian_South
                    ITS2 M Ignore_Italian_South(first_degree_relative)
                    ITS4 F Italian_South
                    ITS5 F Italian_South
                    ITS7 M Italian_South



```bash
%%bash
grep Italian_South v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Italian_south.id && wc -l Italian_south.id
```

    4 Italian_south.id



```bash
%%bash
grep Italian_North v44.3_HO_public.ind
```

               HGDP01147 M Italian_North
               HGDP01149 M Ignore_Italian_North
               HGDP01151 M Italian_North
               HGDP01152 M Italian_North
               HGDP01153 M Italian_North
               HGDP01155 M Italian_North
               HGDP01156 F Italian_North
               HGDP01157 U Italian_North
               HGDP01161 M Italian_North
               HGDP01162 M Italian_North
               HGDP01163 M Italian_North
               HGDP01164 M Italian_North
               HGDP01166 M Italian_North
               HGDP01167 M Italian_North
               HGDP01168 F Italian_North
               HGDP01169 F Italian_North
               HGDP01171 F Italian_North
               HGDP01172 F Italian_North
               HGDP01173 M Italian_North
               HGDP01174 M Italian_North
               HGDP01177 F Italian_North
          S_Bergamo-1.DG M Italian_North.DG
           HGDP01157.SDG U Italian_North.SDG
           HGDP01164.SDG M Italian_North.SDG
           HGDP01156.SDG F Italian_North.SDG
           HGDP01166.SDG M Italian_North.SDG
           HGDP01173.SDG M Italian_North.SDG
           HGDP01171.SDG F Italian_North.SDG
           HGDP01169.SDG F Italian_North.SDG
           HGDP01149.SDG M Ignore_Italian_North.SDG
           HGDP01155.SDG M Italian_North.SDG
           HGDP01161.SDG M Italian_North.SDG
           HGDP01167.SDG M Italian_North.SDG
           HGDP01152.SDG M Italian_North.SDG
           HGDP01162.SDG M Italian_North.SDG
           HGDP01174.SDG M Italian_North.SDG
           HGDP01151.SDG M Italian_North.SDG
           HGDP01177.SDG F Italian_North.SDG
           HGDP01153.SDG M Italian_North.SDG
           HGDP01172.SDG F Italian_North.SDG
           HGDP01163.SDG M Italian_North.SDG
           HGDP01168.SDG F Italian_North.SDG



```bash
%%bash
grep Italian_North v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Italian_north.id && wc -l Italian_north.id
```

    20 Italian_north.id



```bash
%%bash
grep English v44.3_HO_public.ind
```

                 HG00128 F    English
                 HG00129 M    English
                 HG00130 F    English
                 HG00131 M    English
                 HG00233 F    English
                 HG00234 M    English
                 HG00232 F    English
                 HG00231 F    English
                 HG00126 M    English
                 HG00160 M    English
          S_English-1.DG M English.DG
          S_English-2.DG F English.DG



```bash
%%bash
grep English v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > English.id && wc -l English.id
```

    10 English.id



```bash
%%bash
grep Scottish v44.3_HO_public.ind
```

                 HG00103 M   Scottish
                 HG00104 F   Scottish
                 HG00105 M   Scottish
                 HG00106 F   Scottish



```bash
%%bash
grep Scottish v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Scottish.id && wc -l Scottish.id
```

    4 Scottish.id



```bash
%%bash
grep Icelandic v44.3_HO_public.ind
```

                 NA15762 F  Icelandic
                 NA15755 M  Icelandic
                 NA15763 F  Icelandic
                 NA15756 M  Icelandic
                 NA15764 F  Icelandic
                 NA15757 M  Icelandic
                 NA15765 F  Icelandic
                 NA15758 M  Icelandic
                 NA15766 F  Icelandic
                 NA15759 M  Icelandic
                 NA15760 M  Icelandic
                 NA15761 F  Icelandic
        S_Icelandic-2.DG F Icelandic.DG
        S_Icelandic-1.DG F Icelandic.DG



```bash
%%bash
grep Icelandic v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Icelandic.id && wc -l Icelandic.id
```

    12 Icelandic.id


## Eastern European people within EU


```bash
%%bash
grep Polish v44.3_HO_public.ind
```

           S_Polish-1.DG M  Polish.DG


#### Only one Polish present


```bash
%%bash
grep Polish v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Polish.id && wc -l Polish.id
```

    1 Polish.id



```bash
%%bash
grep Czech v44.3_HO_public.ind
```

                 NA15725 M      Czech
                 NA15733 F      Czech
                 NA15726 M      Czech
                 NA15727 M      Czech
                 NA15728 M      Czech
                 NA15729 F      Czech
                 NA15730 F      Czech
                 NA15731 F      Czech
                 NA15724 M      Czech
                 NA15732 F      Czech
                   I7207 M Czech_CordedWare
                   I7208 M Czech_CordedWare
                   I7209 M Czech_CordedWare
                   I7271 M Czech_BellBeaker_brother.I7278
                   I7278 M Czech_BellBeaker
                   I7282 M Czech_BellBeaker
                   I7283 F Czech_BellBeaker_mother.I7282
                   I7289 M Czech_BellBeaker_sibling.I7214
                   I7210 M Czech_BellBeaker_father.or.son.I7212
                   I7212 M Czech_BellBeaker
                   I7214 F Czech_BellBeaker
                   I6695 F Czech_CordedWare
                   I6696 M Czech_CordedWare
                   I7195 F  Czech_EBA
                   I7196 M  Czech_EBA
                   I7197 M   Czech_MN
                   I7198 F  Czech_EBA
                   I7199 M  Czech_EBA
                   I7200 F  Czech_EBA
                   I7201 F  Czech_EBA
                   I7202 M  Czech_EBA
                   I7203 M  Czech_EBA
                   I7205 M Czech_BellBeaker
                   I7211 F Czech_BellBeaker
                   I7213 F Czech_BellBeaker
                   I7249 M Czech_BellBeaker
                   I7251 M Czech_BellBeaker
                   I7250 F Czech_BellBeaker
                   I7269 M Czech_BellBeaker
                   I7272 M Czech_Eneolithic
                   I7270 F Czech_BellBeaker
                   I7275 M Czech_BellBeaker
                   I7276 M Czech_BellBeaker
                   I7279 M Czech_CordedWare
                   I7280 M Czech_CordedWare
                   I7281 F Czech_BellBeaker
                   I7286 M Czech_BellBeaker
                   I7287 M Czech_BellBeaker
                   I7288 M Czech_BellBeaker
                   I7290 F Czech_BellBeaker
            S_Czech-2.DG M   Czech.DG
              RISE566.SG M Czech_BellBeaker_dup.I4145.SG
              RISE567.SG F Czech_BellBeaker_dup.I4136.SG
              RISE568.SG F Czech_EarlySlav.SG
              RISE569.SG F Czech_EarlySlav_dup.I4137.SG
              RISE577.SG F Czech_Unetice_EBA_dup.I4139.SG
              RISE586.SG F Czech_EBA_Unetice_dup.I4130
             Vestonice16 M Czech_Vestonice16
               Pavlov1_d M Czech_Pavlov1
           Vestonice13_d M Czech_Vestonice13
           Vestonice15_d M Czech_Vestonice15
           Vestonice43_d M Czech_Vestonice43
           Vestonice14_d M Czech_Vestonice14_lc
                   I4145 M Czech_BellBeaker
                   I4136 F Czech_BellBeaker
                   I4139 F Czech_EBA_Starounetice_dup.I4139
                   I5037 M Czech_EBA_Protounetice
                   I5042 M Czech_EBA_Protounetice
                   I4141 M Czech_EBA_Unetice
                   I4130 F Czech_EBA_Unetice
                   I4945 F Czech_BellBeaker
                   I4946 F Czech_BellBeaker
                   I4947 M Czech_BellBeaker_lc
                   I4884 M  Czech_EBA
                   I4885 M Czech_BellBeaker
                   I4886 M Czech_BellBeaker
                   I4887 M Czech_BellBeaker_brother.I4888
                   I4888 M Czech_BellBeaker
                   I4889 M Czech_BellBeaker
                   I4890 M Czech_BellBeaker
                   I4891 M Czech_BellBeaker
                   I4892 F  Czech_EBA
                   I4893 M    Czech_N
                   I4894 F    Czech_N
                   I4895 M Czech_BellBeaker
                   I4896 F Czech_BellBeaker
                   I5514 M Czech_BellBeaker
                   I5666 M Czech_BellBeaker
                   I6468 F Czech_BellBeaker
                   I6476 F Czech_BellBeaker
                   I6480 M Czech_BellBeaker
                   I6677 M Czech_Baalberge
                   I7949 M Czech_BA_Veterov_1
                DA111.SG M Czech_IA_Hallstatt.SG
                DA112.SG F Czech_IA_Hallstatt.SG
               kol2BE.SG F Czech_Megalithic.SG
        kol6-ALL_DATA.SG F Czech_Megalithic.SG


### Only the individuals starting with <code>NA</code> are modern ones. Reference: https://www.coriell.org/0/Sections/Search/Panel_Detail.aspx?Ref=HD29&Product=HDP


```bash
%%bash
grep Czech v44.3_HO_public.ind|grep NA|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Czech.id && wc -l Czech.id
```

    10 Czech.id



```bash
%%bash
grep Hungarian v44.3_HO_public.ind
```

                hungary3 M  Hungarian
                hungary6 M  Hungarian
             HungarianC5 M  Hungarian
                hungary7 M  Hungarian
             HungarianE5 M  Hungarian
             HungarianH3 M  Hungarian
               hungary15 M  Hungarian
               hungary20 M  Hungarian
                hungary2 M  Hungarian
             HungarianD1 M  Hungarian
                 NA15202 M  Hungarian
                 NA15203 F  Hungarian
                 NA15204 F  Hungarian
                 NA15205 M  Hungarian
                 NA15206 M  Hungarian
                 NA15199 F  Hungarian
                 NA15207 F  Hungarian
                 NA15200 M  Hungarian
                 NA15208 M  Hungarian
                 NA15201 F  Hungarian
        S_Hungarian-1.DG F Hungarian.DG
        S_Hungarian-2.DG M Hungarian.DG


### All this samples are modern ones


```bash
%%bash
grep Hungarian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Hungarian.id && wc -l Hungarian.id
```

    20 Hungarian.id



```bash
%%bash
grep Slovakian v44.3_HO_public.ind # No Slovakians present
```


```bash
%%bash
grep Dutch v44.3_HO_public.ind # No Dutch present
```


```bash
%%bash
grep Croatian v44.3_HO_public.ind
```

                   CRO53 M   Croatian
                  CRO103 M   Croatian
                  CRO107 M   Croatian
                   CRO47 M   Croatian
                   CRO41 M   Croatian
                   CRO66 M   Croatian
                  CRO153 M   Croatian
                   CRO31 M   Croatian
                   CRO93 M   Croatian
                   CRO48 M   Croatian



```bash
%%bash
grep Croatian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Croatian.id && wc -l Croatian.id
```

    10 Croatian.id



```bash
%%bash
grep Greek v44.3_HO_public.ind
```

                  TLA010 M  Greek.WGA
                  TLA011 M Greek_outlier.WGA
                  TLA012 M  Greek.WGA
                  TLA013 M Greek_outlier.WGA
                  TLA015 M  Greek.WGA
                  TLA017 M  Greek.WGA
                  TLA018 M  Greek.WGA
                  TLA019 M  Greek.WGA
                  TLA020 M  Greek.WGA
                  TLA021 M  Greek.WGA
                  TLA022 M  Greek.WGA
                  TLA023 M  Greek.WGA
                  TLA024 M  Greek.WGA
                  TLA025 M  Greek.WGA
                  TLA026 M  Greek.WGA
                  TLA027 M  Greek.WGA
                  TLA028 M  Greek.WGA
                  TLA029 M  Greek.WGA
                 NA17373 F      Greek
                 NA17374 M      Greek
                 NA17375 F      Greek
                 NA17376 M      Greek
                 NA17377 M      Greek
                 NA17370 M Ignore_Greek
                 NA17371 F Ignore_Greek
                 NA17372 F      Greek
          GREEKGRALPOP18 M      Greek
          GREEKGRALPOP13 M      Greek
          GREEKGRALPOP15 M      Greek
           GREEKGRALPOP5 M      Greek
           GREEKGRALPOP9 M      Greek
           GREEKGRALPOP4 M      Greek
          GREEKGRALPOP17 M      Greek
          GREEKGRALPOP16 M      Greek
          GREEKGRALPOP11 M      Greek
          GREEKGRALPOP10 M      Greek
           GREEKGRALPOP3 M      Greek
           GREEKGRALPOP8 M      Greek
          GREEKGRALPOP12 M      Greek
          GREEKGRALPOP14 M      Greek
            S_Greek-1.DG M Greek_1.DG
            S_Greek-2.DG M Greek_2.DG
                   I8340 F Spain_Greek_oLocal
                   I8341 M Spain_Greek_oLocal
                   I8344 M Spain_Greek_oLocal
                   I8209 M Spain_Greek_oLocal
                   I8210 M Spain_Greek_oLocal
                   I8211 M Spain_Greek_oLocal
                   I8212 M Spain_Greek_oLocal
                   I8213 F Spain_Greek_oLocal_lc
                   I8214 F Spain_Greek_oLocal
                   I8215 F Spain_Greek_oAegean



```bash
%%bash
grep Greek v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|grep -v Spain|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Greek.id && wc -l Greek.id
```

    36 Greek.id



```bash
%%bash
grep Bulgarian v44.3_HO_public.ind
```

             BulgarianD6 M  Bulgarian
             BulgarianA4 M  Bulgarian
             BulgarianE2 M  Bulgarian
             BulgarianB4 M  Bulgarian
             BulgarianA1 M  Bulgarian
             BulgarianB1 M  Bulgarian
             BulgarianC1 M  Bulgarian
             BulgarianF1 M  Bulgarian
             BulgarianH2 F  Bulgarian
             BulgarianF2 F  Bulgarian
        S_Bulgarian-1.DG M Bulgarian.DG
        S_Bulgarian-2.DG M Bulgarian.DG



```bash
%%bash
grep Bulgarian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Bulgarian.id && wc -l Bulgarian.id
```

    10 Bulgarian.id



```bash
%%bash
grep Lithuanian v44.3_HO_public.ind
```

            LithuanianF1 F Lithuanian
              lithuania3 M Lithuanian
             lithuania10 M Lithuanian
              lithuania9 M Lithuanian
            LithuanianA1 M Lithuanian
            LithuanianE2 F Lithuanian
              lithuania1 M Lithuanian
              lithuania8 M Lithuanian
              lithuania2 M Lithuanian
            LithuanianD1 F Lithuanian



```bash
%%bash
grep Lithuanian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Lithuanian.id && wc -l Lithuanian.id
```

    10 Lithuanian.id



```bash
%%bash
grep Latvian v44.3_HO_public.ind ## No Latvians
```


```bash
%%bash
grep Estonian v44.3_HO_public.ind
```

                  Est393 M   Estonian
                  Est375 M   Estonian
                  Est380 M   Estonian
                  Est391 M   Estonian
                  Est377 M   Estonian
                  Est372 M   Estonian
                  Est358 M   Estonian
                  Est400 M   Estonian
                  Est397 M   Estonian
                  Est394 M   Estonian
         S_Estonian-2.DG M Estonian.DG
         S_Estonian-1.DG M Estonian.DG



```bash
%%bash
grep Estonian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Estonian.id && wc -l Estonian.id
```

    10 Estonian.id



```bash
%%bash
grep Ukrainian v44.3_HO_public.ind
```

                UKR-1283 M  Ukrainian
                UKR-1291 M  Ukrainian
                UKR-1292 M  Ukrainian
                UKR-1377 M  Ukrainian
                UKR-1399 M Ukrainian_North
                UKR-1903 M Ukrainian_North
                UKR-1909 M Ukrainian_North
                UKR-1913 M Ukrainian_North
                UKR-1951 M Ukrainian_North
                UKR-1978 M Ukrainian_North
                UKR-1992 M Ukrainian_North
                UKR-2021 M Ukrainian_North
               UkrBel618 M  Ukrainian
                UkrLv240 M  Ukrainian
               UkrBel620 M  Ukrainian
               UkrBel622 M  Ukrainian
               UkrBel733 F  Ukrainian
               UkrBel736 F  Ukrainian
                UkrLv228 M  Ukrainian
               UkrBel614 M  Ukrainian
                UkrLv237 M  Ukrainian



```bash
%%bash
grep Ukrainian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|grep -v North|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Ukrainian.id && wc -l Ukrainian.id
```

    13 Ukrainian.id



```bash
%%bash
grep Ukrainian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|grep North|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Ukrainian_north.id && wc -l Ukrainian_north.id
```

    8 Ukrainian_north.id



```bash
%%bash
grep Belarusian v44.3_HO_public.ind
```

                  bel43s M Belarusian
                  bel30s M Belarusian
                  bel72c M Belarusian
                  bel93c M Belarusian
                 bel110c M Belarusian
                   bel8s M Belarusian
          belarusian23vp M Belarusian
                  bel23s M Belarusian
          belarusian47zp M Belarusian
                  bel82s M Belarusian



```bash
%%bash
grep Belarusian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Belarusian.id && wc -l Belarusian.id
```

    10 Belarusian.id



```bash
%%bash
grep Russian v44.3_HO_public.ind|grep -v [Ii]gnore && grep Russian v44.3_HO_public.ind|grep -v [Ii]gnore|wc -l 
```

                Rakr-203 M Russian_Archangelsk_Krasnoborsky
                Rakr-205 M Russian_Archangelsk_Krasnoborsky
                Rakr-237 M Russian_Archangelsk_Krasnoborsky
                Rakr-248 M Russian_Archangelsk_Krasnoborsky
                Rakr-341 M Russian_Archangelsk_Krasnoborsky
                Rakr-345 M Russian_Archangelsk_Krasnoborsky
             Rakrlsh-002 M Russian_Archangelsk_Leshukonsky
             Rakrlsh-140 F Russian_Archangelsk_Leshukonsky
             Rakrlsh-143 M Russian_Archangelsk_Leshukonsky
             Rakrlsh-144 M Russian_Archangelsk_Leshukonsky
             Rakrlsh-149 M Russian_Archangelsk_Leshukonsky
                Rbgp-200 M    Russian
                Rbgp-201 M    Russian
                Rbgp-203 M    Russian
                Rbgp-205 U    Russian
                 Rkbo-12 F    Russian
                 Rkbo-16 M    Russian
                 Rkbo-44 M    Russian
                 Rkbo-58 M    Russian
                Rksh-402 M    Russian
                Rksh-405 M    Russian
                Rksh-407 M    Russian
                Rksh-412 F    Russian
                Rkuch-03 M    Russian
                Rkuch-05 M    Russian
                Rkuch-53 M    Russian
                Rkuch-58 M    Russian
                Rorl-102 F    Russian
                Rorl-110 F    Russian
                Rorl-114 M    Russian
                Rorl-155 M    Russian
                RPin-114 M Russian_Archangelsk_Pinezhsky
                RPin-123 M Russian_Archangelsk_Pinezhsky
                RPin-143 M Russian_Archangelsk_Pinezhsky
                RPin-145 M Russian_Archangelsk_Pinezhsky
                RPin-151 M Russian_Archangelsk_Pinezhsky
                 Rps-002 M    Russian
                 Rps-004 M    Russian
                 Rps-006 M    Russian
                 Rps-012 M    Russian
                 Rps-090 M    Russian
                 Rps-091 M    Russian
                 Rps-098 M    Russian
                 Rrzm-08 M    Russian
                 Rrzm-10 F    Russian
                 Rrzm-13 M    Russian
                 Rrzm-16 M    Russian
                 Rrzm-83 M    Russian
                  Rrzs-3 F    Russian
                  Rrzs-7 M    Russian
                 Rrzs-11 M    Russian
                 Rrzs-32 F    Russian
                 Rrzs-58 M    Russian
                 Rrzs-66 M    Russian
                 Rrzs-88 M    Russian
                 Rsm-103 M    Russian
                 Rsm-109 M    Russian
                 Rsm-166 M    Russian
                 Rsm-171 M    Russian
                 Rsm-176 M    Russian
                 Rsm-179 M    Russian
                 Rsm-181 M    Russian
                RYAR-173 M    Russian
                RYAR-223 M    Russian
                RYAR-232 M    Russian
               HGDP00879 M    Russian
               HGDP00880 M    Russian
               HGDP00882 M    Russian
               HGDP00883 M    Russian
               HGDP00884 F    Russian
               HGDP00887 M    Russian
               HGDP00888 M    Russian
               HGDP00889 F    Russian
               HGDP00890 M    Russian
               HGDP00891 M    Russian
               HGDP00892 M    Russian
               HGDP00893 M    Russian
               HGDP00894 M    Russian
               HGDP00895 M    Russian
               HGDP00896 M    Russian
               HGDP00897 M    Russian
               HGDP00898 F    Russian
               HGDP00899 F    Russian
               HGDP00900 M    Russian
               HGDP00901 F    Russian
               HGDP00902 F    Russian
               HGDP00903 F    Russian
          S_Russian-2.DG F Russian.DG
          S_Russian-1.DG M Russian.DG
           HGDP00879.SDG M Russian.SDG
           HGDP00902.SDG F Russian.SDG
           HGDP00880.SDG M Russian.SDG
           HGDP00900.SDG M Russian.SDG
           HGDP00898.SDG F Russian.SDG
           HGDP00899.SDG F Russian.SDG
           HGDP00901.SDG F Russian.SDG
           HGDP00886.SDG U Russian.SDG
           HGDP00895.SDG M Russian.SDG
           HGDP00891.SDG M Russian.SDG
           HGDP00883.SDG M Russian.SDG
           HGDP00882.SDG M Russian.SDG
           HGDP00884.SDG F Russian.SDG
           HGDP00889.SDG F Russian.SDG
           HGDP00896.SDG M Russian.SDG
           HGDP00890.SDG M Russian.SDG
           HGDP00893.SDG M Russian.SDG
           HGDP00897.SDG M Russian.SDG
           HGDP00885.SDG U Russian.SDG
           HGDP00892.SDG M Russian.SDG
           HGDP00888.SDG M Russian.SDG
           HGDP00894.SDG M Russian.SDG
           HGDP00903.SDG F Russian.SDG
           HGDP00887.SDG M Russian.SDG
    113



```bash
%%bash
grep Russian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep -v Archangelsk|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Russian.id && wc -l Russian.id
```

    71 Russian.id



```bash
%%bash
ls *.id
```

    Algerian.id    Finnish.id	 Lebanese_Christian.id	Sardinian.id
    Assyrian.id    French.id	 Lebanese_Muslim.id	Saudi.id
    Basque.id      French_south.id	 Lebanese_unknown.id	Scottish.id
    Bedouin.id     Greek.id		 Libyan.id		Spanish.id
    Belarusian.id  Hungarian.id	 Lithuanian.id		Spanish_north.id
    Bulgarian.id   Icelandic.id	 Maltese.id		Syrian.id
    Croatian.id    Iranian_gulf.id	 Moroccan.id		Tunisian.id
    Czech.id       Iranian.id	 Norwegian.id		Tuscan.id
    Druze.id       Italian_north.id  Palestinian.id		Ukrainian.id
    Egyptian.id    Italian_south.id  Pathan.id		Ukrainian_north.id
    English.id     Jordanian.id	 Polish.id
    Estonian.id    Kurd.id		 Russian.id



```bash
%%bash
grep Russian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep Archangelsk|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Russian_nordic.id && wc -l Russian_nordic.id
```

    16 Russian_nordic.id



```bash
%%bash
grep Albanian v44.3_HO_public.ind
```

                  ALB191 F   Albanian
                  ALB213 F   Albanian
                  ALB202 M   Albanian
                  ALB212 F   Albanian
                  ALB220 F   Albanian
                  ALB230 F   Albanian
         S_Albanian-1.DG F Albanian.DG



```bash
%%bash
grep Albanian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Albanian.id && wc -l Albanian.id
```

    6 Albanian.id



```bash
%%bash
grep Romanian v44.3_HO_public.ind|grep -v [Jjew]
```

                    A306 F   Romanian
                    A325 F   Romanian
                    A343 F   Romanian
                    A362 F   Romanian
                    A374 F   Romanian
                    G408 F   Romanian
                    G421 F   Romanian
                    G428 F   Romanian
                    G429 F   Romanian
                    G434 F   Romanian



```bash
%%bash
grep Romanian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Romanian.id && wc -l Romanian.id
```

    10 Romanian.id



```bash
%%bash
grep Moldavian v44.3_HO_public.ind|grep -v [Jjew]
```

                 MOL-005 M  Moldavian
                 MOL-008 M  Moldavian
                 MOL-015 M  Moldavian
                 MOL-024 M  Moldavian
                 MOL-058 M  Moldavian
                 MOL-064 F  Moldavian
                 MOL-065 F  Moldavian
                 MOL-066 F  Moldavian
                 MOL-067 F  Moldavian
                 MOL-069 M  Moldavian



```bash
%%bash
grep Moldavian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Moldavian.id && wc -l Moldavian.id
```

    10 Moldavian.id



```bash
%%bash
grep Mordovian v44.3_HO_public.ind|grep -v [Jjew]
```

                 MOE-001 M  Mordovian
                 MOE-002 M  Mordovian
                 MOE-010 M  Mordovian
                 MOE-014 M  Mordovian
                 MOE-015 M  Mordovian
                 MOE-020 M  Mordovian
                 MOE-025 M  Mordovian
                 MOE-036 M  Mordovian
                 MOE-043 M  Mordovian
                 MOE-045 M  Mordovian
                 MOE-433 M  Mordovian
                 MOE-445 M  Mordovian
                 MOE-450 M  Mordovian
                 MOE-451 M  Mordovian
                 MOE-452 M  Mordovian
                 MOE-455 M  Mordovian
                 MOE-475 M  Mordovian
                 MOE-485 M  Mordovian
                 MOE-491 M  Mordovian
                 MOE-492 M  Mordovian
                 MOE-495 M  Mordovian
                 MOE-497 M  Mordovian
            Mordovians27 M  Mordovian
            Mordovians28 M  Mordovian
            Mordovians30 M  Mordovian
             Mordovians4 M  Mordovian
            Mordovians31 M  Mordovian
            Mordovians17 M  Mordovian
             Mordovians1 M  Mordovian
            Mordovians22 M  Mordovian
            Mordovians32 M  Mordovian
             Mordovians5 M  Mordovian



```bash
%%bash
grep Mordovian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Mordovian.id && wc -l Mordovian.id
```

    32 Mordovian.id


## Turkic people


```bash
%%bash
grep Turkish v44.3_HO_public.ind|grep -v [Ii]gnore|grep -v [Jj]ew
```

            Kayseri24392 M    Turkish
          Balikesir16675 F Turkish_Balikesir
            Turkish7BA57 M    Turkish
          Balikesir16790 M Turkish_Balikesir
            Turkish8BA62 M    Turkish
           Istanbul20010 M    Turkish
            Turkish4BA57 M    Turkish
              Adana23114 F    Turkish
            Trabzon21177 M    Turkish
              Aydin18784 M    Turkish
           Istanbul19810 F    Turkish
            Trabzon21557 M    Turkish
            Trabzon21534 M    Turkish
            Turkish9BA57 M    Turkish
            Kayseri24075 M    Turkish
              Adana23108 M    Turkish
          Balikesir16887 F Turkish_Balikesir
            Kayseri24266 F    Turkish
          Balikesir16653 M Turkish_Balikesir
           Istanbul25095 M    Turkish
              Aydin18112 M    Turkish
            Kayseri24402 F    Turkish
           Istanbul25081 M    Turkish
              Aydin18596 M    Turkish
           Istanbul17778 M    Turkish
           Istanbul15781 M    Turkish
              Adana23136 M    Turkish
              Adana23113 M    Turkish
            Trabzon21575 M    Turkish
              Aydin18636 F    Turkish
           Istanbul20040 F    Turkish
              Adana23144 F    Turkish
            Kayseri23967 M    Turkish
            Trabzon21174 M    Turkish
              Aydin18873 M    Turkish
            Kayseri23892 M    Turkish
           Istanbul19185 M    Turkish
           Istanbul25098 M    Turkish
              Adana23133 F    Turkish
            Trabzon21515 M    Turkish
            Kayseri24032 M    Turkish
            Kayseri23549 M    Turkish
            Trabzon21544 M    Turkish
              Adana23150 M    Turkish
            Kayseri23271 F    Turkish
            Trabzon21645 M    Turkish
              Adana23147 M    Turkish
              Aydin18483 F    Turkish
          Balikesir17006 M Turkish_Balikesir
              Aydin18419 F    Turkish
            Kayseri24276 F    Turkish
           Istanbul19708 F    Turkish
            Trabzon21173 M    Turkish
              Adana23112 M    Turkish
              Adana23117 M    Turkish
          S_Turkish-2.DG F Turkish.DG
          S_Turkish-1.DG M Turkish.DG


### Therefore, let's sort out the geographic location of these Turks


```bash
%%bash
grep Turkish v44.3_HO_public.ind|grep -v [Ii]gnore|grep -v [Jj]ew|sed -r 's/^\s+//'|sed -r 's/\s+/ /g'|\
cut -d ' ' -f 1|sed -r 's/[0-9].+$//g'|sort|uniq
```

    Adana
    Aydin
    Balikesir
    Istanbul
    Kayseri
    S_Turkish-
    Trabzon
    Turkish


#### Western coastal: Balikesir, Aydin
#### Southern coastal: Adana
#### Northeastern coastal: Trabzon
#### Central: Kayseri
#### mixed: Istanbul, S_Turkish-, Turkish


```bash
%%bash
grep Turkish v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep -E '(Balikesir|Aydin)'|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Turkish_west.id && wc -l Turkish_west.id
```

    12 Turkish_west.id



```bash
%%bash
grep Turkish v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep -E '(Istanbul|S_Turkish-|Turkish[0-9])'|sed -r 's/^\s+//'|cut -d ' ' -f 1 > \
Turkish_mixed.id && wc -l Turkish_mixed.id
```

    14 Turkish_mixed.id



```bash
%%bash
grep Turkish v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep Trabzon|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Turkish_northeast.id && wc -l Turkish_northeast.id
```

    9 Turkish_northeast.id



```bash
%%bash
grep Turkish v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep Adana|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Turkish_south.id && wc -l Turkish_south.id
```

    10 Turkish_south.id



```bash
%%bash
grep Turkish v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep Kayseri|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Turkish_central.id && wc -l Turkish_central.id
```

    10 Turkish_central.id



```bash
%%bash
grep Turkish v44.3_HO_public.ind|grep -v [Ii]gnore|grep -v [Jj]ew|grep -v DG$|wc -l #check the sum, correct!
```

    55


### It's possible to sub-group these Turkish individuals from their iid(city present). Since Anatolia is a relatively large geographic region. We would like to perform it later in our IBD analysis


```bash
%%bash
grep Tatar v44.3_HO_public.ind
```

                 STA-003 M Tatar_Siberian
                 STA-004 M Tatar_Siberian
                 STA-005 M Ignore_Tatar_Siberian(PCA_outlier)
                 STA-006 M Tatar_Siberian
                 STA-112 M Tatar_Siberian
                 STA-116 M Tatar_Siberian
                 STA-120 M Tatar_Siberian
                 STA-122 M Tatar_Siberian
                 STA-126 M Tatar_Siberian
                 STA-128 M Tatar_Siberian
                 STA-205 M Tatar_Siberian
                 STA-211 M Tatar_Siberian
                 STA-212 M Tatar_Siberian
                 STA-237 M Tatar_Siberian
                 STA-265 M Tatar_Siberian
                 STA-297 M Tatar_Siberian
                 STA-300 M Tatar_Siberian
                 STA-304 M Tatar_Siberian
                 STA-306 M Tatar_Siberian
                 STA-309 M Tatar_Siberian
                 STA-357 M Tatar_Siberian_Zabolotniye
                 STA-362 M Tatar_Siberian_Zabolotniye
                 STA-366 M Tatar_Siberian_Zabolotniye
                 STA-434 M Tatar_Siberian_Zabolotniye
                 STA-435 M Tatar_Siberian_Zabolotniye
                 TTR-086 M Tatar_Kazan
                 TTR-094 M Tatar_Kazan
                 TTR-097 M Tatar_Kazan
                 TTR-201 M Tatar_Kazan
                 TTR-217 M Tatar_Kazan
                 TTR-241 M Tatar_Mishar
                 TTR-244 M Tatar_Kazan
                 TTR-245 M Tatar_Kazan
                 TTR-249 M Tatar_Kazan
                 TTR-250 M Tatar_Kazan
                 TTR-271 M Tatar_Mishar
                 TTR-272 M Tatar_Mishar
                 TTR-330 M Tatar_Kazan
                 TTR-356 M Tatar_Mishar
                 TTR-359 M Tatar_Mishar
                 TTR-362 M Tatar_Mishar
                 TTR-436 M Tatar_Mishar
                 TTR-460 M Tatar_Mishar
                 TTR-462 M Tatar_Mishar
                 TTR-464 M Tatar_Mishar
                 TTR-493 M Tatar_Kazan
                 TTR-501 M Tatar_Kazan
                 TTR-514 M Tatar_Kazan
    IrtyshBarabinskTatars1.SG U Tatar_Irtysh_Barabinsk.SG
    IrtyshBarabinskTatars2.SG U Tatar_Irtysh_Barabinsk.SG
         TomskTatars1.SG U Tatar_Tomsk.SG
         TomskTatars2.SG U Tatar_Tomsk.SG
         VolgaTatars1.SG U Tatar_Volga.SG
         VolgaTatars2.SG U Tatar_Volga.SG


### We would also need sub-grouping for Tatars


```bash
%%bash
grep Tatar v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep Kazan|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Tatar_Kazan.id && wc -l Tatar_Kazan.id
```

    13 Tatar_Kazan.id


### Mishar is a subgroup of Volga Tatars


```bash
%%bash
grep Tatar v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep Mishar|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Tatar_Volga.id && wc -l Tatar_Volga.id
```

    10 Tatar_Volga.id



```bash
%%bash
grep Tatar v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep Siberian|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Tatar_Siberian.id && wc -l Tatar_Siberian.id
```

    24 Tatar_Siberian.id



```bash
%%bash
grep Tatar v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|wc -l #Correct!
```

    47



```bash
%%bash
grep Bashkir v44.3_HO_public.ind
```

                 BAS-005 M    Bashkir
                 BAS-006 M    Bashkir
                 BAS-008 M    Bashkir
                 BAS-014 M    Bashkir
                 BAS-017 M    Bashkir
                 BAS-021 M    Bashkir
                 BAS-029 M    Bashkir
                 BAS-031 M    Bashkir
                 BAS-033 M    Bashkir
                 BAS-034 M    Bashkir
                 BAS-042 M    Bashkir
                 BAS-045 M    Bashkir
                 BAS-046 M    Bashkir
                 BAS-060 M    Bashkir
                 BAS-062 M    Bashkir
                 BAS-091 M    Bashkir
                 BAS-094 M    Bashkir
                 BAS-096 M    Bashkir
                 BAS-105 M    Bashkir
                 BAS-111 M    Bashkir
                 BAS-120 M    Bashkir
                 BAS-121 M    Bashkir
                 BAS-125 M    Bashkir
                 BAS-135 M    Bashkir
                 BAS-150 M    Bashkir
                 BAS-153 M    Bashkir
                 BAS-156 M    Bashkir
                 BAS-164 M    Bashkir
                 BAS-600 M    Bashkir
                 BAS-622 M    Bashkir
                 BAS-652 M    Bashkir
                 BAS-655 M    Bashkir
                 BAS-661 M    Bashkir
                 BAS-663 M    Bashkir
                 BAS-669 M    Bashkir
                 BAS-670 M    Bashkir
                 BAS-671 M    Bashkir
                 BAS-672 M    Bashkir
                 BAS-683 M    Bashkir
                 BAS-811 M    Bashkir
                 BAS-813 M    Bashkir
                 BAS-822 M    Bashkir
                 BAS-825 M    Bashkir
                 BAS-831 M    Bashkir
                 BAS-833 M    Bashkir
                 BAS-834 M    Bashkir
                 BAS-849 M    Bashkir
                BAS-1392 M    Bashkir
                BAS-1393 M    Bashkir
                BAS-1394 M    Bashkir
                BAS-1396 M    Bashkir
                BAS-1398 M    Bashkir
                BAS-1400 M    Bashkir
            Bashkirs1.SG U Bashkir.SG
            Bashkirs2.SG U Bashkir_o.SG
            Bashkirs3.SG U Bashkir.SG



```bash
%%bash
grep Bashkir v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Bashkir.id && wc -l Bashkir.id
```

    53 Bashkir.id



```bash
%%bash
grep Chuvash v44.3_HO_public.ind
```

                 Ttr-473 M    Chuvash
                 Ttr-474 F    Chuvash
                 Ttr-481 M    Chuvash
                 Ttr-507 M    Chuvash
                 Ttr-568 M    Chuvash
                 Ttr-569 M    Chuvash
               Chuvash33 M    Chuvash
               Chuvash37 M    Chuvash
               Chuvash13 M    Chuvash
               Chuvash20 M    Chuvash
               Chuvash22 M    Chuvash
               Chuvash26 M    Chuvash
               Chuvash24 M    Chuvash
               Chuvash29 M    Chuvash
               Chuvash25 M    Chuvash
               Chuvash31 M    Chuvash



```bash
%%bash
grep Chuvash v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Chuvash.id && wc -l Chuvash.id
```

    16 Chuvash.id



```bash
%%bash
grep Turkmen v44.3_HO_public.ind|grep -v 'I[0-9]' #exclude Ancient indivs
```

          UZB178_turkmen M    Turkmen
          UZB180_turkmen M    Turkmen
          UZB101_turkmen M Turkmen_outlier
          UZB102_turkmen M    Turkmen
          UZB105_turkmen M    Turkmen
          UZB111_turkmen M    Turkmen
          UZB150_turkmen M    Turkmen
          DA379_final.SG M Turkmenistan_C_Namazga_1d.rel.DA380.SG
          DA380_final.SG F Turkmenistan_C_Namazga.SG
          DA381_final.SG M Turkmenistan_C_Namazga.SG
          DA382_final.SG M Turkmenistan_IA.SG
          DA383_final.SG F Turkmenistan_C_Namazga_o.SG
            Turkmens1.SG U Turkmen.SG
            Turkmens2.SG U Turkmen.SG



```bash
%%bash
grep Turkmen v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|grep -v 'I[0-9]'|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Turkmen.id && wc -l Turkmen.id
```

    6 Turkmen.id



```bash
%%bash
grep Uzbek v44.3_HO_public.ind|grep -v 'I[0-9]' #exclude Ancient indivs
```

                 UZB-229 M      Uzbek
                 UZB-232 M      Uzbek
                 UZB-234 M      Uzbek
                 UZB-235 M      Uzbek
                 UZB-236 M      Uzbek
                 UZB-237 M      Uzbek
                 UZB-238 M      Uzbek
                 UZB-244 M      Uzbek
                 UZB-246 M      Uzbek
                 UZB-301 M      Uzbek
                 UZB-302 M      Uzbek
                 UZB-303 M      Uzbek
                 UZB-304 M      Uzbek
                 UZB-305 M      Uzbek
                 UZB-306 M      Uzbek
                 UZB-307 M      Uzbek
                 UZB-308 M      Uzbek
                  Tash02 M Uzbek_outlier.WGA
                   usb24 M      Uzbek
                   usb25 M      Uzbek
                   usb35 F      Uzbek
                   usb40 M      Uzbek
                    usb6 F      Uzbek
                   usb64 F      Uzbek
                   usb72 M      Uzbek
                    usb2 F      Uzbek
                   usb78 M      Uzbek
                    usb8 F      Uzbek
              Uzbeks1.SG U   Uzbek.SG
              Uzbeks2.SG U   Uzbek.SG
              Uzbeks3.SG U   Uzbek.SG



```bash
%%bash
grep Uzbek v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|\
grep -v outlier|grep -v 'I[0-9]'|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Uzbek.id && wc -l
```

## Population from Caucasus


```bash
%%bash
grep Armenian v44.3_HO_public.ind
```

                  ARK-59 M Armenian_Hemsheni
                  ARK-72 M Armenian_Hemsheni
                  ARK-76 M Armenian_Hemsheni
                  ARK-78 M Armenian_Hemsheni
                  ARK-82 M Armenian_Hemsheni
                  ARK-84 M Armenian_Hemsheni
                  ARK-88 M Armenian_Hemsheni
                  ARK-89 M Armenian_Hemsheni
                  ARM012 M Armenian.WGA
                  ARM013 M Armenian.WGA
                  ARM014 M Armenian.WGA
              armenia176 M   Armenian
              armenia191 M   Armenian
               armenia86 M   Armenian
              armenia279 M   Armenian
               armenia91 M   Armenian
              armenia293 M   Armenian
              armenia102 M   Armenian
              armenia106 M   Armenian
              armenia139 M   Armenian
              armenia162 M   Armenian
         S_Armenian-1.DG M Armenian.DG
         S_Armenian-2.DG M Armenian.DG


### Hemsheni is a distinct population in Northeastern Anatolia with Armenian origin. It should be regarded as a seperate ethnic group other than Armenians


```bash
%%bash
grep Hemsheni v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Armenian_Hemsheni.id && wc -l Armenian_Hemsheni.id
```

    8 Armenian_Hemsheni.id



```bash
%%bash
grep Armenian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep -v Hemsheni|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Armenian.id && wc -l Armenian.id
```

    13 Armenian.id



```bash
%%bash
grep Georgian v44.3_HO_public.ind
```

                 GEO-002 M   Georgian
                 GEO-005 M   Georgian
                 GEO-010 M   Georgian
                 GEO-015 M   Georgian
                 GEO-020 M   Georgian
                 GEO-028 M   Georgian
                 GEO-031 M   Georgian
                 GEO-032 M   Georgian
                 GEO-039 M   Georgian
                 GEO-051 M   Georgian
                 GEO-061 M   Georgian
                 GEO-082 M   Georgian
                  GEO001 M   Georgian
                  GEO013 M Georgian.WGA
                  GEO014 M Georgian.WGA
                    mg43 M   Georgian
                    mg47 M   Georgian
                    mg22 M   Georgian
                    mg49 M   Georgian
                    mg23 M   Georgian
                    mg62 M   Georgian
                    mg27 M   Georgian
                    mg31 M   Georgian
                    mg34 M   Georgian
                    mg40 M   Georgian
         GeorgianJew1607 F Jew_Georgian
         GeorgianJew1671 M Jew_Georgian
         GeorgianJew1577 M Jew_Georgian
         GeorgianJew1971 F Jew_Georgian
         GeorgianJew1594 M Jew_Georgian
         GeorgianJew1961 M Ignore_Jew_Georgian
         GeorgianJew1654 F Jew_Georgian
         GeorgianJew1972 F Ignore_Jew_Georgian
         GeorgianJew1883 M Jew_Georgian
         S_Georgian-2.DG M Georgian.DG
         S_Georgian-1.DG M Georgian.DG



```bash
%%bash
grep Georgian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Georgian.id && wc -l Georgian.id
```

    25 Georgian.id



```bash
%%bash
grep Azeri v44.3_HO_public.ind
```

                AZR-0853 M      Azeri
                AZR-0860 M      Azeri
                AZR-0864 M      Azeri
                AZR-0865 M      Azeri
                AZR-0866 M      Azeri
                AZR-0868 M      Azeri
                AZR-0869 M      Azeri
                AZR-0870 M      Azeri
                AZR-1010 M      Azeri
                AZR-1012 M      Azeri
                AZR-1013 M      Azeri
                AZR-1017 M      Azeri
                AZR-1018 M      Azeri
                AZR-1037 M      Azeri
                AZR-1039 M      Azeri
                AZR-1054 M      Azeri
                AZR-1058 M      Azeri
                  Baku01 M  Azeri.WGA
                  Baku02 M  Azeri.WGA
                  Baku03 M  Azeri.WGA



```bash
%%bash
grep Azeri v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Azeri.id && wc -l Azeri.id
```

    20 Azeri.id



```bash
%%bash
grep Chechen v44.3_HO_public.ind
```

                   ch126 M    Chechen
                    ch16 M    Chechen
                   ch174 M    Chechen
                   ch179 M    Chechen
                   ch193 M    Chechen
                    ch21 M    Chechen
                     ch3 M    Chechen
                   ch113 M    Chechen
                    ch31 M    Chechen
          S_Chechen-1.DG M Chechen.DG



```bash
%%bash
grep Chechen v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Chechen.id && wc -l Chechen.id
```

    9 Chechen.id



```bash
%%bash
grep Avar v44.3_HO_public.ind
```

                 DAG-436 M Avar_outlier1
                 DAG-448 M       Avar
                 DAG-461 M       Avar
                 DAG-474 M       Avar
                 DAG-477 M       Avar
                 DAG-485 M       Avar
                 DAG-492 M       Avar
                 DAG-506 M       Avar
                 DAG-515 M       Avar
                 DAG-529 M Avar_outlier2
                 Avar.SG M    Avar.SG
                     AV1 F Hungary_Avar_5
                     AV2 F Hungary_Avar_5_daughter.or.mother.AV1
                  SZ1.SG M Hungary_Avar_1.SG



```bash
%%bash
grep Avar v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
grep -v Hungary|sed -r 's/^\s+//'|cut -d ' ' -f 1 > Avar.id && wc -l Avar.id
```

    8 Avar.id



```bash
%%bash
grep Adygei v44.3_HO_public.ind
```

                 Adg-185 F     Adygei
                 Adg-192 M     Adygei
                 Adg-194 M     Adygei
                 Adg-222 M     Adygei
                 Adg-224 M     Adygei
                 Adg-226 F     Adygei
                 Adg-227 F     Adygei
                 Adg-238 M     Adygei
                lpgh-303 M     Adygei
                 Shap-55 F     Adygei
                Shap-156 M     Adygei
                Shap-225 M     Adygei
                    YY-2 M     Adygei
                   YY-18 M     Adygei
                   YY-36 M     Adygei
               HGDP01381 F     Adygei
               HGDP01382 F     Adygei
               HGDP01383 M     Adygei
               HGDP01384 F Ignore_Adygei
               HGDP01385 M     Adygei
               HGDP01386 F     Adygei
               HGDP01387 F     Adygei
               HGDP01388 F Ignore_Adygei
               HGDP01396 M     Adygei
               HGDP01397 M     Adygei
               HGDP01398 F     Adygei
               HGDP01399 F     Adygei
               HGDP01400 F     Adygei
               HGDP01401 F     Adygei
               HGDP01402 M     Adygei
               HGDP01403 M     Adygei
               HGDP01404 M     Adygei
                 NA13626 M     Adygei
                 NA13619 M Ignore_Adygei(relative.of.HGDP01382)
                 NA13622 M Ignore_Adygei(duplicate)
                 NA13624 F Ignore_Adygei(duplicate)
                 NA13625 F Ignore_Adygei(duplicate)
                 NA13617 F Ignore_Adygei(duplicate)
                 NA13618 F Ignore_Adygei(duplicate)
                 NA13620 M Ignore_Adygei(duplicate)
           S_Adygei-1.DG M  Adygei.DG
           S_Adygei-2.DG F  Adygei.DG
           HGDP01385.SDG M Adygei.SDG
           HGDP01396.SDG M Ignore_Adygei.SDG
           HGDP01383.SDG M Adygei.SDG
           HGDP01404.SDG M Adygei.SDG
           HGDP01384.SDG F Ignore_Adygei.SDG
           HGDP01388.SDG F Ignore_Adygei.SDG
           HGDP01400.SDG F Adygei.SDG
           HGDP01398.SDG F Adygei.SDG
           HGDP01382.SDG F Adygei.SDG
           HGDP01397.SDG M Adygei.SDG
           HGDP01386.SDG F Adygei.SDG
           HGDP01403.SDG M Adygei.SDG
           HGDP01399.SDG F Adygei.SDG
           HGDP01387.SDG F Adygei.SDG
           HGDP01402.SDG M Adygei.SDG
           HGDP01401.SDG F Adygei.SDG



```bash
%%bash
grep Adygei v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Adygei.id && wc -l Adygei.id
```

    31 Adygei.id



```bash
%%bash
grep Ingushian v44.3_HO_public.ind
```

                ING-1154 M  Ingushian
                ING-1167 M  Ingushian
                ING-1168 M  Ingushian
                ING-1169 M  Ingushian
                ING-1176 M  Ingushian
                ING-1177 M  Ingushian
                ING-1178 M  Ingushian
                ING-1179 M  Ingushian
                ING-1206 M  Ingushian
                ING-1330 M  Ingushian



```bash
%%bash
grep Ingushian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Ingushian.id && wc -l Ingushian.id
```

    10 Ingushian.id



```bash
%%bash
grep Darginian v44.3_HO_public.ind
```

                 DAG-146 M  Darginian
                 DAG-161 M  Darginian
                 DAG-165 M  Darginian
                 DAG-206 M  Darginian
                 DAG-225 M  Darginian
                 DAG-231 M  Darginian
                 DAG-232 M  Darginian
                 DAG-350 M  Darginian



```bash
%%bash
grep Darginian v44.3_HO_public.ind|grep -v [Jj]ew|grep -v [Ii]gnore|grep -v SG$|grep -v SDG$|grep -v DG$|grep -v outlier|\
sed -r 's/^\s+//'|cut -d ' ' -f 1 > Darginian.id && wc -l Darginian.id
```

    8 Darginian.id


## Therefore, concatenate all <code>.id</code> files


```bash
%%bash
ls *.id; ls *.id|wc -l
```

    Adygei.id	      French_south.id	     Polish.id
    Albanian.id	      Georgian.id	     Romanian.id
    Algerian.id	      Greek.id		     Russian.id
    Armenian_Hemsheni.id  Hungarian.id	     Russian_nordic.id
    Armenian.id	      Icelandic.id	     Sardinian.id
    Assyrian.id	      Ingushian.id	     Saudi.id
    Avar.id		      Iranian_gulf.id	     Scottish.id
    Azeri.id	      Iranian.id	     Spanish.id
    Bashkir.id	      Italian_north.id	     Spanish_north.id
    Basque.id	      Italian_south.id	     Syrian.id
    Bedouin.id	      Jordanian.id	     Tatar_Kazan.id
    Belarusian.id	      Kurd.id		     Tatar_Siberian.id
    Bulgarian.id	      Lebanese_Christian.id  Tatar_Volga.id
    Chechen.id	      Lebanese_Muslim.id     Tunisian.id
    Chuvash.id	      Lebanese_unknown.id    Turkish_central.id
    Croatian.id	      Libyan.id		     Turkish_mixed.id
    Czech.id	      Lithuanian.id	     Turkish_northeast.id
    Darginian.id	      Maltese.id	     Turkish_south.id
    Druze.id	      Moldavian.id	     Turkish_west.id
    Egyptian.id	      Mordovian.id	     Turkmen.id
    English.id	      Moroccan.id	     Tuscan.id
    Estonian.id	      Norwegian.id	     Ukrainian.id
    Finnish.id	      Palestinian.id	     Ukrainian_north.id
    French.id	      Pathan.id		     Uzbek.id
    72


### 72 population groups are present


```bash
%%bash
cat *.id > indivs_for_ibd.id
```

    cat: indivs_for_ibd.id: input file is output file



```bash
%%bash
cat indivs_for_ibd.id|sort|uniq|wc -l
```

    1305



```bash
%%bash
cat v44.3_HO_public.fam|sed -r 's/^\s+//'|sed -r 's/\s+/ /'|\
cut -d $' ' -f 1,2|tr ' ' '\n' > v44.3_HO_public.fam.keep
```


```bash
%%bash
head -n5 v44.3_HO_public.fam.keep
```

    1
    MAL-005
    2
    MAL-009
    3



```bash
%%bash
grep -B 1 -Fxf indivs_for_ibd.id v44.3_HO_public.fam.keep|sed -E '/^--$/d'|paste - -|head -n 5
```

    47	Adg-185
    48	Adg-192
    49	Adg-194
    50	Adg-222
    51	Adg-224



```bash
%%bash
grep -B 1 -Fxf indivs_for_ibd.id v44.3_HO_public.fam.keep|sed -E '/^--$/d'|paste - - > indivs_for_ibd.keep
```


```bash
%%bash
wc -l indivs_for_ibd.keep
```

    1305 indivs_for_ibd.keep



```bash
%%bash
../plink --bfile v44.3_HO_public --keep indivs_for_ibd.keep --make-bed --out ibd_non_Jews_main
```

    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to ibd_non_Jews_main.log.
    Options in effect:
      --bfile v44.3_HO_public
      --keep indivs_for_ibd.keep
      --make-bed
      --out ibd_non_Jews_main
    
    15952 MB RAM detected; reserving 7976 MB for main workspace.
    597573 variants loaded from .bim file.
    13197 people (8042 males, 4999 females, 156 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to ibd_non_Jews_main.nosex .
    13197 phenotype values loaded from .fam.
    --keep: 1305 people remaining.
    Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
    phenotypes to be ignored, use the --allow-no-sex flag.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 1305 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    treat these as missing.
    Total genotyping rate in remaining samples is 0.985945.
    597573 variants and 1305 people pass filters and QC.
    Among remaining phenotypes, 0 are cases and 1305 are controls.
    --make-bed to ibd_non_Jews_main.bed + ibd_non_Jews_main.bim +
    ibd_non_Jews_main.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.

### Extract German, Latvian, Polish, Swedish and Romani individuals from the high-coerage SNP dataset Pagani et. al., 2016

#### Step 1. Merge the original dataset from 22 seperate sets of files chromosome-wise

```bash
%%bash
curl -sSL https://evolbio.ut.ee/CGgenomes_VCF/EGDP_PaganiEtAl2016_Release_Build37_PLINK.tar.gz -o PaganiEtAl2016Dataset.tar.gz &&
tar -xvzf PaganiEtAl2016Dataset.tar.gz
```

```bash
%%bash
ls PaganiEtAl2016Dataset/EGDP_PaganiEtAl2016_Release_Build37_PLINK/
```

    EGDP_PaganiEtAl2016_Release_Build37_Chr10.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr10.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr10.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr11.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr11.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr11.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr12.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr12.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr12.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr13.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr13.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr13.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr14.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr14.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr14.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr15.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr15.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr15.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr16.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr16.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr16.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr17.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr17.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr17.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr18.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr18.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr18.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr19.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr19.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr19.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr1.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr1.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr1.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr20.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr20.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr20.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr21.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr21.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr21.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr22.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr22.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr22.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr2.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr2.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr2.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr3.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr3.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr3.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr4.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr4.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr4.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr5.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr5.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr5.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr6.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr6.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr6.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr7.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr7.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr7.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr8.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr8.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr8.vcf_plink.fam
    EGDP_PaganiEtAl2016_Release_Build37_Chr9.vcf_plink.bed
    EGDP_PaganiEtAl2016_Release_Build37_Chr9.vcf_plink.bim
    EGDP_PaganiEtAl2016_Release_Build37_Chr9.vcf_plink.fam
    Samples.xlsx



```bash
%%bash
ssconvert PaganiEtAl2016Dataset/EGDP_PaganiEtAl2016_Release_Build37_PLINK/Samples.xlsx PaganiEtAl2016Samples.csv
```

### Extract only the German, Latvian, Polish, Swedish and Romani individuals


```bash
%%bash
grep -E "(,Roma,|,Swedes,|,Germans,|,Latvians,|,Poles,)" PaganiEtAl2016Samples.csv
grep -E "(,Roma,|,Swedes,|,Germans,|,Latvians,|,Poles,)" PaganiEtAl2016Samples.csv|wc -l
```

    Ger1,GS000016892-ASM,1,1,"S&W Europe",0,female,,Germans,Western-Europe,Germany,Europe,"Western Europe",Indo-European/Germanic,approximate_probably_born_in_Germany,"this study","52,5","10,1",blood,normal,,,,,,135,2.2.0.26,"2,2",German13,German13,,GS01565-DNA_D06,GS000021467-DID,GS000016892-ASM,TRUE,"660 000","Human660W-Quad v1.0"
    Ger2,GS000016893-ASM,1,1,"S&W Europe",0,female,,Germans,Western-Europe,Germany,Europe,"Western Europe",Indo-European/Germanic,approximate_germany_actual_birthplace_Järvamaa,"this study","52,5","10,1",blood,normal,,,,,,135,2.2.0.26,"2,2",German16,German16,,GS01565-DNA_E06,GS000021468-DID,GS000016893-ASM,TRUE,"660 000","Human660W-Quad v1.0"
    Ger3,GS000016891-ASM,1,1,"S&W Europe",1,male,R1b3,Germans,Western-Europe,Germany,Europe,"Western Europe",Indo-European/Germanic,approximate_germany_actual_birthplace_Sverdlovskoe_oblast,"this study","52,5","10,1",blood,normal,,,,,,135,2.2.0.26,"2,2",German6,german_V01118,German6,GS01565-DNA_C06,GS000021466-DID,GS000016891-ASM,TRUE,"660 000","Human660W-Quad v1.0"
    Lat1,GS000016903-ASM,1,1,"E&N Europe",0,female,,Latvians,,Latvia,Europe,"Northern Europe",Indo-European/Balto-Slavic/Baltic,approximate,"this study",57,"24,4",blood,normal,,,,,,135,2.2.0.26,"2,2",latvian_V45318,latvian_V45318,,GS01565-DNA_A08,GS000021480-DID,GS000016903-ASM,TRUE,"1 000 000",HumanOmni1-Quad
    Lat2,GS000035027-ASM,1,1,"E&N Europe",0,male,N3a3a,Latvians,Cesis,Latvia,Europe,"Northern Europe",Indo-European/Balto-Slavic/Baltic,"parents from Cesis and Jelgava_Liepaja. Coordinates are approximate of those of parents","this study","56,9","24,5",blood,normal,,,,,,137,2.4.0.43,"2,4",NOT_DONE,NOT_DONE,,GS01548-DNA_G02,GS000035456-DID,GS000035027-ASM,TRUE,,
    Lat3,GS000035148-ASM,1,1,"E&N Europe",0,male,N3a3a,Latvians,Dobele,Latvia,Europe,"Northern Europe",Indo-European/Balto-Slavic/Baltic,"parents from Kurzeme and Zemgale regions. Cooordinates are of birthplace which is inbetween of parents origins","this study","56,6","23,3",blood,normal,,,,,,137,2.4.0.43,"2,4",NOT_DONE,NOT_DONE,,GS01548-DNA_H02,GS000035457-DID,GS000035148-ASM,TRUE,,
    Pole1,GS000016888-ASM,1,1,"E&N Europe",0,female,,Poles,,Poland,Europe,"Eastern Europe",Indo-European/Balto-Slavic/Slavic/West,approximate_but_born_in_Kohtla_Järve,"this study","52,5","20,8",blood,normal,,,,,,135,2.2.0.26,"2,2",Poland137,Poland137,polish_V2873,GS01565-DNA_H05,GS000021463-DID,GS000016888-ASM,TRUE,"660 000","Human660W-Quad v1.0"
    Pole2,GS000016889-ASM,1,1,"E&N Europe",0,female,,Poles,,Poland,Europe,"Eastern Europe",Indo-European/Balto-Slavic/Slavic/West,approximate,"this study","52,5","20,8",blood,normal,,,,,,135,2.2.0.26,"2,2",Poland141,Poland141,,GS01565-DNA_A06,GS000021464-DID,GS000016889-ASM,TRUE,"660 000","Human660W-Quad v1.0"
    Pole3,GS000016890-ASM,1,1,"E&N Europe",0,female,,Poles,,Poland,Europe,"Eastern Europe",Indo-European/Balto-Slavic/Slavic/West,approximate,"this study","52,5","20,8",blood,normal,,,,,,135,2.2.0.26,"2,2",Poland150,Poland150,,GS01565-DNA_B06,GS000021465-DID,GS000016890-ASM,TRUE,"660 000","Human660W-Quad v1.0"
    Pole4,GS000015869-ASM,1,1,"E&N Europe",0,female,,Poles,Wroclaw,Poland,Europe,"Eastern Europe",Indo-European/Balto-Slavic/Slavic/West,,"this study","51,1",17,blood,normal,,,,,,135,2.2.0.26,"2,2",pole12p_1m,pole12p_1m,,GS01562-DNA_B01,GS000020360-DID,GS000015869-ASM,TRUE,"1 000 000",HumanOmni1-Quad
    RomBH1,GS000015870-ASM,1,0,,1,female,,Roma,Zavidovici,Bosnia-Herzegovina,Europe,"Southern Europe",Indo-European/Indo-Iranian,"3rd degree relative with RomaBH5","this study","44,4","18,1",blood,normal,,,,,,135,2.2.0.26,"2,2",romaBH1_1m,romaBH1_1m,,GS01562-DNA_C01,GS000020361-DID,GS000015870-ASM,TRUE,"1 000 000",HumanOmni1-Quad
    RomBH2,GS000014325-ASM,1,0,,0,female,,Roma,Zavidovici,Bosnia-Herzegovina,Europe,"Southern Europe",Indo-European/Indo-Iranian,,"this study","44,4","18,1",blood,normal,,,,,,135,2.2.0.26,"2,2",romaBH2_1m,romaBH2_1m,,GS01567-DNA_D01,GS000018360-DID,GS000014325-ASM,TRUE,"1 000 000",HumanOmni1-Quad
    RomBH5,GS000014352-ASM,1,0,,1,female,,Roma,Zavidovici,Bosnia-Herzegovina,Europe,"Southern Europe",Indo-European/Indo-Iranian,"3rd degree relative with RomaBH1","this study","44,4","18,1",blood,normal,,,,,,135,2.2.0.26,"2,2",romaBH5_1m,romaBH5_1m,,GS01567-DNA_E01,GS000018361-DID,GS000014352-ASM,TRUE,"1 000 000",HumanOmni1-Quad
    Swe1,GS000035109-ASM,1,1,"E&N Europe",0,male,R1a1e,Swedes,,Sweden,Europe,"Northern Europe",Indo-European/Germanic,"approximate location (capital)","this study","59,4",18,blood,normal,,,,,,137,2.4.0.43,"2,4",swede_V49245,swede_V49245,,GS01548-DNA_A05,GS000035474-DID,GS000035109-ASM,TRUE,"1 000 000",HumanOmni1-Quad
    Swe2,GS000035240-ASM,1,1,"E&N Europe",0,male,R1a1e,Swedes,Nyköping,Sweden,Europe,"Northern Europe",Indo-European/Germanic,"parents and grandparents from same region","this study","58,8",17,blood,normal,,,,,,137,2.4.0.43,"2,4",NOT_DONE,NOT_DONE,,GS01548-DNA_E08,GS000035502-DID,GS000035240-ASM,TRUE,,
    15


### Therefore, we would have 15 individuals to extract


```bash
%%bash
grep ",Poles," PaganiEtAl2016Samples.csv|cut -d , -f 2 >> Polish.id && wc -l Polish.id
```

    5 Polish.id



```bash
%%bash
grep ",Roma," PaganiEtAl2016Samples.csv|cut -d , -f 2 > Romani.id && wc -l Romani.id
grep ",Germans," PaganiEtAl2016Samples.csv|cut -d , -f 2 > German.id && wc -l German.id
grep ",Swedes," PaganiEtAl2016Samples.csv|cut -d , -f 2 > Swede.id && wc -l Swede.id
grep ",Latvians," PaganiEtAl2016Samples.csv|cut -d , -f 2 > Latvian.id && wc -l Latvian.id
```

    3 Romani.id
    3 German.id
    2 Swede.id
    3 Latvian.id



```bash
%%bash
ls *.id|grep -E -v "(^HGDP|^indiv)"|wc -l  ## we would have 76 non-Jewish populations to study
```

    76


```bash
%%bash
grep -E "(,Roma,|,Swedes,|,Germans,|,Latvians,|,Poles,)" PaganiEtAl2016Samples.csv|cut -d , -f 2|\
awk -v OFS=" " '{print $0,$0 >> "PaganiEtAl2016_15indivs.keep"}' && cat PaganiEtAl2016_15indivs.keep &&
wc -l PaganiEtAl2016_15indivs.keep
```

    GS000016892-ASM GS000016892-ASM
    GS000016893-ASM GS000016893-ASM
    GS000016891-ASM GS000016891-ASM
    GS000016903-ASM GS000016903-ASM
    GS000035027-ASM GS000035027-ASM
    GS000035148-ASM GS000035148-ASM
    GS000016888-ASM GS000016888-ASM
    GS000016889-ASM GS000016889-ASM
    GS000016890-ASM GS000016890-ASM
    GS000015869-ASM GS000015869-ASM
    GS000015870-ASM GS000015870-ASM
    GS000014325-ASM GS000014325-ASM
    GS000014352-ASM GS000014352-ASM
    GS000035109-ASM GS000035109-ASM
    GS000035240-ASM GS000035240-ASM
    15 PaganiEtAl2016_15indivs.keep


### Merge 22 .bim files at the original dataset of Pagani et. al., 2016 and extract only the overlapping SNps with the 1240K+HO dataset


```bash
%%bash
ls PaganiEtAl2016Dataset/EGDP_PaganiEtAl2016_Release_Build37_PLINK/*.bim|wc -l
cat PaganiEtAl2016Dataset/EGDP_PaganiEtAl2016_Release_Build37_PLINK/*.bim > PaganiEtAl2016_allSNPs.bim &&
wc -l PaganiEtAl2016_allSNPs.bim
```

    22
    42971058 PaganiEtAl2016_allSNPs.bim

```bash
%%bash
head -n5 PaganiEtAl2016_allSNPs.bim &&
cat PaganiEtAl2016_allSNPs.bim|cut -d $'\t' -f 2 > PaganiEtAl2016_allSNPs.txt
head -n3 PaganiEtAl2016_allSNPs.txt && wc -l PaganiEtAl2016_allSNPs.txt
wc -l ibd_non_Jews_main.bim &&
head -n5 ibd_non_Jews_main.bim
cat ibd_non_Jews_main.bim|cut -d $'\t' -f 1,4|tr $'\t' : > 1240K_HO_snps_alt_form.txt &&
head -n3 1240K_HO_snps_alt_form.txt && wc -l 1240K_HO_snps_alt_form.txt
```

    10	10:62010	0	62010	T	C
    10	10:62493	0	62493	T	C
    10	10:68303	0	68303	G	A
    10	10:69967	0	69967	T	C
    10	10:70769	0	70769	A	C
    10:62010
    10:62493
    10:68303
    42971058 PaganiEtAl2016_allSNPs.txt
    597573 ibd_non_Jews_main.bim
    1	rs3094315	0.02013	752566	G	A
    1	rs7419119	0.022518	842013	G	T
    1	rs13302957	0.024116	891021	G	A
    1	rs6696609	0.024457	903426	T	C
    1	rs8997	0.025727	949654	A	G
    1:752566
    1:842013
    1:891021
    597573 1240K_HO_snps_alt_form.txt



```bash
%%bash
grep -Fxf 1240K_HO_snps_alt_form.txt PaganiEtAl2016_allSNPs.txt|wc -l
grep -Fxf 1240K_HO_snps_alt_form.txt PaganiEtAl2016_allSNPs.txt > PaganiEtAl2016_1240K_HO_snps.txt
```

    539759



```bash
%%bash
for chr in {2..22}; do
echo -e "PaganiEtAl2016Dataset/EGDP_PaganiEtAl2016_Release_Build37_PLINK/\
EGDP_PaganiEtAl2016_Release_Build37_Chr${chr}.vcf_plink.bed \
PaganiEtAl2016Dataset/EGDP_PaganiEtAl2016_Release_Build37_PLINK/\
EGDP_PaganiEtAl2016_Release_Build37_Chr${chr}.vcf_plink.bim \
PaganiEtAl2016Dataset/EGDP_PaganiEtAl2016_Release_Build37_PLINK/\
EGDP_PaganiEtAl2016_Release_Build37_Chr${chr}.vcf_plink.fam" >> PaganiEtAl2016_mergelist.txt
done && wc -l PaganiEtAl2016_mergelist.txt
```

    21 PaganiEtAl2016_mergelist.txt



```bash
%%bash
../plink --bfile PaganiEtAl2016Dataset/EGDP_PaganiEtAl2016_Release_Build37_PLINK/\
EGDP_PaganiEtAl2016_Release_Build37_Chr1.vcf_plink --merge-list PaganiEtAl2016_mergelist.txt \
--extract PaganiEtAl2016_1240K_HO_snps.txt --keep PaganiEtAl2016_15indivs.keep --allow-no-sex --autosome \
--make-bed --out PaganiEtAl2016_15indivs_1240K_HO
```

    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to PaganiEtAl2016_15indivs_1240K_HO.log.
    Options in effect:
      --allow-no-sex
      --autosome
      --bfile PaganiEtAl2016Dataset/EGDP_PaganiEtAl2016_Release_Build37_PLINK/EGDP_PaganiEtAl2016_Release_Build37_Chr1.vcf_plink
      --extract PaganiEtAl2016_1240K_HO_snps.txt
      --keep PaganiEtAl2016_15indivs.keep
      --make-bed
      --merge-list PaganiEtAl2016_mergelist.txt
      --out PaganiEtAl2016_15indivs_1240K_HO
    
    15952 MB RAM detected; reserving 7976 MB for main workspace.
    Performing 2-pass merge (402 people, 38236289/42971058 variants per pass).
    Pass 1 complete.
    Merged fileset written to PaganiEtAl2016_15indivs_1240K_HO-merge.bed +
    PaganiEtAl2016_15indivs_1240K_HO-merge.bim +
    PaganiEtAl2016_15indivs_1240K_HO-merge.fam .
    42971058 variants loaded from .bim file.
    402 people (0 males, 0 females, 402 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to PaganiEtAl2016_15indivs_1240K_HO.nosex .
    --extract: 539759 variants remaining.
    --keep: 15 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 15 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate in remaining samples is exactly 1.
    539759 variants and 15 people pass filters and QC.
    Note: No phenotypes present.
    --make-bed to PaganiEtAl2016_15indivs_1240K_HO.bed +
    PaganiEtAl2016_15indivs_1240K_HO.bim + PaganiEtAl2016_15indivs_1240K_HO.fam ...
    0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.


## Modify <code>.fam</code> and <code>.bim</code> files of <code>PaganiEtAl2016_15indivs_1240K_HO</code>. Change all family names to 1 and modify SNP ids to #rs


```bash
%%bash
mv PaganiEtAl2016_15indivs_1240K_HO.bim PaganiEtAl2016_15indivs_1240K_HO.bim.bk
mv PaganiEtAl2016_15indivs_1240K_HO.fam PaganiEtAl2016_15indivs_1240K_HO.fam.bk
```


```python
import pandas as pd

## alter .bim file
bim_header = ['chr', 'snp_id', 'cM', 'bp_pos', 'ref', 'alt']
HO_1240K_bim = pd.read_csv("ibd_non_Jews_main.bim", sep="\t", header=None, names=bim_header)
HO_1240K_bim['alt_snp_id'] = HO_1240K_bim['chr'].astype(str) + ':' + HO_1240K_bim['bp_pos'].astype(str)
pagani2016_HO_bim = pd.read_csv("PaganiEtAl2016_15indivs_1240K_HO.bim.bk", sep="\t", header=None, names=bim_header)
pagani2016_HO_snps = pagani2016_HO_bim[['chr', 'snp_id']]
pagani2016_HO_snps.rename(columns={'snp_id': 'alt_snp_id'}, inplace=True)
pagani2016_HO_new_bim = pagani2016_HO_snps.merge(HO_1240K_bim, how='inner', on=['chr', 'alt_snp_id'])
pagani2016_HO_new_bim = pagani2016_HO_new_bim.drop(['alt_snp_id'], axis=1)
print(pagani2016_HO_new_bim)
pagani2016_HO_new_bim.to_csv("PaganiEtAl2016_15indivs_1240K_HO.bim", sep='\t', header=False, index=False)
## alter .fam file
fam_header = ['fid', 'iid', 'pid', 'mid', 'sex', 'pheno']
pagani2016_HO_fam = pd.read_csv("PaganiEtAl2016_15indivs_1240K_HO.fam.bk", sep=" ", header=None, names=fam_header)
pagani2016_HO_fam['fid'] = "1"
pagani2016_HO_fam.to_csv("PaganiEtAl2016_15indivs_1240K_HO.fam", sep='\t', header=False, index=False)
```

                chr       snp_id        cM    bp_pos ref alt
    0         1    rs7419119  0.022518    842013   G   T
    1         1   rs13302957  0.024116    891021   G   A
    2         1    rs6696609  0.024457    903426   T   C
    3         1    rs9442372  0.026288   1018704   A   G
    4         1  rs147606383  0.026665   1045331   A   G
    ...     ...          ...       ...       ...  ..  ..
    539754   22     rs715586  0.740465  51163138   T   C
    539755   22    rs8137951  0.740478  51165664   A   G
    539756   22    rs3810648  0.740745  51175626   G   A
    539757   22   rs78827609  0.740750  51176347   T   C
    539758   22  rs116026001  0.740758  51177526   G   A
    
    [539759 rows x 6 columns]



### Merge the HRC Europe panel-imputed imputed Jewish genotype dataset with the ones of non-Jews <br />
### The merged dataset will be input for <code>phasedibd</code> analysis

### List all SNPs in the imputed Jewish dataset


```bash
%%bash
../plink --bfile ../Jews_Genetics_Project/HRC_eu_ref_imputed_935Jews_for_IBD/merged_HRC_Eu_panel_imputed_1240K_935Jews \
--write-snplist --out imputed_935Jews
```

## 935 Jews with 1320 selected non-Jews for groupwise IBD analysis)

### Steps for merging datasets:
#### 1. using the default option to create a union, found out same-position conflicts
#### 2. use <code>--list-duplicate-vars</code> option in PLINK to list the confliting SNPs at Step 1.
#### 3. use <code>.dupvar</code> file generated at the 2nd step and <code>../../py_scripts/correct_bim_file.py</code> to correct <code>.bim</code> files of both datasets using the Jewish one as benchmark.
#### 4. Merge the two datasets with corrected <code>.bim</code> files. Sort the overlapping SNPs out using <code>--write-snplist</code> and <code>--extract</code> option in PLINK to obtain the final dataset. 


```bash
%%bash
head -n5 ibd_non_Jews_main.bim
```

    1	rs3094315	0.02013	752566	G	A
    1	rs7419119	0.022518	842013	G	T
    1	rs13302957	0.024116	891021	G	A
    1	rs6696609	0.024457	903426	T	C
    1	rs8997	0.025727	949654	A	G


## Non-Jews with imputed Jews

### There's Morgan annotation in <code>ibd_non_Jews_main.bim</code>. Therefore, take it as reference instead of the Jewish dataset.


```bash
%%bash
../plink --bfile ibd_non_Jews_main --bmerge ../merged_HRC_Eu_panel_imputed_1240K_935Jews --autosome \
--allow-no-sex --make-bed --out merged_imputed935Jews_1305nonJews
```

    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to merged_imputed935Jews_1305nonJews.log.
    Options in effect:
      --allow-no-sex
      --autosome
      --bfile ibd_non_Jews_main
      --bmerge ../merged_HRC_Eu_panel_imputed_1240K_935Jews
      --make-bed
      --out merged_imputed935Jews_1305nonJews
    
    64039 MB RAM detected; reserving 32019 MB for main workspace.
    1305 people loaded from ibd_non_Jews_main.fam.
    935 people to be merged from ../merged_HRC_Eu_panel_imputed_1240K_935Jews.fam.
    Of these, 935 are new, while 0 are present in the base dataset.
    597573 markers loaded from ibd_non_Jews_main.bim.
    1067391 markers to be merged from
    ../merged_HRC_Eu_panel_imputed_1240K_935Jews.bim.
    Of these, 555820 are new, while 511571 are present in the base dataset.
    Performing single-pass merge (2240 people, 1148944 variants).
    Merged fileset written to merged_imputed935Jews_1305nonJews-merge.bed +
    merged_imputed935Jews_1305nonJews-merge.bim +
    merged_imputed935Jews_1305nonJews-merge.fam .
    1148944 variants loaded from .bim file.
    2240 people (999 males, 301 females, 940 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to merged_imputed935Jews_1305nonJews.nosex .
    1305 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 2240 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.660108.
    1148944 variants and 2240 people pass filters and QC.
    Among remaining phenotypes, 0 are cases and 1305 are controls.  (935 phenotypes
    are missing.)
    --make-bed to merged_imputed935Jews_1305nonJews.bed +
    merged_imputed935Jews_1305nonJews.bim + merged_imputed935Jews_1305nonJews.fam
    ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.


    Warning: Variants 'rs146976167' and '1:27112603:A:T' have the same position.
    Warning: Variants 'Affx-8828625' and '1:30763015:A:T' have the same position.
    Warning: Variants 'rs192299508' and '1:32410196:C:G' have the same position.
    154 more same-position warnings: see log file.



```bash
%%bash
../plink  --bfile merged_imputed935Jews_1305nonJews --allow-no-sex \
--list-duplicate-vars --out merged_imputed935Jews_1305nonJews
```

    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to merged_imputed935Jews_1305nonJews.log.
    Options in effect:
      --allow-no-sex
      --bfile merged_imputed935Jews_1305nonJews
      --list-duplicate-vars
      --out merged_imputed935Jews_1305nonJews
    
    64039 MB RAM detected; reserving 32019 MB for main workspace.
    1148944 variants loaded from .bim file.
    2240 people (999 males, 301 females, 940 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to merged_imputed935Jews_1305nonJews.nosex .
    1305 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 2240 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.660108.
    1148944 variants and 2240 people pass filters and QC.
    Among remaining phenotypes, 0 are cases and 1305 are controls.  (935 phenotypes
    are missing.)
    --list-duplicate-vars report written to
    merged_imputed935Jews_1305nonJews.dupvar .



```bash
%%bash
python ../Jews_Genetics_Project/py_scripts/correct_bim_file.py \
--ref-map ../Jews_Genetics_Project/HRC_eu_ref_imputed_935Jews_for_IBD/merged_HRC_Eu_panel_imputed_1240K_935Jews.bim \
--dupvar merged_imputed935Jews_1305nonJews.dupvar \
--bim1 ../Jews_Genetics_Project/HRC_eu_ref_imputed_935Jews_for_IBD/merged_HRC_Eu_panel_imputed_1240K_935Jews.bim \
--bim2 ibd_non_Jews_main.bim
```

### Extract the overlapping SNPs between "1240K+HO" and the imputed Jewish dataset


```bash
%%bash 
../plink --bfile ibd_non_Jews_main --extract imputed_935Jews.snplist --allow-no-sex --make-bed \
--out ibd_non_Jews_main_overlap
```

    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to ibd_non_Jews_main_overlap.log.
    Options in effect:
      --allow-no-sex
      --bfile ibd_non_Jews_main
      --extract imputed_935Jews.snplist
      --make-bed
      --out ibd_non_Jews_main_overlap
    
    64039 MB RAM detected; reserving 32019 MB for main workspace.
    597573 variants loaded from .bim file.
    1305 people (999 males, 301 females, 5 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to ibd_non_Jews_main_overlap.nosex .
    1305 phenotype values loaded from .fam.
    --extract: 511571 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 1305 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.990643.
    511571 variants and 1305 people pass filters and QC.
    Among remaining phenotypes, 0 are cases and 1305 are controls.
    --make-bed to ibd_non_Jews_main_overlap.bed + ibd_non_Jews_main_overlap.bim +
    ibd_non_Jews_main_overlap.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.



```bash
%%bash
../plink --bfile ibd_non_Jews_main_overlap --allow-no-sex --write-snplist --out ibd_non_Jews_main_overlap
```

    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to ibd_non_Jews_main_overlap.log.
    Options in effect:
      --allow-no-sex
      --bfile ibd_non_Jews_main_overlap
      --out ibd_non_Jews_main_overlap
      --write-snplist
    
    64039 MB RAM detected; reserving 32019 MB for main workspace.
    511571 variants loaded from .bim file.
    1305 people (999 males, 301 females, 5 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to ibd_non_Jews_main_overlap.nosex .
    1305 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 1305 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.990643.
    511571 variants and 1305 people pass filters and QC.
    Among remaining phenotypes, 0 are cases and 1305 are controls.
    List of variant IDs written to ibd_non_Jews_main_overlap.snplist .


```bash
%%bash
../plink --bfile ibd_non_Jews_main \
--bmerge ../Jews_Genetics_Project/HRC_eu_ref_imputed_935Jews_for_IBD/merged_HRC_Eu_panel_imputed_1240K_935Jews --autosome \
--allow-no-sex --extract ibd_non_Jews_main_overlap.snplist --make-bed \
--out merged_imputed935Jews_1305nonJews_overlap
```


### Check if all SNPs are annotated by genetic distances (Morgan)


```bash
%%bash
rm merged_imputed935Jews_1305nonJews.*
```


```bash
%%bash
awk 'BEGIN{FS="\t"; count=0} {if($3!='0'){count++;}} END{print count}' \
merged_imputed935Jews_1305nonJews_overlap.bim
```

    511571


### Add the extra <code>PaganiEtAl2016_15indivs_1240K_HO</code> dataset then

### Add the extra <code>PaganiEtAl2016_15indivs_1240K_HO</code> dataset then


```bash
%%bash
../plink --bfile merged_imputed935Jews_1305nonJews_overlap --bmerge PaganiEtAl2016_15indivs_1240K_HO \
--allow-no-sex --make-bed --out merged_imputed935Jews_1320nonJews_overlap
```

    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to merged_imputed935Jews_1320nonJews_overlap.log.
    Options in effect:
      --allow-no-sex
      --bfile merged_imputed935Jews_1305nonJews_overlap
      --bmerge PaganiEtAl2016_15indivs_1240K_HO
      --make-bed
      --out merged_imputed935Jews_1320nonJews_overlap
    
    64039 MB RAM detected; reserving 32019 MB for main workspace.
    2240 people loaded from merged_imputed935Jews_1305nonJews_overlap.fam.
    15 people to be merged from PaganiEtAl2016_15indivs_1240K_HO.fam.
    Of these, 15 are new, while 0 are present in the base dataset.
    511571 markers loaded from merged_imputed935Jews_1305nonJews_overlap.bim.
    539759 markers to be merged from PaganiEtAl2016_15indivs_1240K_HO.bim.
    Of these, 40713 are new, while 499046 are present in the base dataset.
    Performing single-pass merge (2255 people, 552284 variants).
    Merged fileset written to merged_imputed935Jews_1320nonJews_overlap-merge.bed +
    merged_imputed935Jews_1320nonJews_overlap-merge.bim +
    merged_imputed935Jews_1320nonJews_overlap-merge.fam .
    552284 variants loaded from .bim file.
    2255 people (999 males, 301 females, 955 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to merged_imputed935Jews_1320nonJews_overlap.nosex .
    1305 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 2255 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.89311.
    552284 variants and 2255 people pass filters and QC.
    Among remaining phenotypes, 0 are cases and 1305 are controls.  (950 phenotypes
    are missing.)
    --make-bed to merged_imputed935Jews_1320nonJews_overlap.bed +
    merged_imputed935Jews_1320nonJews_overlap.bim +
    merged_imputed935Jews_1320nonJews_overlap.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.


### Apply QC to the merged dataset. min genotyping rate=0.85, min MAF=0.001


```bash
%%bash
../plink --bfile merged_imputed935Jews_1320nonJews_overlap --geno 0.15 --maf 0.001 --allow-no-sex --make-bed \
--out merged_imputed935Jews_1320nonJews_overlap_QC
```

    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to merged_imputed935Jews_1320nonJews_overlap_QC.log.
    Options in effect:
      --allow-no-sex
      --bfile merged_imputed935Jews_1320nonJews_overlap
      --geno 0.15
      --maf 0.001
      --make-bed
      --out merged_imputed935Jews_1320nonJews_overlap_QC
    
    64039 MB RAM detected; reserving 32019 MB for main workspace.
    552284 variants loaded from .bim file.
    2255 people (999 males, 301 females, 955 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to merged_imputed935Jews_1320nonJews_overlap_QC.nosex
    .
    1305 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 2255 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.89311.
    96644 variants removed due to missing genotype data (--geno).
    94 variants removed due to minor allele threshold(s)
    (--maf/--max-maf/--mac/--max-mac).
    455546 variants and 2255 people pass filters and QC.
    Among remaining phenotypes, 0 are cases and 1305 are controls.  (950 phenotypes
    are missing.)
    --make-bed to merged_imputed935Jews_1320nonJews_overlap_QC.bed +
    merged_imputed935Jews_1320nonJews_overlap_QC.bim +
    merged_imputed935Jews_1320nonJews_overlap_QC.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.


### List all SNPs in the original Jewish dataset


```bash
%%bash
../plink --bfile Jews_Genetics_Project/pre_QC_935Jews_original_gt/pre_QC_935Jews_original_gt --write-snplist \
--out original_935Jews
```

### Check how many SNPs in the post-QC merged dataset are present in the original Jewish dataset

```bash
%%bash
cat merged_imputed935Jews_1320nonJews_overlap_QC.bim|cut -d $'\t' -f 2 > \
merged_imputed935Jews_1320nonJews_overlap_QC.snplist &&
grep -Fxf merged_imputed935Jews_1320nonJews_overlap_QC.snplist original_935Jews.snplist|wc -l
```

    218186


### Generate genetic maps from <code>.bim</code> files by removing alleles and converting Morgan to cM


```python
import pandas as pd

bim_df = pd.read_csv("merged_imputed935Jews_1320nonJews_overlap_QC.bim", delimiter='\t', header=None)
bim_df = bim_df.drop(bim_df[[4, 5]], axis=1)
bim_df[2] = bim_df[2].apply(lambda morgan: round(float(morgan)*100, 4))
row_sum = 0
for chro in range(1, 23):
    bim_chr_df = bim_df.loc[bim_df[0]==chro]
    row_count, col_count = bim_chr_df.shape
    row_sum+=row_count
    new_map_name = f'merged_imputed935Jews_1320nonJews_overlap_QC_chr{chro}.map'
    bim_chr_df.to_csv(new_map_name, sep='\t', index=False, header=False)
print(row_sum)
```

    455546

### Split the post-QC merged dataset by chromosome

```bash
%%bash
for chr in {1..22}; do
  ../plink --silent --bfile merged_imputed935Jews_1320nonJews_overlap_QC --allow-no-sex \
  --recode vcf --chr $chr --out merged_imputed935Jews_1320nonJews_overlap_QC_chr${chr} &&
  bgzip merged_imputed935Jews_1320nonJews_overlap_QC_chr${chr}.vcf &&
  tabix -f merged_imputed935Jews_1320nonJews_overlap_QC_chr${chr}.vcf.gz
done
```


<a name="sec4"></a>
# Section IV: IBD inference by <code>phasedibd</code> and Extract of IBD segments


### Reference: https://github.com/23andMe/phasedibd 
### https://www.biorxiv.org/content/10.1101/2020.09.14.296939v1.full.pdf
### https://www.biorxiv.org/content/10.1101/2020.09.14.296939v1

```python
cd ../
```

```bash
%%bash
git clone https://github.com/23andMe/phasedibd
```

    Cloning into 'phasedibd'...
    remote: Enumerating objects: 106, done.[K
    remote: Counting objects: 100% (106/106), done.[K
    remote: Compressing objects: 100% (72/72), done.[K
    remote: Total 106 (delta 41), reused 85 (delta 30), pack-reused 0[K
    Receiving objects: 100% (106/106), 16.47 MiB | 18.02 MiB/s, done.
    Resolving deltas: 100% (41/41), done.


<a name="run-phasedibd"></a>
```python
cd phasedibd/
```

```bash
%%bash
make
python setup.py install 
python tests/unit_tests.py
```


```python
cd phasedibd/
```

### merged_imputed935Jews_1305nonJews_overlap_QC: 455546 variants

```bash
%%bash
for chr in {1..22}; do \
  gunzip -k ../../1240K_HO_pop_of_interest/merged_imputed935Jews_1320nonJews_overlap_QC_chr${chr}.vcf.gz
done
```

### According to the paper, <code>L_f=3</code> and <code>L_m=200</code> were used for the case study (empirical analysis of geographic patterns of haplotype sharing within Mexico)


```python
import phasedibd as ibd
from multiprocessing import Pool

def run_chr_phasedibd(chro):
    prefix_f = f"../../1240K_HO_pop_of_interest/merged_imputed935Jews_1320nonJews_overlap_QC_chr{chro}"
    tpbwt = ibd.TPBWTAnalysis()
    haplotypes_f = ibd.VcfHaplotypeAlignment(f"{prefix_f}.vcf", f"{prefix_f}.map")
    tpbwt.compute_ibd(haplotypes_f, compressed_out_path=f'compressed_haplotypes_1320_chr{chro}/', \
        chromosome=chro, segments_out_path=f"imputed935Jews_1320nonJews_L_m-200_chr{chro}", L_f=3, L_m=200)
    
if __name__ == '__main__':
    # calculate 22 autosomes in parallel
    pool = Pool(22)
    pool.map(run_chr_phasedibd, [str(chro) for chro in range(1, 23)])
    pool.close()
    pool.join()
```

```bash
%%bash
rm -rf compressed_haplotypes_1320_chr*
```

### Rename all these output .csv files


```bash
%%bash
ls imputed*/*.csv|wc -l
ls imputed*/*.csv|while read csv; do 
  newname=$(echo $csv|sed -r 's/\/[0-9]+\.csv/\.csv/') &&
  mv $csv $newname;
done
rmdir imputed935Jews*
ls *.csv
```

### In order to find out the individuals correctly, we would like to load <code>.fam</code> file of two merged dataset and switch <code>id1</code> and <code>id2</code> in the dfs from integers to their real <code>f'{fid}_{iid}'</code>


```python
df_aj_1320_fam = pd.read_csv("../../1240K_HO_pop_of_interest/merged_imputed935Jews_1320nonJews_overlap_QC.fam", delimiter=" ", header=None)
## The first column is $fid, the second is $iid
aj_1320_ids = df_aj_1320_fam[0].str.cat(df_aj_1320_fam[1].astype(str), sep="_").tolist()
```

### Therefore, let's write a script to reformat the output <code>.csv</code> files by converting the IDs, adding a column 'cm_len' and removing same-person IBD segments


```bash
%%bash
mkdir -p imputed935Jews_1320nonJews_reformatted
```

```python
import pandas as pd

def convert_ibd_csv(ibd_f: str, id_list: list, out_fname):
    ibd_df = pd.read_csv(ibd_f)
    ## Remove same-person IBD segments
    ibd_df = ibd_df.drop(ibd_df.loc[ibd_df['id1']==ibd_df['id2']].index).reset_index()
    ibd_df = ibd_df.drop(['index'], axis=1)
    ibd_df['cm_len'] = round(ibd_df['end_cm'] - ibd_df['start_cm'], 4)
    ibd_df['id1'] = ibd_df['id1'].apply(lambda el: id_list[el])
    ibd_df['id2'] = ibd_df['id2'].apply(lambda el: id_list[el])
    ibd_df.to_csv(out_fname, index=False)

if __name__ == "__main__": 
    df_aj_1320_fam = pd.read_csv("../../1240K_HO_pop_of_interest/merged_imputed935Jews_1320nonJews_overlap_QC.fam", delimiter=" ", header=None)
    ## The first column is $fid, the second is $iid
    aj_1320_ids = df_aj_1320_fam[0].str.cat(df_aj_1320_fam[1].astype(str), sep="_").tolist()
    for chro in range(1, 23):
        chro = str(chro)
        ibd_f = f"imputed935Jews_1320nonJews/imputed935Jews_1320nonJews_L_m-200_chr{chro}.csv"
        out_fname= f"imputed935Jews_1320nonJews_reformatted/imputed935Jews_1320nonJews_L_m-200_chr{chro}.csv"
        convert_ibd_csv(ibd_f, aj_1320_ids, out_fname)
```

### Gzip reformatted <code>.csv</code> files into a <code>.tar.gz</code> file


```bash
%%bash
tar -cz imputed935Jews_1320nonJews_reformatted/ > imputed935Jews_1320nonJews_reformatted.tar.gz &&
du -sh *.tar.gz
```

    64M	imputed935Jews_1320nonJews_reformatted.tar.gz
    47M	imputed935Jews_1320nonJews.tar.gz


```bash
%%bash
mv *_reformatted ..
```


```python
cd ../
```

```bash
%%bash
ls ../1240K_HO_pop_of_interest/*.id|grep -vE "(HGDP|indiv)"|wc -l
ls ../1240K_HO_pop_of_interest/*.id|grep -vE "(HGDP|indiv)" > 76pop_non_jews_id_file_list.txt
```

    76


<a name="jew-qc"></a>
### Removing replicate individuals and close relatives(minimally half-cousins) from Jewish individuals

#### Extract IBD within each Jewish community

```bash
%%bash
mkdir -p Jewish_pop_intra_IBD
```


```python
from multiprocessing import Pool
import pandas as pd
import csv 
import numpy as np
from os import path
from functools import partial
import scipy.stats as st
import warnings
import re

'''
There're 935 modern Jewish individuals in total, including 511 Ashkenazic Jews(AJ) and 424 other Jews(OJ).
This script presents average IBD segment length per pair in Centimorgans.
'''

def write_jew_intra_pop_pair_ibd(jew_pop: tuple, jew_geo_id_dict: dict):
    ## jew_pop is a tuple of two elements. 
    ### The first element is the real name of Jewish population. 
    ### The second element is its key at jew_geo_id_dict
    jew_pop_indivs = jew_geo_id_dict[jew_pop[1]]
    ## obtain the number of individuals from this Jewish population
    jew_pop_size = len(jew_pop_indivs)
    ## If n<2, we should not calculate intra-group IBD 
    if jew_pop_size < 2: return
    ## exclude IBD between two haplotypes of the same individual
    num_indiv_pairs = jew_pop_size*(jew_pop_size-1)/2
    jew_indiv_pair_ibd_dict = {}
    for i in range(jew_pop_size):
        for j in range(i+1, jew_pop_size):
            iid1, iid2 = jew_pop_indivs[i], jew_pop_indivs[j]
            jew_indiv_pair_ibd_dict[f'{iid1}__{iid2}'] = 0
    ## Create a dataframe to store sums of intra-population individual pairwise IBD length
    ### Output the dataframe into a csv file later
    jew_pop_intra_df = pd.DataFrame(columns=['IBD_cM'], index=list(jew_indiv_pair_ibd_dict.keys()))
    jew_pop_intra_df['IBD_cM'] = 0
    #print(jew_pop_intra_df.index)
    ## Loop through 22 autosomes
    for chr in range(1, 23):
        ### Sum up IBD segment lengths chromosome-wise and obtain the mean value finally after loop
        ibd_file = f'imputed935Jews_1320nonJews_reformatted/imputed935Jews_1320nonJews_L_m-200_chr{chr}.csv'
        if not path.isfile(ibd_file): continue
        df = pd.read_csv(ibd_file, header=0)
        pop_intra_df = df.loc[\
                (df['id1'].isin({iid: "" for iid in jew_pop_indivs})) & 
                (df['id2'].isin({iid: "" for iid in jew_pop_indivs})) 
            ]
        ibd_row_num = pop_intra_df.shape[0]
        for i in range(ibd_row_num):
            ibd_row = pop_intra_df.iloc[i]
            # id is in the form of f'{fid}_{iid}', $fid is a digit
            iid1=ibd_row['id1']
            iid2=ibd_row['id2']
            ibd_len = ibd_row['cm_len']
            #print(iid1, iid2)
            if f'{iid1}__{iid2}' in jew_indiv_pair_ibd_dict:
                #print("matched")
                jew_pop_intra_df.loc[f'{iid1}__{iid2}', 'IBD_cM']+=round(float(ibd_len), 5)
            elif f'{iid2}__{iid1}' in jew_indiv_pair_ibd_dict:
                #print("matched")
                jew_pop_intra_df.loc[f'{iid2}__{iid1}', 'IBD_cM']+=round(float(ibd_len), 5)
    ## Write to the general CSV file
    with open("Jew_subpops_intra_mean_pairwise_ibd_length.csv", "a") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        mean_v = round(jew_pop_intra_df['IBD_cM'].sum()/num_indiv_pairs, 3)
        ibd_csv_writer.writerow([jew_pop[0], str(mean_v), jew_pop_size, \
                                 round(st.sem(jew_pop_intra_df['IBD_cM'].tolist()), 3)])
    ## Write to an individual CSV file
    ## "jew_pop" is a tuple in the form like ("Berlin", "Germany", nan, nan), 
    ### so the second element annotates country
    jew_pop_intra_df.round(decimals=4).to_csv(f"Jewish_pop_intra_IBD/{jew_pop[0]}_intraIBD.csv", encoding='utf-8')
    

if __name__ == '__main__':
    # warnings.filterwarnings("ignore")
    # load 935 Jews into variables
    df_jews = pd.read_csv('../Jews_Genetics_Project/pre_QC_935_Jews_annotation.csv', header=0)
    df_aj = df_jews.loc[lambda x: x['data_source'].str.contains('AJ$', regex = True)]
    df_otherj = df_jews.loc[lambda x: x['data_source'].str.contains('otherJ$', regex = True)]
    ## add single-origin French Jews to df_otherJ
    df_FrenchJ = df_jews.loc[(df_jews['country1']=="France") & (df_jews['country2'].isnull())]
    ## add single-origin Italian Jews to df_otherJ
    df_ItalianJ = df_jews.loc[df_jews['country1']=="Italy"]
    ## Concat three dfs and drop duplicates
    df_otherj = pd.concat([df_otherj, df_FrenchJ, df_ItalianJ]).drop_duplicates()
    ## Create two dictionaries according to unique geographic origin
    df_aj_geo_groups= df_aj.groupby(['city1', 'country1', 'city2', 'country2']).groups
    df_oj_geo_groups= df_otherj.groupby(['city1', 'country1']).groups
    ##########################
    # print(df_oj_geo_groups)
    
    aj_geo_id_dict, oj_geo_id_dict = {}, {}
    
    ## non-mixed & unknown AJs, city2 & country2 must be nan
    for key, val in df_aj_geo_groups.items():
        for ind in val:
            fid=df_jews['family_id'].iloc[ind]
            iid=df_jews['individual_id'].iloc[ind]
            cat_id = "_".join((fid, iid))
            ## If non nan are present in key(a tuple of 4 elements) it's a mixed indiv, pass
            if not np.nan in key: 
                pass
            ## only NaN in key, thus being geographically unknown
            #elif len(set(key)) == 1:
                #pass
            elif not "Italy" in key and not "France" in key:
                if key not in aj_geo_id_dict: aj_geo_id_dict[key]=[]  
                aj_geo_id_dict[key].append(cat_id)
    
    ## Other-Jews
    for key, val in df_oj_geo_groups.items():
        for ind in val:
            ## avoid mixed individuals
            if not np.nan in key: 
                pass
            fid=df_jews['family_id'].iloc[ind]
            iid=df_jews['individual_id'].iloc[ind]
            if key not in oj_geo_id_dict: oj_geo_id_dict[key]=[]
            cat_id = "_".join((fid, iid))
            oj_geo_id_dict[key].append(cat_id)
    
    ## use this variable to store population names for all Jews
    jewish_pops = []
    for key in aj_geo_id_dict:
        pop_el = (f'{key[1]}AJ', key)
        jewish_pops.append(pop_el)
    for key in oj_geo_id_dict:
        if key[0] == "Mosul":
            pop_el = ('KurdistanOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Ar Raqqah':
            pop_el = ('SyrianKurdistanOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Mumbai':
            pop_el = ('IndiaMumbaiOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Thiruvananthapuram':
            pop_el = ('IndiaCochinOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Tbilisi':
            pop_el = ('GeorgiaOJ', key)
            jewish_pops.append(pop_el)
        else:
            pop_el = (f'{key[1]}OJ', key)
            jewish_pops.append(pop_el)
    
    #print(jewish_pops)
    aj_geo_id_dict.update(oj_geo_id_dict)
    jews_geo_dict = aj_geo_id_dict
    aj_geo_id_dict = None
    #print(jews_geo_dict)
    ## Output .csv file for mean pairwise IBD length
    with open("Jew_subpops_intra_mean_pairwise_ibd_length.csv", "w") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ibd_csv_writer.writerow(['pop1', 'pop2', 'mean_pairwise_IBD_length', \
                                 'pop1_size', 'pop2_size', 'IBD_len_SE'])
    
    ## For shorter computing time, 
    ## Run multi-processing using non-mixed subpop AJs and OJs
    partial_func = partial(write_jew_intra_pop_pair_ibd, jew_geo_id_dict=jews_geo_dict)
    ## the pool occupies 14 threads(14 sub-populations of AJ in parallel, excluding the Italian & French "AJ")
    ## There should be n*(n+1)/2 pairs within each group to compare IBD
    ## However, only one Moldavian AJ and one Slovakian AJ are present. 
    ### Therefore, it does not make sense to investigate these two groups
    pool = Pool(20)
    pool.map(partial_func, jewish_pops)
    pool.close()
    pool.join()
```

### We found some replicate individuals as well as suspected close relatives within same groups of Jews

### Select 448.75 cM as the threshold for filtering since any individual pairs that are IBD at > 448.75 cM are considered to be at least first cousins.


```bash
%%bash
ls Jewish_pop_intra_IBD/*.csv|while read csv; do
  pop=$(echo $csv|cut -d "/" -f 2|cut -d "_" -f 1);
  cat $csv|while read pair; do
    awk -v csv="$csv" -v pair="$pair" -v pop="$pop" -v FS="," '{if(int($2) > 448.75){print $0, pop}}';
  done
done|sort
```

    0_pop_11903_1118__TunisiaJew1118_TunisiaJew1118,5456.9549 TunisiaOJ
    0_pop_11903_1421__TunisiaJew1421_TunisiaJew1421,5447.2472 TunisiaOJ
    0_pop_11903_1511__TunisiaJew1511_TunisiaJew1511,5625.7502 TunisiaOJ
    0_pop_11903_1544__TunisiaJew1544_TunisiaJew1544,5521.0044 TunisiaOJ
    0_pop_11903_5200__TunisiaJew5200_TunisiaJew5200,5672.7867 TunisiaOJ
    0_pop_11907_5159__sephardic17bul_sephardic17bul,7816.7279 BulgariaOJ
    0_pop_11910_1599__eth_jew14_eth_jew14,5441.2251 EthiopiaOJ
    0_pop_11910_1599__eth_jew4_eth_jew4,5443.6669 EthiopiaOJ
    0_pop_11910_1613__eth_jew5_eth_jew5,5420.3634 EthiopiaOJ
    0_pop_11910_1683__eth_jew6_eth_jew6,5468.0983 EthiopiaOJ
    0_pop_11910_1804__eth_jew7_eth_jew7,5470.8836 EthiopiaOJ
    0_pop_11910_1818__eth_jew8_eth_jew8,5487.5635 EthiopiaOJ
    0_pop_11910_4524__eth_jew10_eth_jew10,5476.3477 EthiopiaOJ
    0_pop_11910_4568__eth_jew11_eth_jew11,5486.7092 EthiopiaOJ
    0_pop_11910_4574__eth_jew12_eth_jew12,5415.701 EthiopiaOJ
    0_pop_11910_4656__eth_jew13_eth_jew13,5497.0591 EthiopiaOJ
    0_pop_11915_1161__iran_jew1_iran_jew1,7145.2153 IranOJ
    0_pop_11915_1409__IranJew1409_IranJew1409,5474.5211 IranOJ
    0_pop_11915_1419__IranJew1419_IranJew1419,5083.2672 IranOJ
    0_pop_11915_1557__IranJew1557_IranJew1557,5436.1558 IranOJ
    0_pop_11916_1418__iraq_jew1_iraq_jew1,7116.8662 IraqOJ
    0_pop_11916_1444__iraq_jew6_iraq_jew6,7279.6497 IraqOJ
    0_pop_11931_1953__Yemen_Jew_6_Yemen_Jew_6,7030.2813 YemenOJ
    0_pop_920_1551__KurdJew1551_KurdJew1551,5422.8385 KurdistanOJ
    0_pop_920_1824__KurdJew1824_KurdJew1824,5370.1578 KurdistanOJ
    0_pop_920_4584__KurdJew1580_KurdJew1580,5344.2976 KurdistanOJ
    0_pop_920_4584__KurdJew1592_KurdJew1592,606.723 KurdistanOJ
    0_pop_920_4633__KurdJew4633_KurdJew4633,5351.5627 KurdistanOJ
    0_pop_920_4663__KurdJew4663_KurdJew4663,5423.2492 KurdistanOJ
    0_pop_922_1579__LibyaJew1579_LibyaJew1579,5414.1643 LibyaOJ
    eth_jew4_eth_jew4__eth_jew14_eth_jew14,6911.4999 EthiopiaOJ
    KurdJew1580_KurdJew1580__KurdJew1592_KurdJew1592,879.0088 KurdistanOJ
    MorJew2001_MorJew2001__Morocco_Jew_1_Morocco_Jew_1,6872.4927 MoroccoOJ
    single-Belarus-sample_108__single-Belarus-sample_119,1763.672 BelarusAJ
    single-Poland-sample_061__single-Poland-sample_068,2133.6032 PolandAJ


### Priority for selection from a pair of duplicates judging the number of SNPs:
#### Behar et. al., 2010 ( ~570K SNPs ) > Kopelman et. al., 2020 ( ~470K SNPs ) > Behar et. al., 2013 ( ~300K SNPs )

#### For suspected close relatives, let's select the first sample at each pair


```bash
%%bash
echo -e "
family_id,individual_id
TunisiaJew1118,TunisiaJew1118
TunisiaJew1421,TunisiaJew1421
TunisiaJew1511,TunisiaJew1511
TunisiaJew1544,TunisiaJew1544
TunisiaJew5200,TunisiaJew5200
0,pop_11907_5159
0,pop_11910_1599
0,pop_11910_1599
0,pop_11910_1613
0,pop_11910_1683
0,pop_11910_1804
0,pop_11910_1818
0,pop_11910_4524
0,pop_11910_4568
0,pop_11910_4574
0,pop_11910_4656
0,pop_11915_1161
IranJew1409,IranJew1409
IranJew1419,IranJew1419
IranJew1557,IranJew1557
0,pop_11916_1418
0,pop_11916_1444
0,pop_11931_1953
KurdJew1551,KurdJew1551
KurdJew1824,KurdJew1824
KurdJew1580,KurdJew1580
KurdJew1592,KurdJew1592
KurdJew4633,KurdJew4633
KurdJew4663,KurdJew4663
LibyaJew1579,LibyaJew1579
eth_jew14,eth_jew14
KurdJew1592,KurdJew1592
MorJew2001,MorJew2001
single-Belarus-sample,119
single-Poland-sample,068
" > jews_to_delete.csv
```


```bash
%%bash
cat jews_to_delete.csv|sort|uniq|wc -l
wc -l ../Jews_Genetics_Project/pre_QC_935_Jews_annotation.csv
## check the number of nique individual IDs
cat ../Jews_Genetics_Project/pre_QC_935_Jews_annotation.csv|cut -d , -f 3|sort|uniq|wc -l
```

    35
    936 ../Jews_Genetics_Project/pre_QC_935_Jews_annotation.csv
    936

```python
import pandas as pd

jews_935_df = pd.read_csv("../Jews_Genetics_Project/pre_QC_935_Jews_annotation.csv")
print(jews_935_df.shape)
jews_to_delete_df = pd.read_csv("jews_to_delete.csv")
iids_to_delete_dict = {v: k for k, v in dict(jews_to_delete_df['individual_id']).items()}
print(len(iids_to_delete_dict.items()))
filtered_jews_df = jews_935_df.loc[~(jews_935_df['individual_id'].isin(iids_to_delete_dict))].reset_index()
filtered_jews_df = filtered_jews_df.drop(['index'], axis=1)
print(filtered_jews_df)
filtered_jews_df.to_csv("Bray2010_Behar1013_Gladstein2019_Kopelman2020_902Jews.csv", index=None)
```

    (935, 11)
    33
                    data_source     family_id individual_id  latitude1  \
    0      Kopelman_2020_148_AJ             0  pop_917_5750    41.8960   
    1    Kopelman2020_239otherJ             0   pop_920_155    36.3450   
    2    Kopelman2020_239otherJ             0   pop_920_174    36.3450   
    3    Kopelman2020_239otherJ             0   pop_920_209    36.3450   
    4    Kopelman2020_239otherJ             0  pop_920_1551    36.3450   
    ..                      ...           ...           ...        ...   
    897   Behar_2010_100_otherJ  Yemen_Jew_11  Yemen_Jew_11    15.3547   
    898   Behar_2010_100_otherJ  Yemen_Jew_12  Yemen_Jew_12    15.3547   
    899   Behar_2010_100_otherJ  Yemen_Jew_13  Yemen_Jew_13    15.3547   
    900   Behar_2010_100_otherJ  Yemen_Jew_14  Yemen_Jew_14    15.3547   
    901   Behar_2010_100_otherJ  Yemen_Jew_15  Yemen_Jew_15    15.3547   
    
         longtitude1  city1 country1  latitude2  longtitude2 city2 country2  
    0        12.4833   Rome    Italy        NaN          NaN   NaN      NaN  
    1        43.1450  Mosul     Iraq        NaN          NaN   NaN      NaN  
    2        43.1450  Mosul     Iraq        NaN          NaN   NaN      NaN  
    3        43.1450  Mosul     Iraq        NaN          NaN   NaN      NaN  
    4        43.1450  Mosul     Iraq        NaN          NaN   NaN      NaN  
    ..           ...    ...      ...        ...          ...   ...      ...  
    897      44.2066  Sanaa    Yemen        NaN          NaN   NaN      NaN  
    898      44.2066  Sanaa    Yemen        NaN          NaN   NaN      NaN  
    899      44.2066  Sanaa    Yemen        NaN          NaN   NaN      NaN  
    900      44.2066  Sanaa    Yemen        NaN          NaN   NaN      NaN  
    901      44.2066  Sanaa    Yemen        NaN          NaN   NaN      NaN  
    
    [902 rows x 11 columns]



```bash
%%bash
head -n5 Bray2010_Behar1013_Gladstein2019_Kopelman2020_902Jews.csv
```

    data_source,family_id,individual_id,latitude1,longtitude1,city1,country1,latitude2,longtitude2,city2,country2
    Kopelman_2020_148_AJ,0,pop_917_5750,41.896,12.4833,Rome,Italy,,,,
    Kopelman2020_239otherJ,0,pop_920_155,36.345,43.145,Mosul,Iraq,,,,
    Kopelman2020_239otherJ,0,pop_920_174,36.345,43.145,Mosul,Iraq,,,,
    Kopelman2020_239otherJ,0,pop_920_209,36.345,43.145,Mosul,Iraq,,,,


### Therefore, 33 individuals, including 2 AJs and 31 other Jews were filtered out. 
### After filtering, we have 288 single-origin AJs, 221 mixed-origin AJs and 393 other Jews, i.e. 509 AJs and 393 OJs.

### In addition, let's write two scripts to extract IBD sharing between different AJ sub-populations and between different OJ sub-populations.

```bash
%%bash
mkdir -p aj_pops_inter_IBD -p oj_pops_inter_IBD
```


```python
from multiprocessing import Pool
import pandas as pd
import csv 
import numpy as np
from os import path
from functools import partial
import scipy.stats as st
import warnings
import re

'''
There're 902 modern Jewish individuals in total, including 509 Ashkenazic Jews(AJ) and 393 other Jews(other_aj).
This script presents average IBD segment length per pair in Centimorgans.
'''

def write_aj_aj_pop_pair_ibd(aj_pop: tuple, aj_geo_id_dict: dict, aj_id_geo_dict: dict):
    ## aj_pop is a tuple of two elements. 
    ### The first element is the real name of Jewish population. 
    ### The second element is its key at aj_geo_id_dict
    aj_pop_indivs = aj_geo_id_dict[aj_pop[1]]
    aj_pop_indivs_dict = {iid: "" for iid in aj_pop_indivs}
    ## obtain the number of individuals from this Jewish population
    aj_pop_size = len(aj_pop_indivs)
    other_aj_col_names = []
    for other_aj_id in aj_id_geo_dict.keys(): 
        ## the prefix should be country name = "AJ"
        other_aj_col_prefix = f'{aj_id_geo_dict[other_aj_id][1]}AJ'
        other_aj_col_names.append(f'{other_aj_col_prefix}__{other_aj_id}')
    aj_other_aj_df = pd.DataFrame(columns=aj_pop_indivs, index=other_aj_col_names)
    for aj_id in aj_pop_indivs:
        aj_other_aj_df[aj_id]=0
    # Use this dict the store the following key:value pair
    ## f'{AJ_SUBPOP_NAME}_{OTHER_SUBPOP_NAME}'(key): mean pairwise IBD length(float)
    ajsubpop_other_ajsubpop_mean_pws_ibd_len = {}
    ## Loop through 22 autosomes
    for chr in range(1, 23):
        ### Sum up IBD segment lengths chromosome-wise and obtain the mean value finally after loop
        ibd_file = f'imputed935Jews_1320nonJews_reformatted/imputed935Jews_1320nonJews_L_m-200_chr{chr}.csv'
        if not path.isfile(ibd_file): continue
        df = pd.read_csv(ibd_file, header=0)
        aj_subpop_ibd_df = df.loc[
            (df['id1'].isin(aj_pop_indivs_dict)) | (df['id2'].isin(aj_pop_indivs_dict)) 
        ]
        for other_aj_pop in aj_geo_id_dict.keys():
            #########################################################
            # Filter into the corresponding IBD segments
            if other_aj_pop == aj_pop[1]: continue
            other_aj_subpop_id_dict = {el: '' for el in aj_geo_id_dict[other_aj_pop]}
            filtered_df = aj_subpop_ibd_df.loc[
                (df['id1'].isin(other_aj_subpop_id_dict)) | (df['id2'].isin(other_aj_subpop_id_dict))
            ]
            other_aj_pop_size = len(aj_geo_id_dict[other_aj_pop])
            other_aj_pop_str = other_aj_pop[1] if np.nan not in other_aj_pop else 'unknown'
            pop_pair_str = f'{aj_pop[0]}_{other_aj_pop_str}AJ'
            ## use the condition "aj_pop[0] > other_aj_pop_str" to avoid replicate pairs, i.e. A & B and B & A
            if aj_pop_size > 0 and other_aj_pop_size > 0 and aj_pop[0] > other_aj_pop_str:
                ## Create a string to describe the AJ subgroup-OJ subgroup pair
                if pop_pair_str not in ajsubpop_other_ajsubpop_mean_pws_ibd_len:
                    ## value of the dict is initially a list, [sum_of_ibd_lengths, subpop_size1, subpop_size2]
                    ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str]=[0, aj_pop_size, other_aj_pop_size, {}]
                ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str][0]+=filtered_df['cm_len'].sum()
            ##################################################
            # Output the CSV file storing IBD length between each pair of individuals
            filtered_row_num = filtered_df.shape[0]
            #print(chr, pop_pair_str, filtered_row_num)
            for i in range(filtered_row_num):
                ibd_row = filtered_df.iloc[i]
                # id is in the form of f'{fid}_{iid}', $fid is a digit
                id1=ibd_row['id1']
                id2=ibd_row['id2']
                ibd_len = ibd_row['cm_len']
                ###Test
                #print(ibd_len)
                ###
                if id1 in aj_pop_indivs_dict and id2 in other_aj_subpop_id_dict:
                    ## the prefix should be country name + 'AJ'
                    id2_prefix = f'{aj_id_geo_dict[id2][1]}AJ'
                    aj_other_aj_df.loc[f'{id2_prefix}__{id2}', id1]+=round(float(ibd_len), 5)
                    ###Test
                    #print(aj_other_aj_df.loc[f'{id2_prefix}__{id2}', id1])
                    ###
                    indiv_pair_key = f'{id1}__{id2}'
                elif id2 in aj_pop_indivs_dict and id1 in other_aj_subpop_id_dict:
                    id1_prefix = f'{aj_id_geo_dict[id1][1]}AJ'
                    aj_other_aj_df.loc[f'{id1_prefix}__{id1}', id2]+=round(float(ibd_len), 5)
                    indiv_pair_key = f'{id2}__{id1}'
                if pop_pair_str in ajsubpop_other_ajsubpop_mean_pws_ibd_len:
                    #print(pop_pair_str)
                    #print(ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str][3].keys())
                    if not indiv_pair_key in ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str][3]:
                        ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]=0
                    ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]+=\
                        round(float(ibd_len), 5)
    ##########################################################################################
    ##########################################################################################
    # Finally, calculate the mean pairwise IBD length for AJ versus other AJs
    ## after looping through 22 autosomes
    ### This dict should have 14*48 items
    with open("AJ_subpop_inter_mean_pairwise_ibd_length.csv", "a+") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ### This dict should have 14*13/2 items
        for k, v in ajsubpop_other_ajsubpop_mean_pws_ibd_len.items():
            if type(v) is list:
                if len(v)==4:
                    num_non_zero_ibd_indiv_pairs = len(v[3].keys())
                    num_possible_pairs = int(v[1])*int(v[2])
                    num_zeros_to_add = num_possible_pairs - num_non_zero_ibd_indiv_pairs
                    v[3] = list(v[3].values())
                    if num_zeros_to_add > 0:
                        v[3].extend([0]*num_zeros_to_add)
                    pop1=k.split("_")[0]
                    pop2="_".join(tuple(k.split("_")[1:]))
                    mean_v = round(v[0]/num_possible_pairs, 3)
                    # v[3] is a list
                    #print(v[3])
                    ibd_csv_writer.writerow([pop1,pop2,str(mean_v),str(v[1]),str(v[2]),\
                                             round(st.sem(v[3]), 3)])
    #################################################################################################
    ## Additionally, output 14 .csv files storing total IBD length between every pairs of individuals
    ## after looping through 22 autosomes
    ### Reference for conversion from df to csv: 
    ####: https://towardsdatascience.com/how-to-export-pandas-dataframe-to-csv-2038e43d9c03
    #### Write to an individual CSV file
    #### "aj_pop" is a tuple in the form like ("Berlin", "Germany", nan, nan), 
    #### so the second element annotates country
    aj_other_aj_df.round(decimals=4).to_csv(f"aj_pops_inter_IBD/{aj_pop[0]}_IBD.csv", encoding='utf-8')
    

if __name__ == '__main__':
    # warnings.filterwarnings("ignore")
    # load 902 Jews into variables
    df_jews = pd.read_csv('Bray2010_Behar1013_Gladstein2019_Kopelman2020_902Jews.csv', header=0)
    df_aj = df_jews.loc[lambda x: x['data_source'].str.contains('AJ$', regex = True)]
    ## Create two dictionaries according to unique geographic origin
    df_aj_geo_groups= df_aj.groupby(['city1', 'country1', 'city2', 'country2']).groups
    aj_geo_id_dict, aj_id_geo_dict = {}, {}
    
    ## non-mixed & unknown AJs, city2 & country2 must be nan
    for key, val in df_aj_geo_groups.items():
        for ind in val:
            fid=df_jews['family_id'].iloc[ind]
            iid=df_jews['individual_id'].iloc[ind]
            cat_id = "_".join((fid, iid))
            ## If non nan are present in key(a tuple of 4 elements) it's a mixed indiv, pass
            if not np.nan in key: 
                pass
            ## only NaN in key, thus being geographically unknown
            #elif len(set(key)) == 1:
                #pass
            elif not "Italy" in key and not "France" in key:
                if key not in aj_geo_id_dict: aj_geo_id_dict[key]=[]  
                aj_geo_id_dict[key].append(cat_id)
                aj_id_geo_dict[cat_id]=key
    
    ## use this variable to store population names for all Jews
    jewish_pops = []
    for key in aj_geo_id_dict:
        pop_el = (f'{key[1]}AJ', key)
        jewish_pops.append(pop_el)
    
    #print(jewish_pops)
    ## Output .csv file for mean pairwise IBD length
    with open("AJ_subpop_inter_mean_pairwise_ibd_length.csv", "w") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ibd_csv_writer.writerow(['pop1', 'pop2', 'mean_pairwise_IBD_length', 'pop1_size', 'pop2_size', \
                                 'IBD_len_SE'])
    
    ## For shorter computing time, 
    ## Run multi-processing using non-mixed subpop AJs and other_ajs
    partial_func = partial(write_aj_aj_pop_pair_ibd, aj_geo_id_dict=aj_geo_id_dict, aj_id_geo_dict=aj_id_geo_dict)
    ## the pool occupies 14 threads(14 sub-populations of AJ in parallel, excluding the Italian & French "AJ")
    ## There should be 14*(14-1)/2 pairs of AJs to compare IBD
    ## However, only one Moldavian AJ and one Slovakian AJ are present. 
    pool = Pool(14)
    pool.map(partial_func, jewish_pops)
    pool.close()
    pool.join()

```


### Inter-group IBD sharing between 23 sub-groups of other Jews


```python
from multiprocessing import Pool
import pandas as pd
import csv 
import numpy as np
from os import path
from functools import partial
import scipy.stats as st
import warnings
import re

'''
There're 901 modern Jewish individuals in total, including 508 Ashkenazic Jews(OJ) and 393 other Jews(OJ).
This script presents average IBD segment length per pair in Centimorgans.
Version of refined-ibd: refined-ibd.17Jan20.102
Version of Java: Java/1.8.0_192
refined-ibd analysis was performed using the default parameters, e.g. minimum-lod=3
'''

def write_oj_oj_pop_pair_ibd(oj_pop: tuple, oj_geo_id_dict: dict, oj_id_geo_dict: dict):
    ## oj_pop is a tuple of two elements. 
    ### The first element is the real name of Jewish population. 
    ### The second element is its key at oj_geo_id_dict
    oj_pop_indivs = oj_geo_id_dict[oj_pop[1]]
    oj_pop_indivs_dict = {iid: "" for iid in oj_pop_indivs}
    ## obtain the number of individuals from this Jewish population
    oj_pop_size = len(oj_pop_indivs)
    other_oj_col_names = []
    for other_oj_id in oj_id_geo_dict.keys(): 
        ## the prefix should be country name = "OJ"
        if oj_id_geo_dict[other_oj_id] == oj_pop[1]: continue
        other_oj_col_prefix = f'{oj_id_geo_dict[other_oj_id][0]}-{oj_id_geo_dict[other_oj_id][1]}OJ'
        other_oj_col_names.append(f'{other_oj_col_prefix}__{other_oj_id}')
    oj_other_oj_df = pd.DataFrame(columns=oj_pop_indivs, index=other_oj_col_names)
    for oj_id in oj_pop_indivs:
        oj_other_oj_df[oj_id]=0
    # Use this dict the store the following key:value pair
    ## f'{OJ_SUBPOP_NAME}_{OTHER_SUBPOP_NAME}'(key): mean pairwise IBD length(float)
    ojsubpop_other_ojsubpop_mean_pws_ibd_len = {}
    ## Loop through 22 autosomes
    for chr in range(1, 23):
        ### Sum up IBD segment lengths chromosome-wise and obtain the mean value finally after loop
        ibd_file = f'imputed935Jews_1320nonJews_reformatted/imputed935Jews_1320nonJews_L_m-200_chr{chr}.csv'
        if not path.isfile(ibd_file): continue
        df = pd.read_csv(ibd_file, header=0)
        oj_subpop_ibd_df = df.loc[
            (df['id1'].isin(oj_pop_indivs_dict)) | (df['id2'].isin(oj_pop_indivs_dict)) 
        ]
        for other_oj_pop in oj_geo_id_dict.keys():
            #########################################################
            # Filter into the corresponding IBD segments
            if other_oj_pop == oj_pop[1]: continue
            other_oj_subpop_id_dict = {el: '' for el in oj_geo_id_dict[other_oj_pop]}
            filtered_df = oj_subpop_ibd_df.loc[
                (df['id1'].isin(other_oj_subpop_id_dict)) | (df['id2'].isin(other_oj_subpop_id_dict))
            ]
            other_oj_pop_size = len(oj_geo_id_dict[other_oj_pop])
            oj_pop_str = oj_pop[1][0]
            other_oj_pop_str = other_oj_pop[0] if np.nan not in other_oj_pop else 'unknown'
            pop_pair_str = f'{oj_pop[0]}_{other_oj_pop_str}OJ'
            ## use the condition "oj_pop[0] > other_oj_pop_str" to avoid replicate pairs, i.e. A & B and B & A
            if oj_pop_size > 0 and other_oj_pop_size > 0 and oj_pop_str > other_oj_pop_str:
                ## Create a string to describe the OJ subgroup-OJ subgroup pair
                if pop_pair_str not in ojsubpop_other_ojsubpop_mean_pws_ibd_len:
                    ## value of the dict is initially a list, [sum_of_ibd_lengths, subpop_size1, subpop_size2]
                    ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str]=[0, oj_pop_size, other_oj_pop_size, {}]
                ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str][0]+=filtered_df['cm_len'].sum()
            ##################################################
            # Output the CSV file storing IBD length between each pair of individuals
            filtered_row_num = filtered_df.shape[0]
            #print(chr, pop_pair_str, filtered_row_num)
            for i in range(filtered_row_num):
                ibd_row = filtered_df.iloc[i]
                # id is in the form of f'{fid}_{iid}', $fid is a digit
                id1=ibd_row['id1']
                id2=ibd_row['id2']
                ibd_len = ibd_row['cm_len']
                ###Test
                #print(ibd_len)
                ###
                if id1 in oj_pop_indivs_dict and id2 in other_oj_subpop_id_dict:
                    ## the prefix should be country name + 'OJ'
                    id2_prefix = f'{oj_id_geo_dict[id2][0]}-{oj_id_geo_dict[id2][1]}OJ'
                    oj_other_oj_df.loc[f'{id2_prefix}__{id2}', id1]+=round(float(ibd_len), 5)
                    ###Test
                    #print(oj_other_oj_df.loc[f'{id2_prefix}__{id2}', id1])
                    ###
                    indiv_pair_key = f'{id1}__{id2}'
                elif id2 in oj_pop_indivs_dict and id1 in other_oj_subpop_id_dict:
                    id1_prefix = f'{oj_id_geo_dict[id1][0]}-{oj_id_geo_dict[id1][1]}OJ'
                    oj_other_oj_df.loc[f'{id1_prefix}__{id1}', id2]+=round(float(ibd_len), 5)
                    indiv_pair_key = f'{id2}__{id1}'
                if pop_pair_str in ojsubpop_other_ojsubpop_mean_pws_ibd_len:
                    #print(pop_pair_str)
                    #print(ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str][3].keys())
                    if not indiv_pair_key in ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str][3]:
                        ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]=0
                    ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]+=\
                        round(float(ibd_len), 5)
    ##########################################################################################
    ##########################################################################################
    # Finally, calculate the mean pairwise IBD length for OJ versus other OJs
    ## after looping through 22 autosomes
    ### This dict should have 14*48 items
    with open("OJ_subpop_inter_mean_pairwise_ibd_length.csv", "a+") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ### This dict should have 14*13/2 items
        for k, v in ojsubpop_other_ojsubpop_mean_pws_ibd_len.items():
            if type(v) is list:
                if len(v)==4:
                    num_non_zero_ibd_indiv_pairs = len(v[3].keys())
                    num_possible_pairs = int(v[1])*int(v[2])
                    num_zeros_to_add = num_possible_pairs - num_non_zero_ibd_indiv_pairs
                    v[3] = list(v[3].values())
                    if num_zeros_to_add > 0:
                        v[3].extend([0]*num_zeros_to_add)
                    pop1=k.split("_")[0]
                    pop2="_".join(tuple(k.split("_")[1:]))
                    mean_v = round(v[0]/num_possible_pairs, 3)
                    # v[3] is a list
                    #print(v[3])
                    ibd_csv_writer.writerow([pop1,pop2,str(mean_v),str(v[1]),str(v[2]),\
                                             round(st.sem(v[3]), 3)])
    #################################################################################################
    ## Additionally, output 14 .csv files storing total IBD length between every pairs of individuals
    ## after looping through 22 autosomes
    ### Reference for conversion from df to csv: 
    ####: https://towardsdatascience.com/how-to-export-pandas-dataframe-to-csv-2038e43d9c03
    #### Write to an individual CSV file
    #### "oj_pop" is a tuple in the form like ("Berlin", "Germany", nan, nan), 
    #### so the second element annotates country
    oj_other_oj_df.round(decimals=4).to_csv(f"oj_pops_inter_IBD/{oj_pop[0]}_IBD.csv", encoding='utf-8')
    

if __name__ == '__main__':
    # warnings.filterwarnings("ignore")
    # load 902 Jews into variables
    df_jews = pd.read_csv('Bray2010_Behar1013_Gladstein2019_Kopelman2020_901Jews.csv', header=0)
    df_otherj = df_jews.loc[lambda x: x['data_source'].str.contains('otherJ$', regex = True)]
    ## add single-origin French Jews to df_otherJ
    df_FrenchJ = df_jews.loc[(df_jews['country1']=="France") & (df_jews['country2'].isnull())]
    ## add single-origin Italian Jews to df_otherJ
    df_ItalianJ = df_jews.loc[df_jews['country1']=="Italy"]
    ## Concat three dfs and drop duplicates
    df_otherj = pd.concat([df_otherj, df_FrenchJ, df_ItalianJ]).drop_duplicates()
    ## Create two dictionaries according to unique geographic origin
    df_oj_geo_groups= df_otherj.groupby(['city1', 'country1']).groups
    ##########################
    # print(df_oj_geo_groups)
    oj_geo_id_dict, oj_id_geo_dict = {}, {}
    
    ## Other-Jews
    for key, val in df_oj_geo_groups.items():
        for ind in val:
            ## avoid mixed individuals
            if not np.nan in key: 
                pass
            fid=df_jews['family_id'].iloc[ind]
            iid=df_jews['individual_id'].iloc[ind]
            if key not in oj_geo_id_dict: oj_geo_id_dict[key]=[]
            cat_id = "_".join((fid, iid))
            oj_geo_id_dict[key].append(cat_id)
            oj_id_geo_dict[cat_id]=key
    
    jewish_pops = []
    ## use this variable to store population names for all Jews
    for key in oj_geo_id_dict:
        if key[0] == "Mosul":
            pop_el = ('KurdistanOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Ar Raqqah':
            pop_el = ('SyrianKurdistanOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Mumbai':
            pop_el = ('IndiaMumbaiOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Thiruvananthapuram':
            pop_el = ('IndiaCochinOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Tbilisi':
            pop_el = ('GeorgiaOJ', key)
            jewish_pops.append(pop_el)
        else:
            pop_el = (f'{key[1]}OJ', key)
            jewish_pops.append(pop_el)
    
    #print(jewish_pops)
    #print(jews_geo_dict)
    ## Output .csv file for mean pairwise IBD length
    with open("OJ_subpop_inter_mean_pairwise_ibd_length.csv", "w") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ibd_csv_writer.writerow(['pop1', 'pop2', 'mean_pairwise_IBD_length', 'pop1_size', 'pop2_size', \
                                 'IBD_len_SE'])
    
    ## For shorter computing time, 
    ## Run multi-processing using non-mixed subpop OJs
    partial_func = partial(write_oj_oj_pop_pair_ibd, oj_geo_id_dict=oj_geo_id_dict, oj_id_geo_dict=oj_id_geo_dict)
    ## the pool occupies 23 threads(23 sub-populations of OJ in parallel, including the Italian & French "AJ")
    ## There should be 23*(23-1)/2 pairs of OJs to compare IBD
    pool = Pool(23)
    pool.map(partial_func, jewish_pops)
    pool.close()
    pool.join()
```

### From the analysis above, we found that an "Austro-Hungarian" AJ "ashkenazy1w_ashkenazy1w" from Behar et.al., 2010 is actually the same person as an "Hungarian" AJ "0_pop_11934_5794" from Kopelman et.al., 2020

### Therefore, we follow the earlier annotation from Behar et. al., and remove "0_pop_11934_5794" from our dataset.


```bash
%%bash
echo -e "
family_id,individual_id
TunisiaJew1118,TunisiaJew1118
TunisiaJew1421,TunisiaJew1421
TunisiaJew1511,TunisiaJew1511
TunisiaJew1544,TunisiaJew1544
TunisiaJew5200,TunisiaJew5200
0,pop_11907_5159
0,pop_11910_1599
0,pop_11910_1599
0,pop_11910_1613
0,pop_11910_1683
0,pop_11910_1804
0,pop_11910_1818
0,pop_11910_4524
0,pop_11910_4568
0,pop_11910_4574
0,pop_11910_4656
0,pop_11915_1161
IranJew1409,IranJew1409
IranJew1419,IranJew1419
IranJew1557,IranJew1557
0,pop_11916_1418
0,pop_11916_1444
0,pop_11931_1953
KurdJew1551,KurdJew1551
KurdJew1824,KurdJew1824
KurdJew1580,KurdJew1580
KurdJew1592,KurdJew1592
KurdJew4633,KurdJew4633
KurdJew4663,KurdJew4663
LibyaJew1579,LibyaJew1579
eth_jew14,eth_jew14
KurdJew1592,KurdJew1592
MorJew2001,MorJew2001
single-Belarus-sample,119
single-Poland-sample,068
0,pop_11934_5794
" > jews_to_delete.csv
```

```python
import pandas as pd

jews_935_df = pd.read_csv("../Jews_Genetics_Project/pre_QC_935_Jews_annotation.csv")
print(jews_935_df.shape)
jews_to_delete_df = pd.read_csv("jews_to_delete.csv")
iids_to_delete_dict = {v: k for k, v in dict(jews_to_delete_df['individual_id']).items()}
print(len(iids_to_delete_dict.items()))
filtered_jews_df = jews_935_df.loc[~(jews_935_df['individual_id'].isin(iids_to_delete_dict))].reset_index()
filtered_jews_df = filtered_jews_df.drop(['index'], axis=1)
print(filtered_jews_df)
filtered_jews_df.to_csv("Bray2010_Behar1013_Gladstein2019_Kopelman2020_901Jews.csv", index=None)
```

    (935, 11)
    34
                    data_source     family_id individual_id  latitude1  \
    0      Kopelman_2020_148_AJ             0  pop_917_5750    41.8960   
    1    Kopelman2020_239otherJ             0   pop_920_155    36.3450   
    2    Kopelman2020_239otherJ             0   pop_920_174    36.3450   
    3    Kopelman2020_239otherJ             0   pop_920_209    36.3450   
    4    Kopelman2020_239otherJ             0  pop_920_1551    36.3450   
    ..                      ...           ...           ...        ...   
    896   Behar_2010_100_otherJ  Yemen_Jew_11  Yemen_Jew_11    15.3547   
    897   Behar_2010_100_otherJ  Yemen_Jew_12  Yemen_Jew_12    15.3547   
    898   Behar_2010_100_otherJ  Yemen_Jew_13  Yemen_Jew_13    15.3547   
    899   Behar_2010_100_otherJ  Yemen_Jew_14  Yemen_Jew_14    15.3547   
    900   Behar_2010_100_otherJ  Yemen_Jew_15  Yemen_Jew_15    15.3547   
    
         longtitude1  city1 country1  latitude2  longtitude2 city2 country2  
    0        12.4833   Rome    Italy        NaN          NaN   NaN      NaN  
    1        43.1450  Mosul     Iraq        NaN          NaN   NaN      NaN  
    2        43.1450  Mosul     Iraq        NaN          NaN   NaN      NaN  
    3        43.1450  Mosul     Iraq        NaN          NaN   NaN      NaN  
    4        43.1450  Mosul     Iraq        NaN          NaN   NaN      NaN  
    ..           ...    ...      ...        ...          ...   ...      ...  
    896      44.2066  Sanaa    Yemen        NaN          NaN   NaN      NaN  
    897      44.2066  Sanaa    Yemen        NaN          NaN   NaN      NaN  
    898      44.2066  Sanaa    Yemen        NaN          NaN   NaN      NaN  
    899      44.2066  Sanaa    Yemen        NaN          NaN   NaN      NaN  
    900      44.2066  Sanaa    Yemen        NaN          NaN   NaN      NaN  
    
    [901 rows x 11 columns]



```bash
%%bash
head -n5 Bray2010_Behar1013_Gladstein2019_Kopelman2020_901Jews.csv &&
rm Bray2010_Behar1013_Gladstein2019_Kopelman2020_902Jews.csv
```

    data_source,family_id,individual_id,latitude1,longtitude1,city1,country1,latitude2,longtitude2,city2,country2
    Kopelman_2020_148_AJ,0,pop_917_5750,41.896,12.4833,Rome,Italy,,,,
    Kopelman2020_239otherJ,0,pop_920_155,36.345,43.145,Mosul,Iraq,,,,
    Kopelman2020_239otherJ,0,pop_920_174,36.345,43.145,Mosul,Iraq,,,,
    Kopelman2020_239otherJ,0,pop_920_209,36.345,43.145,Mosul,Iraq,,,,


### Therefore, 33 individuals, including 2 AJs and 31 other Jews were filtered out. 
### After filtering, we have 288 single-origin AJs, 221 mixed-origin AJs and 393 other Jews, i.e. 509 AJs and 393 OJs.


### Hence, let's run IBD analyses using the new annotation file consisting of 901 Jews.
### Additionally,  we would like to extract specified inter-population or intra-population IBD segments at 22 autosomes to a file for each pair.


#### Extracting IBD between AJs and non-Jews
```bash
%%bash
mkdir -p IBD_segments_AJ_subpop_vs_76_nonJ_pop
```


```python
from multiprocessing import Pool
import pandas as pd
import csv 
import numpy as np
from os import path
from functools import partial
import scipy.stats as st
import warnings
import re

'''
There're 901 modern Jewish individuals in total, including 508 Ashkenazic Jews(AJ) and 393 other Jews(OJ).
This script presents average IBD segment length per pair in Centimorgans.
'''

def write_aj_otherpop_pair_ibd(aj_subpop: tuple, aj_geo_id_dict: dict, nj_pop_id_dict: dict, \
                              aj_id_geo_dict: dict, nj_id_pop_dict: dict):
    # One CSV file as output from a pandas.dataframe per AJ subpop (14*2 CSV files in total)
    ## Two dataframes to store total IBD length between every pair of individuals
    nj_row_names = []
    for nj_id in nj_id_pop_dict.keys(): 
        nj_row_names.append(f'{nj_id_pop_dict[nj_id]}__{nj_id}')
    nj_ajsub_df = pd.DataFrame(columns=list(aj_geo_id_dict[aj_subpop]), index=nj_row_names)
    for aj_id in aj_geo_id_dict[aj_subpop]:
        nj_ajsub_df[aj_id]=0
    # Use this dict the store the following key:value pair
    ## f'{AJ_SUBPOP_NAME}_{OTHER_SUBPOP_NAME}'(key): mean pairwise IBD length(float)
    ajsubpop_njpop_mean_pws_ibd_len= {}
    aj_pop_str = aj_subpop[1] if np.nan not in aj_subpop[:2] else 'unknown'
    ## key should be $nj_pop, value should be a dataframe of selected IBd segments
    nj_pop_ibd_seg_df_dict = {
        nj_pop: pd.DataFrame(columns=[str(i) for i in range(12)]) for nj_pop in nj_pop_id_dict.keys()
    }
    for chr in range(1, 23):
        # print(chr)
        # Calculate pairwise IBD length and average lod score of IBD segments for 
        ## each non-mixed & unknown AJ sub-groups (14 in total) against other non-AJ populations
        ### Sum up IBD segment lengths chromosome-wise and obtain the mean value finally after loop
        ibd_file = f'imputed935Jews_1320nonJews_reformatted/imputed935Jews_1320nonJews_L_m-200_chr{chr}.csv'
        if not path.isfile(ibd_file): continue
        df = pd.read_csv(ibd_file, header=0)
        df = df.loc[(df['id1'].isin({el: '' for el in aj_geo_id_dict[aj_subpop]}))|\
                   (df['id2'].isin({el: '' for el in aj_geo_id_dict[aj_subpop]}))]
        if chr == 1: 
            for ibd_seg_df in nj_pop_ibd_seg_df_dict.values(): 
                ibd_seg_df.columns = list(df.columns)
        # obtain the number of individuals from this AJ subgroup
        aj_subpop_size = len(aj_geo_id_dict[aj_subpop])
        ##########################################################################
        ##########################################################################
        ## Against NJs from 76 populations
        for nj_pop in nj_pop_id_dict.keys():
            ## Create an empty dataframe
            filtered_df = pd.DataFrame(columns=list(df.columns))
            #########################################################
            # Filter into the corresponding IBD segments
            for iid in list(nj_pop_id_dict[nj_pop]):
                # 'id' of NJ consists of fid and iid, like f"{fid}_{iid}"
                filtered_df_sub = df.loc[(df['id1'].str.endswith(f'_{iid}'))|(df['id2'].str.endswith(f'_{iid}'))]
                filtered_df = pd.concat([filtered_df, filtered_df_sub])
            nj_pop_ibd_seg_df_dict[nj_pop] = pd.concat([nj_pop_ibd_seg_df_dict[nj_pop], filtered_df])
            nj_pop_size = len(nj_pop_id_dict[nj_pop])
            pop_pair_str = f'{aj_pop_str}AJ_{nj_pop}'
            num_possible_pairs = aj_subpop_size * nj_pop_size
            if aj_subpop_size > 0 and nj_pop_size > 0:
                pairwise_avg_ibd_len = round(filtered_df['cm_len'].sum()/num_possible_pairs, 3)
                ## Create a string to describe the AJ subgroup-NJ group pair
                if pop_pair_str not in ajsubpop_njpop_mean_pws_ibd_len:
                    ## value of the dict is initially a list, [sum_of_ibd_lengths, subpop_size1, subpop_size2]
                    ### The 4th element records IBD length of each pair between 2 pops
                    ajsubpop_njpop_mean_pws_ibd_len[pop_pair_str]=[0, aj_subpop_size, nj_pop_size, {}]
                ajsubpop_njpop_mean_pws_ibd_len[pop_pair_str][0]+=filtered_df['cm_len'].sum()
            ##################################################
            # Output the CSV file storing IBD length between each pair of individuals
            filtered_row_num, col_num = filtered_df.shape
            # Loop through every found individual pair
            for i in range(filtered_row_num):
                ibd_row = filtered_df.iloc[i]
                # id is in the form of f'{fid}_{iid}', $fid is a digit
                id1=ibd_row['id1']
                iid1 = "_".join(tuple(id1.split("_")[1:]))
                id2=ibd_row['id2']
                iid2 = "_".join(tuple(id2.split("_")[1:]))
                ibd_len = ibd_row['cm_len']
                if id1 in aj_id_geo_dict and iid2 in nj_id_pop_dict:
                    nj_ajsub_df.loc[f'{nj_id_pop_dict[iid2]}__{iid2}', id1]+=round(float(ibd_len), 5)
                    indiv_pair_key = f'{id1}_{id2}'
                    if not indiv_pair_key in ajsubpop_njpop_mean_pws_ibd_len[pop_pair_str][3]:
                        ajsubpop_njpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]=0
                    ajsubpop_njpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]+=round(float(ibd_len), 5)
                elif iid1 in nj_id_pop_dict and id2 in aj_id_geo_dict:
                    nj_ajsub_df.loc[f'{nj_id_pop_dict[iid1]}__{iid1}', id2]+=round(float(ibd_len), 5)
                    indiv_pair_key = f'{id2}_{id1}'
                    if not indiv_pair_key in ajsubpop_njpop_mean_pws_ibd_len[pop_pair_str][3]:
                        ajsubpop_njpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]=0
                    ajsubpop_njpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]+=round(float(ibd_len), 5)
    ##########################################################################################
    # Finally, output IBD segments between AJ pop and each NJ pop to 76 CSV files
    for nj_pop, ibd_seg_df in nj_pop_ibd_seg_df_dict.items():
        ibd_seg_df.to_csv(f"IBD_segments_AJ_subpop_vs_76_nonJ_pop/{aj_subpop[1]}AJ__{nj_pop}_IBD_segs.csv", 
                           encoding='utf-8', index=False)
    # Calculate the mean pairwise IBD length for AJ versus NJ pairs
    ## after looping through 22 autosomes
    ### This dict should have 14*76 items
    with open("AJ_subpop_vs_76_nonJ_pop_mean_pairwise_ibd_length.csv", "a+") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        for k, v in ajsubpop_njpop_mean_pws_ibd_len.items():
            if type(v) is list:
                # print(v)
                if len(v)==4: 
                    # Finally, convert ajsubpop_njpop_mean_pws_ibd_len[pop_pair_str][3] to a list of values
                    ## And extend a list of zeros according to size of keys in this dict 
                    ### for calculating standard deviation
                    num_non_zero_ibd_indiv_pairs = len(v[3].keys())
                    num_possible_pairs = v[1]*v[2]
                    num_zeros_to_add = num_possible_pairs - num_non_zero_ibd_indiv_pairs
                    v[3]=list(v[3].values())
                    if num_zeros_to_add > 0:
                        v[3].extend([0]*num_zeros_to_add)
                    pop1=k.split("_")[0]
                    pop2="_".join(tuple(k.split("_")[1:]))
                    mean_v=round(v[0]/num_possible_pairs, 3)
                    # v[3] is a list
                    ibd_csv_writer.writerow([pop1,pop2,str(mean_v),str(v[1]),str(v[2]), round(st.sem(v[3]), 3)])
    ##########################################################################################
    # Additionally, output 14*2 .csv files storing total IBD length between every pairs of individuals
    ## after looping through 22 autosomes
    ### Reference for conversion from df to csv: 
    ####: https://towardsdatascience.com/how-to-export-pandas-dataframe-to-csv-2038e43d9c03
    nj_ajsub_df.to_csv(f"{aj_subpop[1]}AJ_nonJew.csv", encoding='utf-8')
    

if __name__ == '__main__':
    #warnings.filterwarnings("ignore")
    # load 901 Jews into variables
    df_jews = pd.read_csv('Bray2010_Behar1013_Gladstein2019_Kopelman2020_901Jews.csv', header=0)
    df_aj = df_jews.loc[lambda x: x['data_source'].str.contains('AJ$', regex = True)]
    ## Create two dictionaries according to unique geographic origin
    df_aj_geo_groups= df_aj.groupby(['city1', 'country1', 'city2', 'country2']).groups
    aj_geo_id_dict, aj_id_geo_dict= {}, {}
    
    ## non-mixed & unknown AJs, city2 & country2 must be nan
    for key, val in df_aj_geo_groups.items():
        for ind in val:
            fid=df_jews['family_id'].iloc[ind]
            iid=df_jews['individual_id'].iloc[ind]
            cat_id = "_".join((fid, iid))
            ## If non nan are present in key(a tuple of 4 elements) it's a mixed indiv, pass
            if not np.nan in key: 
                pass
            ## only NaN in key, thus being geographically unknown
            #elif len(set(key)) == 1:
                #pass
            elif not "Italy" in key and not "France" in key:
                if key not in aj_geo_id_dict: aj_geo_id_dict[key]=[]  
                aj_geo_id_dict[key].append(cat_id)
                aj_id_geo_dict[cat_id]=key
    
    ## Non-Jews from 76 extra ethnic groups
    nj_pop_id_dict, nj_id_pop_dict = {}, {}
    with open("76pop_non_jews_id_file_list.txt", "r") as f:
        for line in f:
            pop_name=line.strip().split("/")[2].split(".")[0]
            if pop_name not in nj_pop_id_dict: nj_pop_id_dict[pop_name]=[]
            with open(line.strip(), "r") as id_f:
                for line in id_f:
                    iid=line.strip()
                    nj_pop_id_dict[pop_name].append(iid)
                    if iid not in nj_id_pop_dict: nj_id_pop_dict[iid]=pop_name
    
    ## Output .csv file for mean pairwise IBD length
    with open("AJ_subpop_vs_76_nonJ_pop_mean_pairwise_ibd_length.csv", "w") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ibd_csv_writer.writerow(['pop1', 'pop2', 'mean_pairwise_IBD_length', 'pop1_size', 'pop2_size', \
                                 'IBD_len_SE'])
    
    # For shorter computing time, 
    ## Run multi-processing using non-mixed & unknown subpop groups in AJs
    partial_func = partial(write_aj_otherpop_pair_ibd, aj_geo_id_dict=aj_geo_id_dict, \
        nj_pop_id_dict=nj_pop_id_dict, aj_id_geo_dict=aj_id_geo_dict, nj_id_pop_dict=nj_id_pop_dict)
    # the pool occupies 14 threads(14 sub-populations of AJ in parallel, excluding the Italian & French "AJ")
    pool = Pool(14)
    pool.map(partial_func, list(aj_geo_id_dict.keys()))
    pool.close()
    pool.join()
```


#### Extracting IBD between AJs and other Jews
```bash
%%bash
mkdir -p IBD_segments_AJ_subpop_vs_23_otherJ_pop
```

```python
from multiprocessing import Pool
import pandas as pd
import csv 
import re
import numpy as np
from os import path
from functools import partial
import scipy.stats as st
import warnings

'''
There're 901 modern Jewish individuals in total, including 508 Ashkenazic Jews(AJ) and 393 other Jews(OJ).
This script presents average IBD segment length per pair in Centimorgans.
'''

def write_aj_otherpop_pair_ibd(aj_subpop: tuple, aj_geo_id_dict: dict, oj_geo_id_dict: dict, \
                               aj_id_geo_dict: dict, oj_id_geo_dict: dict):
    # One CSV file as output from a pandas.dataframe per AJ subpop (14*2 CSV files in total)
    ## Two dataframes to store total IBD length between every pair of individuals
    aj_subpop_size = len(aj_geo_id_dict[aj_subpop])
    oj_row_names = []
    for oj_id in oj_id_geo_dict.keys(): 
        oj_col_prefix = f'{oj_id_geo_dict[oj_id][0]}-{oj_id_geo_dict[oj_id][1]}'
        oj_row_names.append(f'{oj_col_prefix}__{oj_id}')
    oj_ajsub_df = pd.DataFrame(columns=list(aj_geo_id_dict[aj_subpop]), index=oj_row_names)
    for aj_id in aj_geo_id_dict[aj_subpop]:
        oj_ajsub_df[aj_id]=0
    # Use this dict the store the following key:value pair
    ## f'{AJ_SUBPOP_NAME}_{OTHER_SUBPOP_NAME}'(key): mean pairwise IBD length(float)
    ajsubpop_ojsubpop_mean_pws_ibd_len = {}
    aj_pop_str = aj_subpop[1] if np.nan not in aj_subpop[:2] else 'unknown'
    oj_pop_ibd_seg_df_dict = {
        f'{oj_pop[0]}-{oj_pop[1]}OJ': 
            pd.DataFrame(columns=[str(i) for i in range(12)]) for oj_pop in oj_geo_id_dict.keys()
    }
    for chr in range(1, 23):
        # Calculate pairwise IBD length and average lod score of IBD segments for 
        ## each non-mixed & unknown AJ sub-groups (14 in total, including unknown ones) against other non-AJ populations
        ### Sum up IBD segment lengths chromosome-wise and obtain the mean value finally after loop
        ibd_file = f'imputed935Jews_1320nonJews_reformatted/imputed935Jews_1320nonJews_L_m-200_chr{chr}.csv'
        if not path.isfile(ibd_file): continue
        df = pd.read_csv(ibd_file, header=0)
        df = df.loc[(df['id1'].isin({el: '' for el in aj_geo_id_dict[aj_subpop]}))|\
                   (df['id2'].isin({el: '' for el in aj_geo_id_dict[aj_subpop]}))]
        # obtain the number of individuals from this AJ subgroup
        ##########################################################################
        ##########################################################################
        ## Against 393 OJs  
        for oj_pop in oj_geo_id_dict.keys():
            #########################################################
            # Filter into the corresponding IBD segments
            filtered_df = df.loc[(df['id1'].isin({el: '' for el in oj_geo_id_dict[oj_pop]}))|\
              (df['id2'].isin({el: '' for el in oj_geo_id_dict[oj_pop]}))]
            if chr == 1: 
                for ibd_seg_df in oj_pop_ibd_seg_df_dict.values(): 
                    ibd_seg_df.columns = list(filtered_df.columns)
            oj_pop_size = len(oj_geo_id_dict[oj_pop])
            oj_pop_str = f'{oj_pop[0]}-{oj_pop[1]}' if np.nan not in oj_pop else 'unknown'
            oj_pop_ibd_seg_df_dict[f'{oj_pop_str}OJ'] = \
                pd.concat([oj_pop_ibd_seg_df_dict[f'{oj_pop_str}OJ'], filtered_df])
            pop_pair_str = f'{aj_pop_str}AJ_{oj_pop_str}OJ'
            if aj_subpop_size > 0 and oj_pop_size > 0:
                ## Create a string to describe the AJ subgroup-OJ subgroup pair
                if pop_pair_str not in ajsubpop_ojsubpop_mean_pws_ibd_len:
                    ## value of the dict is initially a list, [sum_of_ibd_lengths, subpop_size1, subpop_size2]
                    ajsubpop_ojsubpop_mean_pws_ibd_len[pop_pair_str]=[0, aj_subpop_size, oj_pop_size, {}]
                ajsubpop_ojsubpop_mean_pws_ibd_len[pop_pair_str][0]+=filtered_df['cm_len'].sum()
            ##################################################
            # Output the CSV file storing IBD length between each pair of individuals
            filtered_row_num, col_num = filtered_df.shape
            for i in range(filtered_row_num):
                ibd_row = filtered_df.iloc[i]
                # id is in the form of f'{fid}_{iid}', $fid is a digit
                id1=ibd_row['id1']
                id2=ibd_row['id2']
                ibd_len = ibd_row['cm_len']
                ###Test
                #print(ibd_len)
                ###
                if id1 in aj_id_geo_dict and id2 in oj_id_geo_dict:
                    id2_prefix = f'{oj_id_geo_dict[id2][0]}-{oj_id_geo_dict[id2][1]}'
                    oj_ajsub_df.loc[f'{id2_prefix}__{id2}', id1]+=round(float(ibd_len), 5)
                    ###Test
                    #print(oj_ajsub_df.loc[f'{id2_prefix}_{id2}', id1])
                    ###
                    indiv_pair_key = f'{id1}_{id2}'
                    if not indiv_pair_key in ajsubpop_ojsubpop_mean_pws_ibd_len[pop_pair_str][3]:
                        ajsubpop_ojsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]=0
                    ajsubpop_ojsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]+=round(float(ibd_len), 5)
                elif id1 in oj_id_geo_dict and id2 in aj_id_geo_dict:
                    id1_prefix = f'{oj_id_geo_dict[id1][0]}-{oj_id_geo_dict[id1][1]}'
                    oj_ajsub_df.loc[f'{id1_prefix}__{id1}', id2]+=round(float(ibd_len), 5)
                    #print(oj_ajsub_df.loc[f'{id1_prefix}_{id1}', id2])
                    indiv_pair_key = f'{id2}_{id1}'
                    if not indiv_pair_key in ajsubpop_ojsubpop_mean_pws_ibd_len[pop_pair_str][3]:
                        ajsubpop_ojsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]=0
                    ajsubpop_ojsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]+=round(float(ibd_len), 5)
                    ###Test
                    #print(ajsubpop_ojsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key])
                    ###
    ##########################################################################################
    ##########################################################################################
    for oj_pop, ibd_seg_df in oj_pop_ibd_seg_df_dict.items():
        ibd_seg_df.to_csv(f"IBD_segments_AJ_subpop_vs_23_otherJ_pop/{aj_subpop[1]}AJ__{oj_pop}_IBD_segs.csv", 
                           encoding='utf-8', index=False)
    # Finally, calculate the mean pairwise IBD length for AJ versus OJ and AJ versus NJ pairs
    ## after looping through 22 autosomes
    ### This dict should have 14*23 items
    with open("AJ_subpop_vs_otherJew_23_subpop_mean_pairwise_ibd_length.csv", "a+") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ### This dict should have 14*22 items
        for k, v in ajsubpop_ojsubpop_mean_pws_ibd_len.items():
            if type(v) is list:
                if len(v)==4:
                    # Finally, convert ajsubpop_ojsubpop_mean_pws_ibd_len[pop_pair_str][3] to a list of values
                    ## And extend a list of zeros according to size of keys in this dict 
                    ### for calculating standard deviation
                    #print(v[3])
                    num_non_zero_ibd_indiv_pairs = len(v[3].keys())
                    num_possible_pairs = int(v[1])*int(v[2])
                    num_zeros_to_add = num_possible_pairs - num_non_zero_ibd_indiv_pairs
                    v[3] = list(v[3].values())
                    if num_zeros_to_add > 0:
                        v[3].extend([0]*num_zeros_to_add)
                    pop1=k.split("_")[0]
                    pop2="_".join(tuple(k.split("_")[1:]))
                    mean_v = round(v[0]/num_possible_pairs, 3)
                    # v[3] is a list
                    #print(v[3])
                    ibd_csv_writer.writerow([pop1,pop2,str(mean_v),str(v[1]),str(v[2]),round(st.sem(v[3]), 3)])
    ##########################################################################################
    # Additionally, output 14*2 .csv files storing total IBD length between every pairs of individuals
    ## after looping through 22 autosomes
    ### Reference for conversion from df to csv: 
    ####: https://towardsdatascience.com/how-to-export-pandas-dataframe-to-csv-2038e43d9c03
    oj_ajsub_df.to_csv(f"{aj_subpop[1]}AJ_otherJew.csv", encoding='utf-8')
    

if __name__ == '__main__':
    # warnings.filterwarnings("ignore")
    # load 901 Jews into variables
    df_jews = pd.read_csv('Bray2010_Behar1013_Gladstein2019_Kopelman2020_901Jews.csv', header=0)
    df_aj = df_jews.loc[lambda x: x['data_source'].str.contains('AJ$', regex = True)]
    df_otherj = df_jews.loc[lambda x: x['data_source'].str.contains('otherJ$', regex = True)]
    ## add single-origin French Jews to df_otherJ
    df_FrenchJ = df_jews.loc[(df_jews['country1']=="France") & (df_jews['country2'].isnull())]
    ## add single-origin Italian Jews to df_otherJ
    df_ItalianJ = df_jews.loc[df_jews['country1']=="Italy"]
    ## Concat three dfs and drop duplicates
    df_otherj = pd.concat([df_otherj, df_FrenchJ, df_ItalianJ]).drop_duplicates()
    ## Create two dictionaries according to unique geographic origin
    df_aj_geo_groups= df_aj.groupby(['city1', 'country1', 'city2', 'country2']).groups
    df_oj_geo_groups= df_otherj.groupby(['city1', 'country1']).groups
    ##########################
    # print(df_oj_geo_groups)
    
    aj_geo_id_dict, oj_geo_id_dict = {}, {}
    aj_id_geo_dict, oj_id_geo_dict = {}, {}
    
    ## non-mixed & unknown AJs, city2 & country2 must be nan
    for key, val in df_aj_geo_groups.items():
        for ind in val:
            fid=df_jews['family_id'].iloc[ind]
            iid=df_jews['individual_id'].iloc[ind]
            cat_id = "_".join((fid, iid))
            ## If non nan are present in key(a tuple of 4 elements) it's a mixed indiv, pass
            if not np.nan in key: 
                pass
            ## only NaN in key, thus being geographically unknown
            #elif len(set(key)) == 1:
                #pass
            elif not "Italy" in key and not "France" in key:
                if key not in aj_geo_id_dict: aj_geo_id_dict[key]=[]  
                aj_geo_id_dict[key].append(cat_id)
                aj_id_geo_dict[cat_id]=key
    
    ## Other-Jews
    for key, val in df_oj_geo_groups.items():
        for ind in val:
            ## avoid mixed individuals
            if not np.nan in key: 
                pass
            fid=df_jews['family_id'].iloc[ind]
            iid=df_jews['individual_id'].iloc[ind]
            if key not in oj_geo_id_dict: oj_geo_id_dict[key]=[]
            cat_id = "_".join((fid, iid))
            oj_geo_id_dict[key].append(cat_id)
            oj_id_geo_dict[cat_id]=key
    
    ## Output .csv file for mean pairwise IBD length
    with open("AJ_subpop_vs_otherJew_23_subpop_mean_pairwise_ibd_length.csv", "w") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ibd_csv_writer.writerow(['pop1', 'pop2', 'mean_pairwise_IBD_length', 'pop1_size', 'pop2_size', \
                                 'IBD_len_SE'])
    
    # For shorter computing time, 
    ## Run multi-processing using non-mixed & unknown subpop groups in AJs
    partial_func = partial(write_aj_otherpop_pair_ibd, aj_geo_id_dict=aj_geo_id_dict, \
        oj_geo_id_dict=oj_geo_id_dict, aj_id_geo_dict=aj_id_geo_dict, oj_id_geo_dict=oj_id_geo_dict)
    # the pool occupies 14 threads(14 sub-populations of AJ in parallel)
    pool = Pool(14)
    pool.map(partial_func, list(aj_geo_id_dict.keys()))
    pool.close()
    pool.join()
```

#### Extracting IBD within each Jewish community
```bash
%%bash
mkdir -p Jewish_pop_intra_IBD -p IBD_segments_intra_Jewish_pops
```


```python
from multiprocessing import Pool
import pandas as pd
import csv 
import numpy as np
from os import path
from functools import partial
import scipy.stats as st
import warnings
import re

'''
There're 901 modern Jewish individuals in total, including 508 Ashkenazic Jews(AJ) and 393 other Jews(OJ).
This script presents average IBD segment length per pair in Centimorgans.
'''

def write_jew_intra_pop_pair_ibd(jew_pop: tuple, jew_geo_id_dict: dict):
    ## jew_pop is a tuple of two elements. 
    ### The first element is the real name of Jewish population. 
    ### The second element is its key at jew_geo_id_dict
    jew_pop_indivs = jew_geo_id_dict[jew_pop[1]]
    ## obtain the number of individuals from this Jewish population
    jew_pop_size = len(jew_pop_indivs)
    ## If n<3, we should not calculate intra-group IBD 
    if jew_pop_size < 2: return
    ## exclude IBD between two haplotypes of the same individual
    num_indiv_pairs = jew_pop_size*(jew_pop_size-1)/2
    jew_indiv_pair_ibd_dict = {}
    jew_pop_intra_ibd_seg_df = pd.DataFrame(columns=[str(i) for i in range(12)])
    for i in range(jew_pop_size):
        for j in range(i+1, jew_pop_size):
            iid1, iid2 = jew_pop_indivs[i], jew_pop_indivs[j]
            jew_indiv_pair_ibd_dict[f'{iid1}__{iid2}'] = 0
    ## Create a dataframe to store sums of intra-population individual pairwise IBD length
    ### Output the dataframe into a csv file later
    jew_pop_intra_df = pd.DataFrame(columns=['IBD_cM'], index=list(jew_indiv_pair_ibd_dict.keys()))
    jew_pop_intra_df['IBD_cM'] = 0
    #print(jew_pop_intra_df.index)
    ## Loop through 22 autosomes
    for chr in range(1, 23):
        ### Sum up IBD segment lengths chromosome-wise and obtain the mean value finally after loop
        ibd_file = f'imputed935Jews_1320nonJews_reformatted/imputed935Jews_1320nonJews_L_m-200_chr{chr}.csv'
        if not path.isfile(ibd_file): continue
        df = pd.read_csv(ibd_file, header=0)
        pop_intra_df = df.loc[\
                (df['id1'].isin({iid: "" for iid in jew_pop_indivs})) & 
                (df['id2'].isin({iid: "" for iid in jew_pop_indivs})) 
            ]
        if chr == 1: 
            jew_pop_intra_ibd_seg_df.columns = list(pop_intra_df.columns)
        jew_pop_intra_ibd_seg_df = pd.concat([jew_pop_intra_ibd_seg_df, pop_intra_df])
        ibd_row_num = pop_intra_df.shape[0]
        for i in range(ibd_row_num):
            ibd_row = pop_intra_df.iloc[i]
            # id is in the form of f'{fid}_{iid}', $fid is a digit
            iid1=ibd_row['id1']
            iid2=ibd_row['id2']
            ibd_len = ibd_row['cm_len']
            #print(iid1, iid2)
            if f'{iid1}__{iid2}' in jew_indiv_pair_ibd_dict:
                #print("matched")
                jew_pop_intra_df.loc[f'{iid1}__{iid2}', 'IBD_cM']+=round(float(ibd_len), 5)
            elif f'{iid2}__{iid1}' in jew_indiv_pair_ibd_dict:
                #print("matched")
                jew_pop_intra_df.loc[f'{iid2}__{iid1}', 'IBD_cM']+=round(float(ibd_len), 5)
    ## Write the dataframe of selected IBD segments to CSV file
    jew_pop_intra_ibd_seg_df.to_csv(f"IBD_segments_intra_Jewish_pops/{jew_pop[0]}_IBD_segs.csv", 
                           encoding='utf-8', index=False)
    ## Write to the general CSV file
    with open("all_Jewish_subpops_intra_mean_pairwise_ibd_length.csv", "a") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        mean_v = round(jew_pop_intra_df['IBD_cM'].sum()/num_indiv_pairs, 3)
        ibd_csv_writer.writerow([jew_pop[0], str(mean_v), jew_pop_size, \
                                 round(st.sem(jew_pop_intra_df['IBD_cM'].tolist()), 3)])
    ## Write to an individual CSV file
    ## "jew_pop" is a tuple in the form like ("Berlin", "Germany", nan, nan), 
    ### so the second element annotates country
    jew_pop_intra_df.round(decimals=4).to_csv(f"Jewish_pop_intra_IBD/{jew_pop[0]}_intraIBD.csv", encoding='utf-8')
    

if __name__ == '__main__':
    # warnings.filterwarnings("ignore")
    # load 901 Jews into variables
    df_jews = pd.read_csv('Bray2010_Behar1013_Gladstein2019_Kopelman2020_901Jews.csv', header=0)
    df_aj = df_jews.loc[lambda x: x['data_source'].str.contains('AJ$', regex = True)]
    df_otherj = df_jews.loc[lambda x: x['data_source'].str.contains('otherJ$', regex = True)]
    ## add single-origin French Jews to df_otherJ
    df_FrenchJ = df_jews.loc[(df_jews['country1']=="France") & (df_jews['country2'].isnull())]
    ## add single-origin Italian Jews to df_otherJ
    df_ItalianJ = df_jews.loc[df_jews['country1']=="Italy"]
    ## Concat three dfs and drop duplicates
    df_otherj = pd.concat([df_otherj, df_FrenchJ, df_ItalianJ]).drop_duplicates()
    ## Create two dictionaries according to unique geographic origin
    df_aj_geo_groups= df_aj.groupby(['city1', 'country1', 'city2', 'country2']).groups
    df_oj_geo_groups= df_otherj.groupby(['city1', 'country1']).groups
    ##########################
    # print(df_oj_geo_groups)
    
    aj_geo_id_dict, oj_geo_id_dict = {}, {}
    
    ## non-mixed & unknown AJs, city2 & country2 must be nan
    for key, val in df_aj_geo_groups.items():
        for ind in val:
            fid=df_jews['family_id'].iloc[ind]
            iid=df_jews['individual_id'].iloc[ind]
            cat_id = "_".join((fid, iid))
            ## If non nan are present in key(a tuple of 4 elements) it's a mixed indiv, pass
            if not np.nan in key: 
                pass
            ## only NaN in key, thus being geographically unknown
            #elif len(set(key)) == 1:
                #pass
            elif not "Italy" in key and not "France" in key:
                if key not in aj_geo_id_dict: aj_geo_id_dict[key]=[]  
                aj_geo_id_dict[key].append(cat_id)
    
    ## Other-Jews
    for key, val in df_oj_geo_groups.items():
        for ind in val:
            ## avoid mixed individuals
            if not np.nan in key: 
                pass
            fid=df_jews['family_id'].iloc[ind]
            iid=df_jews['individual_id'].iloc[ind]
            if key not in oj_geo_id_dict: oj_geo_id_dict[key]=[]
            cat_id = "_".join((fid, iid))
            oj_geo_id_dict[key].append(cat_id)
    
    ## use this variable to store population names for all Jews
    jewish_pops = []
    for key in aj_geo_id_dict:
        pop_el = (f'{key[1]}AJ', key)
        jewish_pops.append(pop_el)
    for key in oj_geo_id_dict:
        if key[0] == "Mosul":
            pop_el = ('KurdistanOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Ar Raqqah':
            pop_el = ('SyrianKurdistanOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Mumbai':
            pop_el = ('IndiaMumbaiOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Thiruvananthapuram':
            pop_el = ('IndiaCochinOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Tbilisi':
            pop_el = ('GeorgiaOJ', key)
            jewish_pops.append(pop_el)
        else:
            pop_el = (f'{key[1]}OJ', key)
            jewish_pops.append(pop_el)
    
    #print(jewish_pops)
    aj_geo_id_dict.update(oj_geo_id_dict)
    jews_geo_dict = aj_geo_id_dict
    aj_geo_id_dict = None
    #print(jews_geo_dict)
    ## Output .csv file for mean pairwise IBD length
    with open("all_Jewish_subpops_intra_mean_pairwise_ibd_length.csv", "w") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ibd_csv_writer.writerow(['pop', 'mean_pairwise_IBD_length', 'pop_size', 'IBD_len_SE'])
    
    ## For shorter computing time, 
    ## Run multi-processing using non-mixed subpop AJs and OJs
    partial_func = partial(write_jew_intra_pop_pair_ibd, jew_geo_id_dict=jews_geo_dict)
    ## the pool occupies 14 threads(14 sub-populations of AJ in parallel, excluding the Italian & French "AJ")
    ## There should be n*(n-1)/2 pairs within each group to compare IBD
    ## However, only one Moldavian AJ and one Slovakian AJ are present. 
    ### Therefore, it does not make sense to investigate these two groups
    pool = Pool(20)
    pool.map(partial_func, jewish_pops)
    pool.close()
    pool.join()
```

### In addition, let's write two scripts to extract IBD sharing between different AJ sub-populations and between different OJ sub-populations.

```bash
%%bash
mkdir -p aj_pops_inter_IBD -p oj_pops_inter_IBD -p IBD_segments_inter_aj_pops -p IBD_segments_inter_oj_pops
```


```python
from multiprocessing import Pool
import pandas as pd
import csv 
import numpy as np
from os import path
from functools import partial
import scipy.stats as st
import warnings
import re

'''
There're 901 modern Jewish individuals in total, including 508 Ashkenazic Jews(AJ) and 393 other Jews(other_aj).
This script presents average IBD segment length per pair in Centimorgans.
'''

def write_aj_aj_pop_pair_ibd(aj_pop: tuple, aj_geo_id_dict: dict, aj_id_geo_dict: dict):
    ## aj_pop is a tuple of two elements. 
    ### The first element is the real name of Jewish population. 
    ### The second element is its key at aj_geo_id_dict
    aj_pop_indivs = aj_geo_id_dict[aj_pop[1]]
    aj_pop_indivs_dict = {iid: "" for iid in aj_pop_indivs}
    ## obtain the number of individuals from this Jewish population
    aj_pop_size = len(aj_pop_indivs)
    other_aj_row_names = []
    for other_aj_id in aj_id_geo_dict.keys(): 
        ## the prefix should be country name = "AJ"
        if aj_id_geo_dict[other_aj_id] == aj_pop[1]: continue
        other_aj_col_prefix = f'{aj_id_geo_dict[other_aj_id][1]}AJ'
        other_aj_row_names.append(f'{other_aj_col_prefix}__{other_aj_id}')
    aj_other_aj_df = pd.DataFrame(columns=aj_pop_indivs, index=other_aj_row_names)
    for aj_id in aj_pop_indivs:
        aj_other_aj_df[aj_id]=0
    # Use this dict the store the following key:value pair
    ## f'{AJ_SUBPOP_NAME}_{OTHER_SUBPOP_NAME}'(key): mean pairwise IBD length(float)
    ajsubpop_other_ajsubpop_mean_pws_ibd_len = {}
    other_aj_pop_ibd_seg_df_dict = {
        other_aj_pop[1]: pd.DataFrame(columns=[str(i) for i in range(12)]) 
        for other_aj_pop in aj_geo_id_dict.keys() if not other_aj_pop == aj_pop[1] 
    }
    ## Loop through 22 autosomes
    for chr in range(1, 23):
        ### Sum up IBD segment lengths chromosome-wise and obtain the mean value finally after loop
        ibd_file = f'imputed935Jews_1320nonJews_reformatted/imputed935Jews_1320nonJews_L_m-200_chr{chr}.csv'
        if not path.isfile(ibd_file): continue
        df = pd.read_csv(ibd_file, header=0)
        aj_subpop_ibd_df = df.loc[
            (df['id1'].isin(aj_pop_indivs_dict)) | (df['id2'].isin(aj_pop_indivs_dict)) 
        ]
        if chr == 1: 
            for ibd_seg_df in other_aj_pop_ibd_seg_df_dict.values(): 
                ibd_seg_df.columns = list(aj_subpop_ibd_df.columns)
        for other_aj_pop in aj_geo_id_dict.keys():
            #########################################################
            # Filter into the corresponding IBD segments
            if other_aj_pop == aj_pop[1]: continue
            other_aj_subpop_id_dict = {el: '' for el in aj_geo_id_dict[other_aj_pop]}
            filtered_df = aj_subpop_ibd_df.loc[
                (df['id1'].isin(other_aj_subpop_id_dict)) | (df['id2'].isin(other_aj_subpop_id_dict))
            ]
            other_aj_pop_size = len(aj_geo_id_dict[other_aj_pop])
            other_aj_pop_str = other_aj_pop[1] if np.nan not in other_aj_pop else 'unknown'
            other_aj_pop_ibd_seg_df_dict[other_aj_pop_str] = \
            pd.concat([other_aj_pop_ibd_seg_df_dict[other_aj_pop_str], filtered_df])
            pop_pair_str = f'{aj_pop[0]}_{other_aj_pop_str}AJ'
            ## use the condition "aj_pop[0] > other_aj_pop_str" to avoid replicate pairs, i.e. A & B and B & A
            if aj_pop_size > 0 and other_aj_pop_size > 0 and aj_pop[0] > other_aj_pop_str:
                ## Create a string to describe the AJ subgroup-OJ subgroup pair
                if pop_pair_str not in ajsubpop_other_ajsubpop_mean_pws_ibd_len:
                    ## value of the dict is initially a list, [sum_of_ibd_lengths, subpop_size1, subpop_size2]
                    ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str]=[0, aj_pop_size, other_aj_pop_size, {}]
                ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str][0]+=filtered_df['cm_len'].sum()
            ##################################################
            # Output the CSV file storing IBD length between each pair of individuals
            filtered_row_num = filtered_df.shape[0]
            #print(chr, pop_pair_str, filtered_row_num)
            for i in range(filtered_row_num):
                ibd_row = filtered_df.iloc[i]
                # id is in the form of f'{fid}_{iid}', $fid is a digit
                id1=ibd_row['id1']
                id2=ibd_row['id2']
                ibd_len = ibd_row['cm_len']
                ###Test
                #print(ibd_len)
                ###
                if id1 in aj_pop_indivs_dict and id2 in other_aj_subpop_id_dict:
                    ## the prefix should be country name + 'AJ'
                    id2_prefix = f'{aj_id_geo_dict[id2][1]}AJ'
                    aj_other_aj_df.loc[f'{id2_prefix}__{id2}', id1]+=round(float(ibd_len), 5)
                    ###Test
                    #print(aj_other_aj_df.loc[f'{id2_prefix}__{id2}', id1])
                    ###
                    indiv_pair_key = f'{id1}__{id2}'
                elif id2 in aj_pop_indivs_dict and id1 in other_aj_subpop_id_dict:
                    id1_prefix = f'{aj_id_geo_dict[id1][1]}AJ'
                    aj_other_aj_df.loc[f'{id1_prefix}__{id1}', id2]+=round(float(ibd_len), 5)
                    indiv_pair_key = f'{id2}__{id1}'
                if pop_pair_str in ajsubpop_other_ajsubpop_mean_pws_ibd_len:
                    #print(pop_pair_str)
                    #print(ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str][3].keys())
                    if not indiv_pair_key in ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str][3]:
                        ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]=0
                    ajsubpop_other_ajsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]+=\
                        round(float(ibd_len), 5)
    ##########################################################################################
    ##########################################################################################
    for other_aj_pop, ibd_seg_df in other_aj_pop_ibd_seg_df_dict.items():
        ibd_seg_df.to_csv(f"IBD_segments_inter_aj_pops/{aj_pop[0]}__{other_aj_pop}AJ_IBD_segs.csv", 
                           encoding='utf-8', index=False)
    # Finally, calculate the mean pairwise IBD length for AJ versus other AJs
    ## after looping through 22 autosomes
    ### This dict should have 14*48 items
    with open("AJ_subpop_inter_mean_pairwise_ibd_length.csv", "a+") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ### This dict should have 14*13/2 items
        for k, v in ajsubpop_other_ajsubpop_mean_pws_ibd_len.items():
            if type(v) is list:
                if len(v)==4:
                    num_non_zero_ibd_indiv_pairs = len(v[3].keys())
                    num_possible_pairs = int(v[1])*int(v[2])
                    num_zeros_to_add = num_possible_pairs - num_non_zero_ibd_indiv_pairs
                    v[3] = list(v[3].values())
                    if num_zeros_to_add > 0:
                        v[3].extend([0]*num_zeros_to_add)
                    pop1=k.split("_")[0]
                    pop2="_".join(tuple(k.split("_")[1:]))
                    mean_v = round(v[0]/num_possible_pairs, 3)
                    # v[3] is a list
                    #print(v[3])
                    ibd_csv_writer.writerow([pop1,pop2,str(mean_v),str(v[1]),str(v[2]),\
                                             round(st.sem(v[3]), 3)])
    #################################################################################################
    ## Additionally, output 14 .csv files storing total IBD length between every pairs of individuals
    ## after looping through 22 autosomes
    ### Reference for conversion from df to csv: 
    ####: https://towardsdatascience.com/how-to-export-pandas-dataframe-to-csv-2038e43d9c03
    #### Write to an individual CSV file
    #### "aj_pop" is a tuple in the form like ("Berlin", "Germany", nan, nan), 
    #### so the second element annotates country
    aj_other_aj_df.round(decimals=4).to_csv(f"aj_pops_inter_IBD/{aj_pop[0]}_IBD.csv", encoding='utf-8')
    

if __name__ == '__main__':
    # warnings.filterwarnings("ignore")
    # load 901 Jews into variables
    df_jews = pd.read_csv('Bray2010_Behar1013_Gladstein2019_Kopelman2020_901Jews.csv', header=0)
    df_aj = df_jews.loc[lambda x: x['data_source'].str.contains('AJ$', regex = True)]
    ## Create two dictionaries according to unique geographic origin
    df_aj_geo_groups= df_aj.groupby(['city1', 'country1', 'city2', 'country2']).groups
    aj_geo_id_dict, aj_id_geo_dict = {}, {}
    
    ## non-mixed & unknown AJs, city2 & country2 must be nan
    for key, val in df_aj_geo_groups.items():
        for ind in val:
            fid=df_jews['family_id'].iloc[ind]
            iid=df_jews['individual_id'].iloc[ind]
            cat_id = "_".join((fid, iid))
            ## If non nan are present in key(a tuple of 4 elements) it's a mixed indiv, pass
            if not np.nan in key: 
                pass
            ## only NaN in key, thus being geographically unknown
            #elif len(set(key)) == 1:
                #pass
            elif not "Italy" in key and not "France" in key:
                if key not in aj_geo_id_dict: aj_geo_id_dict[key]=[]  
                aj_geo_id_dict[key].append(cat_id)
                aj_id_geo_dict[cat_id]=key
    
    ## use this variable to store population names for all Jews
    jewish_pops = []
    for key in aj_geo_id_dict:
        pop_el = (f'{key[1]}AJ', key)
        jewish_pops.append(pop_el)
    
    #print(jewish_pops)
    ## Output .csv file for mean pairwise IBD length
    with open("AJ_subpop_inter_mean_pairwise_ibd_length.csv", "w") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ibd_csv_writer.writerow(['pop1', 'pop2', 'mean_pairwise_IBD_length', 'pop1_size', 'pop2_size', \
                                 'IBD_len_SE'])
    
    ## For shorter computing time, 
    ## Run multi-processing using non-mixed subpop AJs and other_ajs
    partial_func = partial(write_aj_aj_pop_pair_ibd, aj_geo_id_dict=aj_geo_id_dict, aj_id_geo_dict=aj_id_geo_dict)
    ## the pool occupies 14 threads(14 sub-populations of AJ in parallel, excluding the Italian & French "AJ")
    ## There should be 14*(14-1)/2 pairs of AJs to compare IBD
    ## However, only one Moldavian AJ and one Slovakian AJ are present. 
    pool = Pool(14)
    pool.map(partial_func, jewish_pops)
    pool.close()
    pool.join()
```

```python
from multiprocessing import Pool
import pandas as pd
import csv 
import numpy as np
from os import path
from functools import partial
import scipy.stats as st
import warnings
import re

'''
There're 901 modern Jewish individuals in total, including 508 Ashkenazic Jews(OJ) and 393 other Jews(OJ).
This script presents average IBD segment length per pair in Centimorgans.
'''

def write_oj_oj_pop_pair_ibd(oj_pop: tuple, oj_geo_id_dict: dict, oj_id_geo_dict: dict):
    ## oj_pop is a tuple of two elements. 
    ### The first element is the real name of Jewish population. 
    ### The second element is its key at oj_geo_id_dict
    oj_pop_indivs = oj_geo_id_dict[oj_pop[1]]
    oj_pop_indivs_dict = {iid: "" for iid in oj_pop_indivs}
    ## obtain the number of individuals from this Jewish population
    oj_pop_size = len(oj_pop_indivs)
    other_oj_row_names = []
    for other_oj_id in oj_id_geo_dict.keys(): 
        ## the prefix should be country name = "OJ"
        if oj_id_geo_dict[other_oj_id] == oj_pop[1]: continue
        other_oj_col_prefix = f'{oj_id_geo_dict[other_oj_id][0]}-{oj_id_geo_dict[other_oj_id][1]}OJ'
        other_oj_row_names.append(f'{other_oj_col_prefix}__{other_oj_id}')
    oj_other_oj_df = pd.DataFrame(columns=oj_pop_indivs, index=other_oj_row_names)
    for oj_id in oj_pop_indivs:
        oj_other_oj_df[oj_id]=0
    # Use this dict the store the following key:value pair
    ## f'{OJ_SUBPOP_NAME}_{OTHER_SUBPOP_NAME}'(key): mean pairwise IBD length(float)
    ojsubpop_other_ojsubpop_mean_pws_ibd_len = {}
    other_oj_pop_ibd_seg_df_dict = {
        f'{other_oj_pop[0]}-{other_oj_pop[1]}': pd.DataFrame(columns=[str(i) for i in range(12)]) 
        for other_oj_pop in oj_geo_id_dict.keys() if not other_oj_pop == oj_pop[1]
    }
    ## Loop through 22 autosomes
    for chr in range(1, 23):
        ### Sum up IBD segment lengths chromosome-wise and obtain the mean value finally after loop
        ibd_file = f'imputed935Jews_1320nonJews_reformatted/imputed935Jews_1320nonJews_L_m-200_chr{chr}.csv'
        if not path.isfile(ibd_file): continue
        df = pd.read_csv(ibd_file, header=0)
        oj_subpop_ibd_df = df.loc[
            (df['id1'].isin(oj_pop_indivs_dict)) | (df['id2'].isin(oj_pop_indivs_dict)) 
        ]
        if chr == 1: 
            for ibd_seg_df in other_oj_pop_ibd_seg_df_dict.values(): 
                ibd_seg_df.columns = list(oj_subpop_ibd_df.columns)
        for other_oj_pop in oj_geo_id_dict.keys():
            #########################################################
            # Filter into the corresponding IBD segments
            if other_oj_pop == oj_pop[1]: continue
            other_oj_subpop_id_dict = {el: '' for el in oj_geo_id_dict[other_oj_pop]}
            filtered_df = oj_subpop_ibd_df.loc[
                (df['id1'].isin(other_oj_subpop_id_dict)) | (df['id2'].isin(other_oj_subpop_id_dict))
            ]
            other_oj_pop_size = len(oj_geo_id_dict[other_oj_pop])
            oj_pop_str = f'{oj_pop[1][0]}-{oj_pop[1][1]}'
            other_oj_pop_str = f'{other_oj_pop[0]}-{other_oj_pop[1]}' if np.nan not in other_oj_pop else 'unknown'
            other_oj_pop_ibd_seg_df_dict[other_oj_pop_str] = \
            pd.concat([other_oj_pop_ibd_seg_df_dict[other_oj_pop_str], filtered_df])
            pop_pair_str = f'{oj_pop_str}OJ_{other_oj_pop_str}OJ'
            ## use the condition "oj_pop[0] > other_oj_pop_str" to avoid replicate pairs, i.e. A & B and B & A
            if oj_pop_size > 0 and other_oj_pop_size > 0 and oj_pop_str > other_oj_pop_str:
                ## Create a string to describe the OJ subgroup-OJ subgroup pair
                if pop_pair_str not in ojsubpop_other_ojsubpop_mean_pws_ibd_len:
                    ## value of the dict is initially a list, [sum_of_ibd_lengths, subpop_size1, subpop_size2]
                    ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str]=[0, oj_pop_size, other_oj_pop_size, {}]
                ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str][0]+=filtered_df['cm_len'].sum()
            ##################################################
            # Output the CSV file storing IBD length between each pair of individuals
            filtered_row_num = filtered_df.shape[0]
            #print(chr, pop_pair_str, filtered_row_num)
            for i in range(filtered_row_num):
                ibd_row = filtered_df.iloc[i]
                # id is in the form of f'{fid}_{iid}', $fid is a digit
                id1=ibd_row['id1']
                id2=ibd_row['id2']
                ibd_len = ibd_row['cm_len']
                ###Test
                #print(ibd_len)
                ###
                if id1 in oj_pop_indivs_dict and id2 in other_oj_subpop_id_dict:
                    ## the prefix should be country name + 'OJ'
                    id2_prefix = f'{oj_id_geo_dict[id2][0]}-{oj_id_geo_dict[id2][1]}OJ'
                    oj_other_oj_df.loc[f'{id2_prefix}__{id2}', id1]+=round(float(ibd_len), 5)
                    ###Test
                    #print(oj_other_oj_df.loc[f'{id2_prefix}__{id2}', id1])
                    ###
                    indiv_pair_key = f'{id1}__{id2}'
                elif id2 in oj_pop_indivs_dict and id1 in other_oj_subpop_id_dict:
                    id1_prefix = f'{oj_id_geo_dict[id1][0]}-{oj_id_geo_dict[id1][1]}OJ'
                    oj_other_oj_df.loc[f'{id1_prefix}__{id1}', id2]+=round(float(ibd_len), 5)
                    indiv_pair_key = f'{id2}__{id1}'
                if pop_pair_str in ojsubpop_other_ojsubpop_mean_pws_ibd_len:
                    #print(pop_pair_str)
                    #print(ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str][3].keys())
                    if not indiv_pair_key in ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str][3]:
                        ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]=0
                    ojsubpop_other_ojsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]+=\
                        round(float(ibd_len), 5)
    ##########################################################################################
    ##########################################################################################
    for other_oj_pop, ibd_seg_df in other_oj_pop_ibd_seg_df_dict.items():
        ibd_seg_df.to_csv(f"IBD_segments_inter_oj_pops/{oj_pop[0]}__{other_oj_pop}OJ_IBD_segs.csv", 
                           encoding='utf-8', index=False)
    # Finally, calculate the mean pairwise IBD length for OJ versus other OJs
    ## after looping through 22 autosomes
    ### This dict should have 14*48 items
    with open("OJ_subpop_inter_mean_pairwise_ibd_length.csv", "a+") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ### This dict should have 14*13/2 items
        for k, v in ojsubpop_other_ojsubpop_mean_pws_ibd_len.items():
            if type(v) is list:
                if len(v)==4:
                    num_non_zero_ibd_indiv_pairs = len(v[3].keys())
                    num_possible_pairs = int(v[1])*int(v[2])
                    num_zeros_to_add = num_possible_pairs - num_non_zero_ibd_indiv_pairs
                    v[3] = list(v[3].values())
                    if num_zeros_to_add > 0:
                        v[3].extend([0]*num_zeros_to_add)
                    pop1=k.split("_")[0]
                    pop2="_".join(tuple(k.split("_")[1:]))
                    mean_v = round(v[0]/num_possible_pairs, 3)
                    # v[3] is a list
                    #print(v[3])
                    ibd_csv_writer.writerow([pop1,pop2,str(mean_v),str(v[1]),str(v[2]),\
                                             round(st.sem(v[3]), 3)])
    #################################################################################################
    ## Additionally, output 14 .csv files storing total IBD length between every pairs of individuals
    ## after looping through 22 autosomes
    ### Reference for conversion from df to csv: 
    ####: https://towardsdatascience.com/how-to-export-pandas-dataframe-to-csv-2038e43d9c03
    #### Write to an individual CSV file
    #### "oj_pop" is a tuple in the form like ("Berlin", "Germany", nan, nan), 
    #### so the second element annotates country
    oj_other_oj_df.round(decimals=4).to_csv(f"oj_pops_inter_IBD/{oj_pop[0]}_IBD.csv", encoding='utf-8')
    

if __name__ == '__main__':
    # warnings.filterwarnings("ignore")
    # load 902 Jews into variables
    df_jews = pd.read_csv('Bray2010_Behar1013_Gladstein2019_Kopelman2020_901Jews.csv', header=0)
    df_otherj = df_jews.loc[lambda x: x['data_source'].str.contains('otherJ$', regex = True)]
    ## add single-origin French Jews to df_otherJ
    df_FrenchJ = df_jews.loc[(df_jews['country1']=="France") & (df_jews['country2'].isnull())]
    ## add single-origin Italian Jews to df_otherJ
    df_ItalianJ = df_jews.loc[df_jews['country1']=="Italy"]
    ## Concat three dfs and drop duplicates
    df_otherj = pd.concat([df_otherj, df_FrenchJ, df_ItalianJ]).drop_duplicates()
    ## Create two dictionaries according to unique geographic origin
    df_oj_geo_groups= df_otherj.groupby(['city1', 'country1']).groups
    ##########################
    # print(df_oj_geo_groups)
    oj_geo_id_dict, oj_id_geo_dict = {}, {}
    
    ## Other-Jews
    for key, val in df_oj_geo_groups.items():
        for ind in val:
            ## avoid mixed individuals
            if not np.nan in key: 
                pass
            fid=df_jews['family_id'].iloc[ind]
            iid=df_jews['individual_id'].iloc[ind]
            if key not in oj_geo_id_dict: oj_geo_id_dict[key]=[]
            cat_id = "_".join((fid, iid))
            oj_geo_id_dict[key].append(cat_id)
            oj_id_geo_dict[cat_id]=key
    
    jewish_pops = []
    ## use this variable to store population names for all Jews
    for key in oj_geo_id_dict:
        if key[0] == "Mosul":
            pop_el = ('KurdistanOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Ar Raqqah':
            pop_el = ('SyrianKurdistanOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Mumbai':
            pop_el = ('IndiaMumbaiOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Thiruvananthapuram':
            pop_el = ('IndiaCochinOJ', key)
            jewish_pops.append(pop_el)
        elif key[0] == 'Tbilisi':
            pop_el = ('GeorgiaOJ', key)
            jewish_pops.append(pop_el)
        else:
            pop_el = (f'{key[1]}OJ', key)
            jewish_pops.append(pop_el)
    
    #print(jewish_pops)
    #print(jews_geo_dict)
    ## Output .csv file for mean pairwise IBD length
    with open("OJ_subpop_inter_mean_pairwise_ibd_length.csv", "w") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ibd_csv_writer.writerow(['pop1', 'pop2', 'mean_pairwise_IBD_length', 'pop1_size', 'pop2_size', \
                                 'IBD_len_SE'])
    
    ## For shorter computing time, 
    ## Run multi-processing using non-mixed subpop OJs
    partial_func = partial(write_oj_oj_pop_pair_ibd, oj_geo_id_dict=oj_geo_id_dict, oj_id_geo_dict=oj_id_geo_dict)
    ## the pool occupies 23 threads(23 sub-populations of OJ in parallel, including the Italian & French "AJ")
    ## There should be 23*(23-1)/2 pairs of OJs to compare IBD
    pool = Pool(23)
    pool.map(partial_func, jewish_pops)
    pool.close()
    pool.join()
```

### In addition, use the following script to detect inter-populatioj individual-wise IBD length between 76 non-Jewish populations


```bash
%%bash
mkdir -p 76_nonJ_pop_inter_IBD -p IBD_segments_inter_76_nonJ_pops
```

```python
from multiprocessing import Pool
import pandas as pd
import csv 
import numpy as np
from os import path
from functools import partial
import scipy.stats as st
import warnings
import re

'''
There're 1320 non-non-Jewish individuals from 76 sub-populations.
This script presents average IBD segment length per pair in Centimorgans.
'''

def write_nj_nj_pop_pair_ibd(nj_pop: str, nj_pop_id_dict: dict, nj_id_pop_dict: dict):
    ## Example of $nj_pop: "French_south" 
    nj_pop_indivs = nj_pop_id_dict[nj_pop]
    nj_pop_indivs_dict = {iid: "" for iid in nj_pop_indivs}
    ## obtain the number of individuals from this non-Jewish population
    nj_pop_size = len(nj_pop_indivs)
    other_nj_row_names = []
    for other_nj_id in nj_id_pop_dict.keys(): 
        if nj_id_pop_dict[other_nj_id] == nj_pop: continue
        other_nj_pop = nj_id_pop_dict[other_nj_id]
        other_nj_row_names.append(f'{other_nj_pop}__{other_nj_id}')
    nj_other_nj_df = pd.DataFrame(columns=nj_pop_indivs, index=other_nj_row_names)
    for nj_id in nj_pop_indivs:
        nj_other_nj_df[nj_id]=0
    # Use this dict the store the following key:value pair
    ## f'{NJ_POP_NAME}__{OTHER_NJ_POP_NAME}'(key): mean pairwise IBD length(float)
    njsubpop_other_njsubpop_mean_pws_ibd_len = {}
    nj_pop_ibd_seg_df_dict = {
        other_nj_pop: pd.DataFrame(columns=[str(i) for i in range(12)]) 
        for other_nj_pop in nj_pop_id_dict.keys() if not nj_pop == other_nj_pop 
    }
    ## Loop through 22 autosomes
    for chr in range(1, 23):
        ### Sum up IBD segment lengths chromosome-wise and obtain the mean value finally after loop
        ibd_file = f'imputed935Jews_1320nonJews_reformatted/imputed935Jews_1320nonJews_L_m-200_chr{chr}.csv'
        if not path.isfile(ibd_file): continue
        df = pd.read_csv(ibd_file, header=0)
        df['id1'] = df['id1'].apply(lambda el: re.sub('^[0-9]+_', '', el))
        df['id2'] = df['id2'].apply(lambda el: re.sub('^[0-9]+_', '', el))
        nj_subpop_ibd_df = df.loc[
            (df['id1'].isin(nj_pop_indivs_dict)) | (df['id2'].isin(nj_pop_indivs_dict)) 
        ]
        if chr == 1: 
            for ibd_seg_df in nj_pop_ibd_seg_df_dict.values(): 
                ibd_seg_df.columns = list(nj_subpop_ibd_df.columns)
        for other_nj_pop in nj_pop_id_dict.keys():
            #########################################################
            # Filter into the corresponding IBD segments
            if other_nj_pop == nj_pop: continue
            other_nj_subpop_id_dict = {el: '' for el in nj_pop_id_dict[other_nj_pop]}
            filtered_df = nj_subpop_ibd_df.loc[
                (df['id1'].isin(other_nj_subpop_id_dict)) | (df['id2'].isin(other_nj_subpop_id_dict))
            ]
            nj_pop_ibd_seg_df_dict[other_nj_pop] = pd.concat([nj_pop_ibd_seg_df_dict[other_nj_pop], filtered_df])
            other_nj_pop_size = len(nj_pop_id_dict[other_nj_pop])
            pop_pair_str = f'{nj_pop}__{other_nj_pop}'
            ## use the condition "nj_pop > other_nj_pop" to avoid replicate pairs, i.e. A & B and B & A
            if nj_pop_size > 0 and other_nj_pop_size > 0 and nj_pop > other_nj_pop:
                ## Create a string to describe the AJ subgroup-OJ subgroup pair
                if pop_pair_str not in njsubpop_other_njsubpop_mean_pws_ibd_len:
                    ## value of the dict is initially a list, [sum_of_ibd_lengths, subpop_size1, subpop_size2]
                    njsubpop_other_njsubpop_mean_pws_ibd_len[pop_pair_str]=[0, nj_pop_size, other_nj_pop_size, {}]
                njsubpop_other_njsubpop_mean_pws_ibd_len[pop_pair_str][0]+=filtered_df['cm_len'].sum()
            ##################################################
            # Output the CSV file storing IBD length between each pair of individuals
            filtered_row_num = filtered_df.shape[0]
            #print(chr, pop_pair_str, filtered_row_num)
            for i in range(filtered_row_num):
                ibd_row = filtered_df.iloc[i]
                # id is in the form of f'{fid}_{iid}', $fid is a digit
                ibd_len = ibd_row['cm_len']
                id1, id2 = ibd_row['id1'], ibd_row['id2']
                if id1 in nj_pop_indivs_dict and id2 in other_nj_subpop_id_dict:
                    ## the prefix should be country name + 'AJ'
                    id2_prefix = nj_id_pop_dict[id2]
                    nj_other_nj_df.loc[f'{id2_prefix}__{id2}', id1]+=round(float(ibd_len), 5)
                    ###Test
                    #print(nj_other_nj_df.loc[f'{id2_prefix}__{id2}', id1])
                    ###
                    indiv_pair_key = f'{id1}__{id2}'
                elif id2 in nj_pop_indivs_dict and id1 in other_nj_subpop_id_dict:
                    id1_prefix = nj_id_pop_dict[id1]
                    nj_other_nj_df.loc[f'{id1_prefix}__{id1}', id2]+=round(float(ibd_len), 5)
                    indiv_pair_key = f'{id2}__{id1}'
                if pop_pair_str in njsubpop_other_njsubpop_mean_pws_ibd_len:
                    #print(pop_pair_str)
                    #print(njsubpop_other_njsubpop_mean_pws_ibd_len[pop_pair_str][3].keys())
                    if not indiv_pair_key in njsubpop_other_njsubpop_mean_pws_ibd_len[pop_pair_str][3]:
                        njsubpop_other_njsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]=0
                    njsubpop_other_njsubpop_mean_pws_ibd_len[pop_pair_str][3][indiv_pair_key]+=\
                        round(float(ibd_len), 5)
    ##########################################################################################
    ##########################################################################################
    for other_nj_pop, ibd_seg_df in nj_pop_ibd_seg_df_dict.items():
        ibd_seg_df.to_csv(f"IBD_segments_inter_76_nonJ_pops/{nj_pop}__{other_nj_pop}_IBD_segs.csv", 
                           encoding='utf-8', index=False)
    # Finally, calculate the mean pairwise IBD length for AJ versus other AJs
    ## after looping through 22 autosomes
    ### This dict should have 14*48 items
    with open("76_nonJ_inter_pop_mean_pairwise_ibd_length.csv", "a+") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ### This dict should have 14*13/2 items
        for k, v in njsubpop_other_njsubpop_mean_pws_ibd_len.items():
            if type(v) is list:
                if len(v)==4:
                    num_non_zero_ibd_indiv_pairs = len(v[3].keys())
                    num_possible_pairs = int(v[1])*int(v[2])
                    num_zeros_to_add = num_possible_pairs - num_non_zero_ibd_indiv_pairs
                    v[3] = list(v[3].values())
                    if num_zeros_to_add > 0:
                        v[3].extend([0]*num_zeros_to_add)
                    pop1=k.split("__")[0]
                    pop2=k.split("__")[1]
                    mean_v = round(v[0]/num_possible_pairs, 3)
                    # v[3] is a list
                    #print(v[3])
                    ibd_csv_writer.writerow([pop1,pop2,str(mean_v),str(v[1]),str(v[2]),\
                                             round(st.sem(v[3]), 3)])
    #################################################################################################
    ## Additionally, output 14 .csv files storing total IBD length between every pairs of individuals
    ## after looping through 22 autosomes
    ### Reference for conversion from df to csv: 
    ####: https://towardsdatascience.com/how-to-export-pandas-dataframe-to-csv-2038e43d9c03
    #### Write to an individual CSV file
    #### "nj_pop" is a string like "French_south"
    nj_other_nj_df.round(decimals=4).to_csv(f"76_nonJ_pop_inter_IBD/{nj_pop}_IBD.csv", encoding='utf-8')
    

if __name__ == '__main__':
    # warnings.filterwarnings("ignore")
    
    ## load 1320 non-Jews from regions of interest to our dict
    nj_pop_id_dict, nj_id_pop_dict = {}, {}
    with open("76pop_non_jews_id_file_list.txt", "r") as f:
        for line in f:
            pop_name=line.strip().split("/")[2].split(".")[0]
            if pop_name not in nj_pop_id_dict: nj_pop_id_dict[pop_name]=[]
            with open(line.strip(), "r") as id_f:
                for line in id_f:
                    iid=line.strip()
                    nj_pop_id_dict[pop_name].append(iid)
                    if iid not in nj_id_pop_dict: nj_id_pop_dict[iid]=pop_name
    
    #print(jewish_pops)
    ## Output .csv file for mean pairwise IBD length
    with open("76_nonJ_inter_pop_mean_pairwise_ibd_length.csv", "w") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ibd_csv_writer.writerow(['pop1', 'pop2', 'mean_pairwise_IBD_length', 'pop1_size', 'pop2_size', \
                                 'IBD_len_SE'])
    
    ## For shorter computing time, 
    partial_func = partial(write_nj_nj_pop_pair_ibd, nj_pop_id_dict=nj_pop_id_dict, nj_id_pop_dict=nj_id_pop_dict)
    ## There should be 76*(76-1)/2 pairs of NJs to compare IBD
    pool = Pool(38)
    pool.map(partial_func, list(nj_pop_id_dict.keys()))
    pool.close()
    pool.join()
```


### In addition, let's investigate intra-population IBD sharing of 76 non-Jewish populations groups using the following python script


```bash
%%bash
mkdir -p 76_nonJ_pop_intra_IBD -p IBD_segments_intra_76_nonJ_pops
```

```python
from multiprocessing import Pool
import pandas as pd
import csv 
import numpy as np
from os import path
from functools import partial
import scipy.stats as st
import warnings
import re

def write_nj_intra_pop_pair_ibd(nj_pop: str, nj_pop_id_dict: dict):
    nj_pop_indivs = nj_pop_id_dict[nj_pop]
    ## obtain the number of individuals from this non-Jewish population
    nj_pop_size = len(nj_pop_indivs)
    ## exclude IBD between two haplotypes of the same individual
    num_indiv_pairs = nj_pop_size*(nj_pop_size-1)/2
    nj_indiv_pair_ibd_dict = {}
    nj_pop_intra_ibd_seg_df = pd.DataFrame(columns=[str(i) for i in range(12)])
    for i in range(nj_pop_size):
        for j in range(i+1, nj_pop_size):
            iid1, iid2 = nj_pop_indivs[i], nj_pop_indivs[j]
            nj_indiv_pair_ibd_dict[f'{iid1}__{iid2}'] = 0
    ## Create a dataframe to store sums of intra-population individual pairwise IBD length
    ### Output the dataframe into a csv file later
    nj_pop_intra_df = pd.DataFrame(columns=['IBD_cM'], index=list(nj_indiv_pair_ibd_dict.keys()))
    nj_pop_intra_df['IBD_cM'] = 0
    #print(nj_pop_intra_df.index)
    ## Loop through 22 autosomes
    for chr in range(1, 23):
        ### Sum up IBD segment lengths chromosome-wise and obtain the mean value finally after loop
        ibd_file = f'imputed935Jews_1320nonJews_reformatted/imputed935Jews_1320nonJews_L_m-200_chr{chr}.csv'
        if not path.isfile(ibd_file): continue
        df = pd.read_csv(ibd_file, header=0)
        df['id1'] = df['id1'].apply(lambda iid: re.sub("[0-9]+_", "", iid))
        df['id2'] = df['id2'].apply(lambda iid: re.sub("[0-9]+_", "", iid))
        pop_intra_df = df.loc[\
                (df['id1'].isin({iid: "" for iid in nj_pop_indivs})) & 
                (df['id2'].isin({iid: "" for iid in nj_pop_indivs})) 
            ]
        if chr == 1: 
            nj_pop_intra_ibd_seg_df.columns = list(pop_intra_df.columns)
        nj_pop_intra_ibd_seg_df = pd.concat([nj_pop_intra_ibd_seg_df, pop_intra_df])
        ibd_row_num = pop_intra_df.shape[0]
        for i in range(ibd_row_num):
            ibd_row = pop_intra_df.iloc[i]
            # id is in the form of f'{fid}_{iid}', $fid is a digit
            iid1=ibd_row['id1']
            iid2=ibd_row['id2']
            ibd_len = ibd_row['cm_len']
            #print(iid1, iid2)
            if f'{iid1}__{iid2}' in nj_indiv_pair_ibd_dict:
                #print("matched")
                nj_pop_intra_df.loc[f'{iid1}__{iid2}', 'IBD_cM']+=round(float(ibd_len), 5)
            elif f'{iid2}__{iid1}' in nj_indiv_pair_ibd_dict:
                #print("matched")
                nj_pop_intra_df.loc[f'{iid2}__{iid1}', 'IBD_cM']+=round(float(ibd_len), 5)
    ## Write the dataframe of selected IBD segments to CSV file
    nj_pop_intra_ibd_seg_df.to_csv(f"IBD_segments_intra_76_nonJ_pops/{nj_pop}_IBD_segs.csv", 
                           encoding='utf-8', index=False)
    ## Write to the general CSV file
    with open("76_nonJ_intra_pop_mean_pairwise_ibd_length.csv", "a") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        mean_v = round(nj_pop_intra_df['IBD_cM'].sum()/num_indiv_pairs, 3)
        ibd_csv_writer.writerow([nj_pop, str(mean_v), nj_pop_size, \
                                 round(st.sem(nj_pop_intra_df['IBD_cM'].tolist()), 3)])
    ## Write to an individual CSV file
    nj_pop_intra_df.round(decimals=4).to_csv(f"76_nonJ_pop_intra_IBD/{nj_pop}_intraIBD.csv", encoding='utf-8')


if __name__ == '__main__':
    #warnings.filterwarnings("ignore")
    ## Non-Jewish ethnic groups in the given 
    nj_pop_id_dict = {}
    with open("76pop_non_jews_id_file_list.txt", "r") as f:
        for line in f:
            pop_name=line.strip().split("/")[2].split(".")[0]
            if pop_name not in nj_pop_id_dict: nj_pop_id_dict[pop_name]=[]
            with open(line.strip(), "r") as id_f:
                for line in id_f:
                    iid=line.strip()
                    nj_pop_id_dict[pop_name].append(iid)
    
    ## Output .csv file for mean pairwise IBD length
    with open("76_nonJ_intra_pop_mean_pairwise_ibd_length.csv", "w") as mean_pairwise_ibd_f:
        ibd_csv_writer = csv.writer(mean_pairwise_ibd_f)
        ibd_csv_writer.writerow(['pop', 'mean_pairwise_IBD_length', 'pop_size', 'IBD_len_SE'])
    
    ## For shorter computing time, 
    ### run multi-processing for 76 populations occuyping 36 threads
    partial_func = partial(write_nj_intra_pop_pair_ibd, nj_pop_id_dict=nj_pop_id_dict)
    pool = Pool(38)
    pool.map(partial_func, list(nj_pop_id_dict.keys()))
    pool.close()
    pool.join()
```
