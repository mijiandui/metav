
metav is a small pipeline to perfom virus identification from raw metagenomics NGS data. 

1. Quality control raw reads with fastp & remove host contamination with BWA.  
2. Assembly use either megahit or MetaSpades.
3. Viral identification with virsorter2&checkv SOP and DeepVirFinder.

# Dependency

* install prerequisites (ubuntu20.04)
```
apt update 
apt install -yyq make cmake gcc g++ libisal-dev libdeflate-dev zlib1g-dev wget git libncurses5-dev libbz2-dev liblzma-dev libcurl4-openssl-dev parallel 
```
* install miniconda3 (Ignore this if you already have `conda` installed)
```
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh
sh ./Miniconda3-py38_4.12.0-Linux-x86_64.sh
conda config --add channels https://mirrors.aliyun.com/anaconda/pkgs/main/
conda config --add channels https://mirrors.aliyun.com/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.aliyun.com/anaconda/cloud/bioconda/
```

# Installation

```
tar -xvzf metav.tgz 
cd metav
make clean && make -j
```

# Setup virsorter2 and checv databases

```
conda activate vs2

mkdir database
cd database
virsorter setup -d db-vs2 -j 4
checkv download_database .
cd -
conda deactivate


# this is optional, you can specify the db path from command line
export CHECKV_DB=$PWD/database/checkv-db-v1.4 
```

# Change default setting of virsorter2

```
conda activate vs2

# please replace /ssd-cache with the mounting point of ssd disk in your system.
# this will modify file: $HOME/miniconda3/envs/vs2/lib/python3.10/site-packages/virsorter/template-config.yaml
virsorter config LOCAL_SCRATCH=/ssd-cache
virsorter config --set HMMSEARCH_THREADS=4
virsorter config --set CONTIG_BP_PER_SPLIT=4000000
virsorter config --set FAA_BP_PER_SPLIT=1000000
virsorter config --set GFF_SEQNUM_PER_SPLIT=50000
conda deactivate
```

# Quick start
```
export PATH=$PWD/src:$PATH
metav --help
```


## Build BWA index for host genome
```
./third_party_tools/bwa/bwa index host_genome.fasta
```
This may take a while (~1h for human genome). You can add any sequences you want to remove to the host_genome FASTA file.


## Run `metav` using assembler ``megahit`` (default)
```
metav -1 ../data/QH616240FF5032_NDSW51265_1.fq.gz \
      -2 ../data/QH616240FF5032_NDSW51265_2.fq.gz \
      --host ../database/GCF_000003025.6_Sscrofa11.1_genomic.fna \
      --min-contig-len 200 \            
      -d ../database/checkv-db-v1.4 \    # path of checkv database 
      -t 72 \                            # number of cpus 
      -o QH616240FF5032_NDSW51265_k1 \   # output directory
      --k-list 59 \  # specify kmer list eg: --k-list 21,39,59,79
      -p QH616240FF5032_NDSW51265   # outprefix
```

## Run `metav` using assembler ``MetaSpades``
```
metav -1 ../data/QH616240FF5032_NDSW51265_1.fq.gz \
      -2 ../data/QH616240FF5032_NDSW51265_2.fq.gz \
      --host ../database/GCF_000003025.6_Sscrofa11.1_genomic.fna \
      -a metaspades \ # specify assembler
      -d ../database/checkv-db-v1.4 \
      -t 72 \
      -o QH616240FF5032_NDSW51265_spades \
      -p QH616240FF5032_NDSW51265
``` 
* If you want to run MetaSpades without read correction, please specify option ``--only-assembly ``; 
* You can change default kmer used by MetaSpades via `--k-list`.
  

## Continue mode
Resume an interrupted job from out:
```
metav -o out --continue -f
```

# Usage and Options 
```
Usage:
  metav [options] {{-1 <pe1> -2 <pe2> }} {{--host <host reference genome>}} {{-d <checkv_db>}} [-o <out_dir>] [-p <output_prefix>] [Options]

Input options:
    -1                       <pe1>          fasta/q paired-end #1 file, paired with file in <pe2> (supporting plain text and gz extensions)
    -2                       <pe2>          fasta/q paired-end #2 file, paired with file in <pe1> (supporting plain text and gz extensions)
    --host                                  host genome in FASTA format(with BWA index available)
    -d/--db-dir              <dir>          path of checkv database(can also be set using export CHECKV_DB=<db_path>)
   

Optional Arguments:

  Basic quality control options: 
    -k/--min-seed-len        <int>          minimun seed length for a bwa hit [31]                     

  Basic assembly options:
    -a/--assembler           <str>          assembler used, selected from <megahit,metaspades>, [megahit]
    
    --k-list                 <int,int,..>   comma-separated list of kmer size used by megahit
                                            all must be odd, in the range 15-{0}, increment <= 28)
                                            [21,39,59,79,99]
  Megahit options:
    --min-count              <int>          minimum multiplicity for filtering (k_min+1)-mers used by megahit [2]
    --min-contig-len         <int>          minimum length of contigs to output [500]
    --no-mercy                              do not add mercy kmers in megahit assembling
    
  MetaSpades options:
    --only-assembly                         run only assembling (without read error correction) with metaspades
  
  Basic viral identification options:
    -l/--min-length           <int>         minimal contig length to viral identification, [1500]
    --include-groups          <str>         classifier of viral groups to include(comma seperated and no space in between);
                                            available options are dsDNAphage,NCLDV,RNA,ssDNA,ssDNA, [dsDNAphage,ssDNA]
    --min-score               <float>       minimal score to be identified as viral, [0.5]
    --dvf-score               <float>       DeepVirFinder score to be identified as viral [0.95]

  Hardware options:
    -t/--num-cpu-threads     <int>          number of CPU threads [# of logical processors]


  Output options:
    -o/--out-dir             <string>       output directory [./megahit_out]
    -p/--out-prefix          <string>       output prefix, eg. sample name
    -f/--force                              overwritten exist directory
    --keep-tmp-files                        keep all temporary files
    --tmp-dir                <string>       set temp directory

Other Arguments:
    --continue                              continue a METAV run from its last available check point.
                                            please set the output directory correctly when using this option.
    -h/--help                               print the usage message
    -v/--version                            print version

```
