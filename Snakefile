# Snakefile for generating the input datasets.
# It downloads reference FASTAs and simulates reads from them.
#
# This creates the downloads/, datasets/ and genomes/ directories.
# (downloads/ can be deleted.)

# Additional 'genomes' that can be added to this list:
# - ecoli50 (fifty E. coli genomes)
# - chrY (chromosome Y of CHM13)

GENOMES = ("fruitfly", "maize", "CHM13", "rye", "chrY")
N_READS = {
    50: 1_000_000,
    75: 1_000_000,
    100: 1_000_000,
    150: 1_000_000,
    200: 1_000_000,
    300: 1_000_000,
    500: 1_000_000,
    1000: 500_000,
    5000: 100_000,
    10000: 50_000,
}

GENOMES = ("ecoli")
N_READS = {rl: n // 100 for rl, n in N_READS.items()}

LONG_READ_LENGTHS = tuple(n for n in N_READS if n >= 1000)  # single-end only
READ_LENGTHS = tuple(n for n in N_READS if n < 1000)
MODELS = {"clr": "data/pbsim3/QSHMM-RSII.model", "ont": "data/pbsim3/QSHMM-ONT-HQ.model", "hifi": "data/pbsim3/QSHMM-RSII.model"}
DATASETS = expand("{genome}-{read_length}", genome=GENOMES, read_length=READ_LENGTHS)
LONG_DATASETS = expand("{genome}-{read_length}", genome=GENOMES, read_length=LONG_READ_LENGTHS)
ENDS = ("pe", "se")

VARIATION_SETTINGS = {
    "sim1": "",
    "sim3": "--snp-rate 0.001 --small-indel-rate 0.0001 --max-small-indel-size 50",
    "sim4": "--snp-rate 0.005 --small-indel-rate 0.0005 --max-small-indel-size 50",
    "sim5": "--snp-rate 0.005 --small-indel-rate 0.001 --max-small-indel-size 100",
    "sim6": "--snp-rate 0.05 --small-indel-rate 0.002 --max-small-indel-size 100",
}
SIM = ["sim0", "sim0p1"] + list(VARIATION_SETTINGS)
SIM = ["sim1", "sim3"]
LONG_SIM = ["ont", "hifi", "clr"]


wildcard_constraints:
    read_length=r"\d{2,3}",
    long_read_length=r"\d{4,5}",
    sim01=r"sim(0|0p1)",
    read_type=r"illumina|ont|clr|hifi",
    long_read_type=r"ont|clr|hifi"


localrules:
    download_fruitfly, download_maize, download_chm13, download_rye, download_ecoli50, filter_ecoli50, filter_fruitfly, clone_seqan, samtools_faidx


rule:
    input:
        expand("datasets/{sim}-illumina/{ds}/{r}.fastq.gz", sim=SIM, ds=DATASETS, r=(1, 2)),
        expand("datasets/{sim}-illumina/{ds}/truth.bam", sim=SIM, ds=DATASETS + LONG_DATASETS),
        expand("datasets/{sim}-{read_type}/{ds}/1.fastq.gz", sim=SIM, read_type=["ont", "clr"], ds=LONG_DATASETS),
        # expand("datasets/{sim}-{read_type}/{ds}/truth.maf.gz", sim=list(VARIATION_SETTINGS), read_type=["ont", "clr"], ds=LONG_DATASETS)
        expand("datasets/{sim}-{read_type}/{ds}/truth.bed", sim=SIM, read_type=["ont", "clr", "hifi"], ds=LONG_DATASETS)

# Download genomes

rule download_fruitfly:
    output: "downloads/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz"
    shell:
        "curl ftp://ftp.ensembl.org/pub/release-97/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz > {output}"

rule download_maize:
    output: "downloads/Zm-B73-REFERENCE-NAM-5.0.fa.gz"
    shell:
        "curl https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz > {output}"

rule download_chm13:
    output: "downloads/chm13v2.0.fa.gz"
    shell:
        "curl https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz > {output}"

rule download_rye:
    output: "downloads/GCA_016097815.1_HAU_Weining_v1.0_genomic.fna.gz"
    shell:
        "curl https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Secale_cereale/latest_assembly_versions/GCA_016097815.1_HAU_Weining_v1.0/GCA_016097815.1_HAU_Weining_v1.0_genomic.fna.gz > {output}"

rule download_ecoli50:
    output:
        "downloads/ecoli50/done.txt"
    input: "ecoli-accessions.txt"
    shell:
        "mkdir -p downloads/ecoli50; "
        "head -n 50 ecoli-accessions.txt > downloads/ecoli50/accessions.txt; "
        "ncbi-genome-download -o downloads/ecoli50 -A downloads/ecoli50/accessions.txt --formats fasta --flat-output bacteria; "
        "touch downloads/ecoli50/done.txt"

rule filter_ecoli50:
    output: "genomes/ecoli50.fa"
    input: "downloads/ecoli50/done.txt"
    shell:
        "python noplasmids.py downloads/ecoli50/*.fna.gz > {output}"


# Uncompress and filter downloaded genomes

rule filter_fruitfly:
    output: "genomes/fruitfly.fa"
    input: rules.download_fruitfly.output
    shell:
        """
        zcat {input} > {output}.tmp.fa
        samtools faidx {output}.tmp.fa
        # Discard contigs shorter than 16 kbp
        awk '$2>=16000 {{print $1}}' {output}.tmp.fa.fai > {output}.tmp.regions.txt
        samtools faidx -r {output}.tmp.regions.txt {output}.tmp.fa > {output}.tmp2.fa
        mv {output}.tmp2.fa {output}

        rm {output}.tmp.fa {output}.tmp.fa.fai {output}.tmp.regions.txt
        """

rule uncompress_maize:
    output: "genomes/maize.fa"
    input: rules.download_maize.output
    shell:
        "zcat {input} > {output}"

rule uncompress_chm13:
    output: "genomes/CHM13.fa"
    input: rules.download_chm13.output
    shell:
        "zcat {input} > {output}"

rule filter_rye:
    output: "genomes/rye_raw.fa"
    input: rules.download_rye.output
    shell:
        """
        zcat {input} > {output}.tmp.fa
        samtools faidx {output}.tmp.fa
        # Discard contigs shorter than 50000
        awk '$2>=50000 {{print $1}}' {output}.tmp.fa.fai > {output}.tmp.regions.txt
        samtools faidx -r {output}.tmp.regions.txt {output}.tmp.fa > {output}.tmp2.fa
        mv {output}.tmp2.fa {output}

        rm {output}.tmp.fa {output}.tmp.fa.fai {output}.tmp.regions.txt
        """

rule chunk_rye:
    input:
        ref="genomes/rye_raw.fa"
    output:
        chunked="genomes/rye.fa"
    params:
        max_chunk_size = 1_000_000_000
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        with open(output.chunked, 'w') as out_handle:
            for record in SeqIO.parse(input.ref, "fasta"):
                seq_len = len(record.seq)

                if seq_len <= params.max_chunk_size:
                    SeqIO.write(record, out_handle, "fasta")
                else:
                    num_chunks = (seq_len + params.max_chunk_size - 1) // params.max_chunk_size

                    for i in range(num_chunks):
                        start = i * params.max_chunk_size
                        end = min((i + 1) * params.max_chunk_size, seq_len)

                        chunk_seq = record.seq[start:end]
                        chunk_id = f"{record.id}_chunk{i+1}of{num_chunks}"
                        chunk_desc = f"{record.description} | chunk {i+1}/{num_chunks} | pos {start+1}-{end}"

                        chunk_record = SeqRecord(
                            chunk_seq,
                            id=chunk_id,
                            description=chunk_desc
                        )

                        SeqIO.write(chunk_record, out_handle, "fasta")


rule extract_chry:
    output: "genomes/chrY.fa"
    input: fasta="genomes/CHM13.fa", fai="genomes/CHM13.fa.fai"
    shell:
        "samtools faidx {input.fasta} chrY > {output}"


# Generate simulated reads and BAM files with expected alignments (truth)

rule mason_variator:
    output:
        vcf="variants/{sim,sim[1-9]}-{genome}.vcf",
        fasta="variants/{sim,sim[1-9]}-{genome}.tmp.fa"
    input:
        fasta="genomes/{genome}.fa",
        fai="genomes/{genome}.fa.fai",
        mason_variator="bin/mason_variator",
        mason_materializer="bin/mason_materializer"
    params:
        variation_settings=lambda wildcards: VARIATION_SETTINGS[wildcards.sim]
    shell:
        """
        if [ "{wildcards.sim}" = "sim1" ]; then
            touch {output.vcf}
            ln -sfr {input.fasta} {output.fasta}
        else
            {input.mason_variator} -ir {input.fasta} {params.variation_settings} -ov {output.vcf}.tmp.vcf
            mv -v {output.vcf}.tmp.vcf {output.vcf}
            {input.mason_materializer} -ir {input.fasta} -iv {output.vcf} -o {output.fasta}
        fi
        """

rule allele_removal:
    input:
        fasta="variants/{sim,sim[2-9]}-{genome}.tmp.fa"
    output:
        fasta="variants/{sim,sim[2-9]}-{genome}.fa"
    run:
        from Bio import SeqIO

        with open(input.fasta) as fa_in, open(output.fasta, "w") as fa_out:
            for record in SeqIO.FastaIO.FastaIterator(fa_in):
                if not record.id.endswith("/1"):
                    raise ValueError(f"{record.id} does not end with /1, check if {input.fasta} was produced by mason_materializer")
                fa_out.write(f">{record.id[:-2]}\n{record.seq}\n")


rule sim1_allele_removal:
    input:
        fasta="genomes/{genome}.fa"
    output:
        fasta="variants/{sim,sim1}-{genome}.fa"
    shell:
        "ln -sfr {input.fasta} {output.fasta}"


def mason_simulator_parameters(wildcards):
    read_length = int(wildcards.read_length)
    result = f"--illumina-read-length {read_length}"
    if read_length >= 250:
        result += " --fragment-mean-size 700"
    return result


rule mason_simulator:
    output:
        r1_fastq="datasets/{sim,sim[1-9]}-illumina/{genome}-{read_length}/1.fastq.gz",
        r2_fastq="datasets/{sim,sim[1-9]}-illumina/{genome}-{read_length}/2.fastq.gz",
        bam="datasets/{sim,sim[1-9]}-illumina/{genome}-{read_length}/truth.bam"
    input:
        fasta="genomes/{genome}.fa",
        vcf="variants/{sim}-{genome}.vcf",
        mason_simulator="bin/mason_simulator"
    params:
        extra=mason_simulator_parameters,
        n_reads=lambda wildcards: N_READS[int(wildcards.read_length)],
        vcf_arg=lambda wildcards: "" if wildcards.sim == "sim1" else f"-iv variants/{wildcards.sim}-{wildcards.genome}.vcf"
    log: "logs/mason_simulator/{sim}-{genome}-{read_length}.log"
    shell:
        "ulimit -n 16384"  # Avoid "Uncaught exception of type MasonIOException: Could not open right/single-end output file."
        "\n{input.mason_simulator}"
        " --num-threads 1"  # Output depends on number of threads, leave at 1 for reproducibility
        " -ir {input.fasta}"
        " -n {params.n_reads}"
        " {params.vcf_arg}"
        " {params.extra}"
        " -o {output.r1_fastq}.tmp.fastq.gz"
        " -or {output.r2_fastq}.tmp.fastq.gz"
        " -oa {output.bam}.tmp.bam"
        " 2>&1 | tee {log}"
        "\nmv -v {output.r1_fastq}.tmp.fastq.gz {output.r1_fastq}"
        "\nmv -v {output.r2_fastq}.tmp.fastq.gz {output.r2_fastq}"
        "\nmv -v {output.bam}.tmp.bam {output.bam}"


rule mason_simulator_long:
    output:
        fastq="datasets/{sim,sim[1-9]}-illumina/{genome}-{long_read_length}/1.fastq.gz",
        bam="datasets/{sim,sim[1-9]}-illumina/{genome}-{long_read_length}/truth.bam"
    input:
        fasta="genomes/{genome}.fa",
        vcf="variants/{sim}-{genome}.vcf",
        mason_simulator="bin/mason_simulator"
    params:
        n_reads=lambda wildcards: N_READS[int(wildcards.long_read_length)],
        fragment_length=lambda wildcards: int(int(wildcards.long_read_length) * 1.5),
        vcf_arg=lambda wildcards: "" if wildcards.sim == "sim1" else f"-iv variants/{wildcards.sim}-{wildcards.genome}.vcf"
    log: "logs/mason_simulator/{sim}-{genome}-{long_read_length}.log"
    shell:
        "ulimit -n 16384"  # Avoid "Uncaught exception of type MasonIOException: Could not open right/single-end output file."
        "\n{input.mason_simulator}"
        " --num-threads 1"  # Output depends on number of threads, leave at 1 for reproducibility
        " --illumina-read-length {wildcards.long_read_length}"
        " --fragment-mean-size {params.fragment_length}"
        " -ir {input.fasta}"
        " -n {params.n_reads}"
        " {params.vcf_arg}"
        " -o {output.fastq}.tmp.fastq.gz"
        " -oa {output.bam}.tmp.bam"
        " 2>&1 | tee {log}"
        "\nmv -v {output.fastq}.tmp.fastq.gz {output.fastq}"
        "\nmv -v {output.bam}.tmp.bam {output.bam}"


def readsimulator_parameters(wildcards):
    read_length = int(wildcards.read_length)
    if read_length >= 250:
        return " --mean-insert-size 700"
    return ""


rule sim01:
    output:
        r1_fastq="datasets/{sim01}-illumina/{genome}-{read_length}/1.fastq.gz",
        r2_fastq="datasets/{sim01}-illumina/{genome}-{read_length}/2.fastq.gz",
        bam="datasets/{sim01}-illumina/{genome}-{read_length}/truth.bam"
    input:
        fasta="genomes/{genome}.fa",
    params:
        extra=readsimulator_parameters,
        n_reads=lambda wildcards: N_READS[int(wildcards.read_length)],
        error_rate=lambda wildcards: {"sim0": 0.0, "sim0p1": 0.1}[wildcards.sim01]
    shell:
        "python readsimulator.py{params.extra} -e {params.error_rate} -n {params.n_reads} --read-length {wildcards.read_length} {input.fasta} | samtools view -o {output.bam}.tmp.bam"
        "\nsamtools fastq -N -1 {output.r1_fastq} -2 {output.r2_fastq} {output.bam}.tmp.bam"
        "\nmv {output.bam}.tmp.bam {output.bam}"


rule sim01_long:
    output:
        fastq="datasets/{sim01}-illumina/{genome}-{long_read_length}/1.fastq.gz",
        bam="datasets/{sim01}-illumina/{genome}-{long_read_length}/truth.bam"
    input:
        fasta="genomes/{genome}.fa",
    params:
        n_reads=lambda wildcards: N_READS[int(wildcards.long_read_length)],
        error_rate=lambda wildcards: {"sim0": 0.0, "sim0p1": 0.1}[wildcards.sim01]
    shell:
        "python readsimulator.py --se -e {params.error_rate} -n {params.n_reads} --read-length {wildcards.long_read_length} {input.fasta} | samtools view -o {output.bam}.tmp.bam"
        "\nsamtools fastq -N -0 {output.fastq} {output.bam}.tmp.bam"
        "\nmv {output.bam}.tmp.bam {output.bam}"


def pbsim_parameters(wildcards):
    mean_read_length = int(wildcards.long_read_length)
    result = "--length-mean {}".format(mean_read_length)

    reference_path = "genomes/" + wildcards.genome + ".fa"
    ref_len = 0
    with open(reference_path, "r") as ref:
        for line in ref:
            if not line.startswith('>'):
                ref_len += len(line.strip())
    num_reads = N_READS[mean_read_length]
    depth = float(num_reads * mean_read_length) / float(ref_len)
    result += " --depth {}".format(depth)

    if wildcards.read_type == "hifi":
        result += " --pass-num 10"
    return result


def pbsim_outprefix(wildcards):
    return f"datasets/{wildcards.sim}-{wildcards.read_type}/{wildcards.genome}-{wildcards.long_read_length}/tmp"

def first_bam_name(wildcards):
    if wildcards.read_type == "hifi":
        return pbsim_outprefix(wildcards) + "_0001.bam"
    return []


rule pbsim:
    output:
        maf="datasets/{sim,sim[1-9]}-{read_type,clr|ont}/{genome}-{long_read_length}/truth.maf.gz",
        fastq="datasets/{sim,sim[1-9]}-{read_type,clr|ont}/{genome}-{long_read_length}/1.fastq.gz"
    input:
        fasta="variants/{sim,sim[1-9]}-{genome}.fa",
        model=lambda wildcards: MODELS[wildcards.read_type]
    params:
        extra=pbsim_parameters,
        outprefix=pbsim_outprefix,
        outid="S"
    log: "logs/pbsim3/{sim}-{read_type}-{genome}-{long_read_length}.log"
    shell:
        "pbsim"
        " --strategy wgs"
        " --genome {input.fasta}"
        " --method qshmm"
        " --qshmm {input.model}"
        " --prefix {params.outprefix}"
        " --id-prefix {params.outid}"
        " {params.extra}"
        " --length-sd 0"
        "\ncat {params.outprefix}_*.fq.gz > {output.fastq}"
        "\ncat {params.outprefix}_*.maf.gz > {output.maf}"
        "\nrm {params.outprefix}_*.ref"
        "\nrm {params.outprefix}_*.maf.gz"
        "\nrm {params.outprefix}_*.fq.gz"


rule pbsim_hifi:
    output:
        maf="datasets/{sim,sim[1-9]}-{read_type,hifi}/{genome}-{long_read_length}/truth.maf.gz",
        bam=temp("datasets/{sim,sim[1-9]}-{read_type,hifi}/{genome}-{long_read_length}/1.bam")
    input:
        fasta="variants/{sim,sim[1-9]}-{genome}.fa",
        model=lambda wildcards: MODELS[wildcards.read_type]
    params:
        extra=pbsim_parameters,
        outprefix=pbsim_outprefix,
        outid="S"
    log: "logs/pbsim3/{sim}-{read_type}-{genome}-{long_read_length}.log"
    shell:
        "pbsim"
        " --strategy wgs"
        " --genome {input.fasta}"
        " --method qshmm"
        " --qshmm {input.model}"
        " --prefix {params.outprefix}"
        " --id-prefix {params.outid}"
        " {params.extra}"
        " --length-sd 0"
        "\ncat {params.outprefix}_*.maf.gz > {output.maf}"
        "\nrm {params.outprefix}_*.ref"
        "\nrm {params.outprefix}_*.maf.gz"
        "\nsamtools merge -o {output.bam} {params.outprefix}_*.bam"
        "\nrm {params.outprefix}_*.bam"


rule ccs:
    output:
        fastq="datasets/{sim,sim[1-9]}-hifi/{genome}-{long_read_length}/1.fastq.gz"
    input:
        bam="datasets/{sim}-hifi/{genome}-{long_read_length}/1.bam"
    log:
        "datasets/{sim}-hifi/{genome}-{long_read_length}/ccs.log"
    threads:
        32
    shell:
        """
        ccs --log-file {log} -j {threads} {input.bam} {output.fastq}
        """


# Add allele fields to vcf for g2gtools
rule convert_vcf:
    input: 
        vcf="variants/{sim,sim[1-9]}-{genome}.vcf"
    output: 
        vcf="variants/{sim,sim[1-9]}-{genome}-g2g.vcf"
    run:
        from pysam import VariantFile

        vcf_in = VariantFile(input.vcf, 'r')
        vcf_in.header.formats.add("GT", "1", "String", "Genotype")
        vcf_in.header.formats.add("GQ", "1", "Integer", "Genotype Quality")
        vcf_in.header.formats.add("DP", "1", "Integer", "Read Depth")

        vcf_out = VariantFile(output.vcf, 'w', header=vcf_in.header)

        for record in vcf_in:
            for sample in record.samples:
                record.samples[sample]['GT'] = (1, 1)  # 1/1 genotype
                record.samples[sample]['GQ'] = 99
                record.samples[sample]['DP'] = 10

            vcf_out.write(record)


# Make inverse vcf from the variated reference to the original for the ground truth liftover
# Taken from https://github.com/samtools/bcftools/issues/2096
rule invert_vcf:
    input:
        vcf="variants/{sim,sim[1-9]}-{genome}-g2g.vcf"
    output:
        vcf="variants/{sim,sim[1-9]}-{genome}-g2g-inverse.vcf"
    run:
        from collections import defaultdict

        with open(output.vcf, "w") as f_out, open(input.vcf, "r") as f_in:
            offset_counter = defaultdict(int)
            for line in f_in:
                if line.startswith("#"):
                    print(line.strip(), file=f_out)
                    continue

                fields = line.strip().split("\t")
                ref = fields[3]
                alt = fields[4]
                pos = int(fields[1])
                chrom = fields[0]
                offset = offset_counter[chrom]
                new_alt = ref
                new_ref = alt
                new_pos = pos + offset
                record_offset = len(alt) - len(ref)
                offset_counter[chrom] += record_offset
                fields[1] = str(new_pos)
                fields[3] = new_ref
                fields[4] = new_alt
                print("\t".join(fields), file=f_out)


rule sort_vcf:
    input:
        vcf="variants/{sim,sim[1-9]}-{genome}-g2g-inverse.vcf"
    output:
        vcf="variants/{sim,sim[1-9]}-{genome}-g2g-inverse.vcf.gz",
        tbi="variants/{sim,sim[1-9]}-{genome}-g2g-inverse.vcf.gz.tbi"
    shell:
        """
        bcftools sort -m 2G {input.vcf} | bgzip --stdout > {output.vcf}
        tabix {output.vcf}
        """


def parse_maf_alignment(alignment, ref_index_to_name: dict, consensus_names=None):
    if len(alignment) != 2:
        raise ValueError(f"{maf_path} should contain 2 alignments per entry, please make sure that it was produced by pbsim3")

    ref_seq = None
    read_seq = None

    for seq_record in alignment:
        if seq_record.id == "ref":
            ref_seq = seq_record
        elif seq_record.id.startswith("S"):
            read_seq = seq_record

    if ref_seq is None or read_seq is None:
        raise ValueError(f"""{alignment} entry in {maf_path} should contain entries for \"ref\" 
                                and \"S<>_<>\", please make sure that {maf_path} was produced by pbsim3""")

    read_name = read_seq.id
    ref_index = None
    read_index = None
    if not consensus_names:
        match = re.match(r'S(\d+)_(\d+)', read_name)
        if not match:
            raise ValueError(f"Read name {read_name} does not match expected format S<reference_index>_<read id>")
        ref_index = int(match.group(1))
    else:
        match = re.match(r'S(\d+)/(\d+)/(\d+)', read_name)
        if not match:
            raise ValueError(f"Read name {read_name} does not match expected format S<reference id>/<read id>/<pass id>")
        ref_index = int(match.group(1))
        read_index = int(match.group(2))
    
    if ref_index - 1 not in ref_index_to_name:
        raise ValueError(f"Reference index {ref_index - 1} not found in reference index file")
    ref_name = ref_index_to_name[ref_index - 1]

    ref_start = ref_seq.annotations["start"]
    ref_end = ref_start + ref_seq.annotations["size"]

    if not consensus_names:
        query_name = read_name
    else:
        query_name = f"S{ref_index}/{read_index}/ccs"
        if query_name not in consensus_names:
            return None
    # if not query_name.endswith("/1"):
    #     query_name += "/1"
    # print(read_seq.__dict__)
    query_or = "+" if read_seq.annotations["strand"] == 1 else "-"

    return query_name, query_or, ref_name, ref_start, ref_end


rule long_read_truth:
    input:
        fai="variants/{sim,sim[1-9]}-{genome}.fa.fai",
        maf="datasets/{sim,sim[1-9]}-{read_type,clr|ont|hifi}/{genome}-{long_read_length}/truth.maf.gz",
        fastq="datasets/{sim,sim[1-9]}-{read_type,clr|ont|hifi}/{genome}-{long_read_length}/1.fastq.gz"
    output:
        bed="datasets/{sim,sim[1-9]}-{read_type,clr|ont|hifi}/{genome}-{long_read_length}/truth-unlifted.bed"
    params:
        consensus_names=lambda wildcards: 0 if wildcards.read_type != "hifi" else 1
    run:
        from xopen import xopen
        from Bio import AlignIO
        from pysam import FastxFile

        read_positions = {}

        ref_index_to_name = {}
        with open(input.fai) as fai_file:
            for line_num, line in enumerate(fai_file):
                ref_name = line.strip().split('\t')[0]
                ref_index_to_name[line_num] = ref_name
        print(ref_index_to_name)

        consensus_names = None
        if params.consensus_names == 1:
            with FastxFile(input.fastq) as fastq_file:
                consensus_names = {read.name for read in fastq_file}

        with xopen(input.maf) as maf_file, open(output.bed, "w") as bed_file:
            alignments = AlignIO.parse(maf_file, "maf")
            for alignment in alignments:
                parsed_maf_alignment = parse_maf_alignment(alignment, ref_index_to_name, consensus_names)
                if parsed_maf_alignment:
                    query_name, query_or, ref_name, ref_start, ref_end = parsed_maf_alignment
                    bed_file.write(f"{ref_name}\t{ref_start}\t{ref_end}\t{query_name}\t0\t{query_or}\n")


rule sim1_gt:
    input:
        bed="datasets/{sim,sim1}-{long_read_type}/{genome}-{long_read_length}/truth-unlifted.bed"
    output:
        bed="datasets/{sim,sim1}-{long_read_type}/{genome}-{long_read_length}/truth.bed"
    shell:
        """
        cp {input.bed} {output.bed}
        """


rule ref_transform:
    input:
        fasta="variants/{sim,sim[2-9]}-{genome}.fa",
        vcf="variants/{sim,sim[2-9]}-{genome}-g2g-inverse.vcf.gz",
        tbi="variants/{sim,sim[2-9]}-{genome}-g2g-inverse.vcf.gz.tbi"
    output:
        # TODO check if the double-inverted fasta is the same as the original
        fasta="variants/{sim,sim[2-9]}-{genome}-inverted.fa",
        vci="variants/{sim,sim[2-9]}-{genome}-inverse.vci.gz"
    params:
        outprefix="variants/{sim,sim[2-9]}-{genome}-inverted",
        outvci="variants/{sim,sim[2-9]}-{genome}-inverse.vci"
    singularity:
        "docker://churchilllab/g2gtools:3.0.0"
    shell:
        """
        g2gtools vcf2vci --vci {params.outvci} --vcf {input.vcf} --fasta {input.fasta} --strain simulated -v
        g2gtools patch --fasta {input.fasta} --vci {output.vci} --out {params.outprefix}.patched.fa
        g2gtools transform --fasta {params.outprefix}.patched.fa --vci {output.vci} --out {output.fasta}
        """


rule gt_liftover:
    input:
        gt="datasets/{sim,sim[2-9]}-{long_read_type}/{genome}-{long_read_length}/truth-unlifted.bed",
        vci="variants/{sim,sim[2-9]}-{genome}-inverse.vci.gz"
        # fasta="variants/{sim,sim[1-9]}-{genome}.fa"
    output:
        gt="datasets/{sim,sim[2-9]}-{long_read_type}/{genome}-{long_read_length}/truth.bed"
        # fasta="datasets/{dataset_id}/ref.fa"
    # params:
        # unmapped="datasets/{dataset_id}/ref_ground_truth.bed.unmapped"
    singularity:
        "docker://churchilllab/g2gtools:3.0.0"
    shell:
        """
        g2gtools convert --in {input.gt} --vci {input.vci} -f BED --out {output.gt}
        """

# Misc

rule samtools_faidx:
    output: "{genome}.fa.fai"
    input: "{genome}.fa"
    shell: "samtools faidx {input}"


# Build our own Mason binaries because the one from Conda crashes
rule clone_seqan:
    output: "seqan/cloned"
    shell:
        "git clone https://github.com/seqan/seqan.git"
        "; ( cd seqan && git checkout seqan-v2.5.0rc2 )"
        "; touch seqan/cloned"


rule build_mason:
    output: "bin/mason_variator", "bin/mason_simulator", "bin/mason_materializer"
    input: "seqan/cloned"
    threads: 99
    shell:
        "cmake -DSEQAN_BUILD_SYSTEM=APP:mason2 -DSEQAN_ARCH_SSE4=1 -B build-seqan seqan; "
        "cmake --build build-seqan -j {threads}; "
        "mv build-seqan/bin/mason_simulator build-seqan/bin/mason_variator build-seqan/bin/mason_materializer bin/"
        #"; rm -r seqan"


rule download_pbsim_models:
    output: "data/pbsim3/QSHMM-ONT-HQ.model", "data/pbsim3/QSHMM-RSII.model"
    params:
        outdir = "data/pbsim3"
    threads: 99
    shell:
        "git clone https://github.com/yukiteruono/pbsim3"
        "; mv pbsim3/data/* {params.outdir}"
        "; rm -rf pbsim3"