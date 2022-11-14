import os

rule parse_prot_fa:
    input:
        prot_fa = os.path.join("data", "prot", "spikeprot0705","spikeprot0705.fasta")
    output:
        seq_info = os.path.join("analysis", "seq_info", "seq_info.tsv"),
        uniq_fa = os.path.join("analysis", "seq_info", "uniq_prot.fa"),
        pass_fa = os.path.join("analysis", "seq_info", "PassProt.fa"),
        fail_fa = os.path.join("analysis", "seq_info", "FailProt.fa"),
        seq2uniq = os.path.join("analysis", "seq_info", "seq2uniq.tsv"),
    params:
        parse_prot_fa = "python script/parse_prot_fa.py",
        min_aa_num = 1173,
    shell:
        """
{params.parse_prot_fa} -i {input.prot_fa} --out-info {output.seq_info} --out-uniq-fa {output.uniq_fa} --out-pass-fa {output.pass_fa} --out-fail-fa {output.fail_fa} --out-seq2uniq {output.seq2uniq} --min-aa-num {params.min_aa_num}
        """
        
rule prot2mut:
    input:
        uniq_fa = rules.parse_prot_fa.output.uniq_fa,
        ref_fa = os.path.join("data", "prot", "spike.ref.fasta"),
    output:
        mut = os.path.join("analysis", "seq_info", "uniq_prot.mut.tsv"),
    params:
        prot2mut = "python script/prot2mut.py",
    shell:
        """
{params.prot2mut} -i {input.uniq_fa} --ref {input.ref_fa} -o {output.mut} -p 100
        """
        
rule find_ann_mut_type:
    input:
        mut = rules.prot2mut.output.mut,
        mut_info = os.path.join("data", "MutInfo.xlsx"),
        ref_fa = os.path.join("data", "prot", "spike.ref.fasta"),
    output:
        mut_type = os.path.join("analysis", "seq_info", "uniq_prot.ann_mut_type.tsv"),
    params:
        mut2type = "Rscript script/mut2type.R",
    shell:
        """
{params.mut2type} --mut-pos-info {input.mut_info} --ref-fa {input.ref_fa} --prot-mut {input.mut} -o {output.mut_type}
        """

rule stat_variants:
    input:
        variant_info = os.path.join("data", "VariantInfo.tsv"),
        ref_fa = os.path.join("data", "prot", "spike.ref.fasta"),
        mut = rules.prot2mut.output.mut,
        seq_info = rules.parse_prot_fa.output.seq_info,
        seq2uniq = rules.parse_prot_fa.output.seq2uniq,
    output:
        mut_time = os.path.join("analysis", "variant_stat", "MutTimeNum.tsv"),
    params:
        stat_variants = "Rscript script/stat_variants.R",
        out_dir = os.path.join("analysis", "variant_stat"),
    shell:
        """
{params.stat_variants} --variant-info {input.variant_info} --ref-fa {input.ref_fa} --prot-mut {input.mut} --seq-info {input.seq_info} --seq2uniq {input.seq2uniq} -o {params.out_dir}
        """

rule all:
    input:
        rules.stat_variants.output,
