import pysam
import argparse
import sys
from collections import defaultdict
import tempfile
from subprocess import Popen, PIPE
from multiprocessing import Pool


def iter_args(prot_fa, ref_seq):
    for rec in prot_fa:
        yield {"UniqProtID": rec.name, "ProtSeq": rec.sequence, "RefSeq": ref_seq}


def shell(cmd):
    p = Popen(cmd, shell=True, stdout=PIPE, close_fds=True)
    return p.stdout.readlines()
    
    
def find_mut(args):
    UniqProtID = args["UniqProtID"]
    ProtSeq = args["ProtSeq"]
    RefSeq = args["RefSeq"]
    ref_seq_align = ""
    seq_align = ""
    start_indx = 0
    with tempfile.NamedTemporaryFile() as f:
        f_tmp_file = f.name
        cmd = """
echo ">Ref" > {f_tmp}
echo "{ref_seq}" >> {f_tmp}
echo ">Seq" >> {f_tmp}
echo "{seq}" >> {f_tmp}
mafft --quiet --thread 1 {f_tmp}
        """.format(f_tmp=f_tmp_file, ref_seq=RefSeq, seq=ProtSeq)
        for line in shell(cmd):
            line = line.decode().rstrip("\n")
            if line.startswith(">"):
                start_indx += 1
            elif start_indx == 1:
                ref_seq_align += line
            elif start_indx == 2:
                seq_align += line
            else:
                raise ValueError()
    # Ref: TLLALH
    # Mut: T---LH
    # Mut: TL---H
    ins_len = ref_seq_align[:239].count("-")
    if seq_align[(239+ins_len): (245+ins_len)] == "T---LH":
        seq_align = seq_align[:(239+ins_len)] + "TL---H" + seq_align[(245+ins_len):]
    # Ref: HRSYLTPGDS
    # Mut: HN-------S
    # Mut: H-------NS
    ins_len = ref_seq_align[:244].count("-")
    if seq_align[(244+ins_len): (254+ins_len)] == "H-------NS":
        seq_align = seq_align[:(244+ins_len)] + "HN-------S" + seq_align[(254+ins_len):]
    # Ref: INLV
    # Mut: I-IV
    # Mut: II-V
    ins_len = ref_seq_align[:209].count("-")
    if seq_align[(209+ins_len):(213+ins_len)] == "I-IV":
        seq_align = seq_align[:(209+ins_len)] + "II-V" + seq_align[(213+ins_len):]

    pos = 0
    res_li = list()
    for (ref_aa, prot_aa) in zip(ref_seq_align, seq_align):
        if ref_aa == "-":
            continue
        pos = pos + 1
        if ref_aa == prot_aa:
            continue
        # if prot_aa == "X":
        #     continue
        res_li.append({"UniqProtID": UniqProtID, "Position": str(pos), "OriAa": ref_aa, "MutAa": prot_aa})
    return res_li
    
def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.fa", required=True)
    base_group.add_argument("--ref", type=str, dest="ref", metavar="ref.fa", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.tsv", required=True)
    base_group.add_argument("-p", "--process", type=int, dest="process", metavar="process", required=False, default=1)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_ref = args.ref
    f_out = args.output
    
    prot_fa = pysam.FastxFile(f_in)
    ref_fa = pysam.FastxFile(f_ref)
    ref_seq = next(ref_fa).sequence

    with open(f_out, "w") as f:
        header = ["UniqProtID", "Position", "OriAa", "MutAa"]
        f.write("\t".join(header) + "\n")
        with Pool(args.process) as p:
            for res_li in p.imap_unordered(find_mut, iter_args(prot_fa, ref_seq)):
                for res in res_li:
                    f.write("\t".join([res[x] for x in header])+"\n")


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
