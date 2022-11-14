import pysam
import argparse
import sys
from collections import defaultdict


def compute_aa_num(seq):
    return len(seq) - seq.count("X")


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.tsv", required=True)
    base_group.add_argument("--out-info", type=str, dest="out_info", metavar="SeqInfo.tsv", required=True)
    base_group.add_argument("--out-seq2uniq", type=str, dest="out_seq2uniq", metavar="seq2uniq.tsv", required=True)
    base_group.add_argument("--out-uniq-fa", type=str, dest="out_uniq_fa", metavar="UniqFa.fa", required=True)
    base_group.add_argument("--out-pass-fa", type=str, dest="out_pass_fa", metavar="PassFa.fa", required=True)
    base_group.add_argument("--out-fail-fa", type=str, dest="out_fail_fa", metavar="FailFa.fa", required=True)
    base_group.add_argument("--min-aa-num", type=int, dest="min_aa_num", metavar="min_aa_num", required=True)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out_info = args.out_info
    min_aa_num = args.min_aa_num
    f_out_seq2uniq = open(args.out_seq2uniq, "w")
    f_out_uniq_fa = open(args.out_uniq_fa, "w")
    f_out_pass_fa = open(args.out_pass_fa, "w")
    f_out_fail_fa = open(args.out_fail_fa, "w")

    uniq_prot_dict = defaultdict(list)
    prof_fa = pysam.FastxFile(f_in)

    with open(f_out_info, "w") as f:
        header = "GeneName\tIsolateName\tYear\tMonth\tDay\tIsolateID\tPassageDetails\tType\tState\tHost\tOriginatingLab\tSubmittingLab\tSubmitter\tLocation\n"
        f.write(header)
        for rec in prof_fa:
            name = rec.name
            if rec.comment is not None:
                name = name + " " + rec.comment
            name = name.rstrip(" ")
            prot_seq = rec.sequence.upper().rstrip("*")
            
            try:
                gene_name, iso_name, time, iso_id, detail, type_state, host, ori_lab, submit_lab, Submitter, Location = name.split("|")
                year, month, day = time.split("-")
                virus_type, state = type_state.split("^^")
                f_out_pass_fa.write(">{0}\n{1}\n".format(name, prot_seq))
            except Exception:
                f_out_fail_fa.write(">{0}\n{1}\n".format(name, prot_seq))
            aa_num = compute_aa_num(prot_seq)
            if aa_num < min_aa_num:
                continue
            uniq_prot_dict[prot_seq].append(iso_name)
            data = [gene_name, iso_name, year, month, day, iso_id, detail, virus_type, state, host, ori_lab, submit_lab, Submitter, Location]
            f.write("\t".join(list(map(str, data)))+"\n")
    
    f_out_seq2uniq.write("IsolateName\tUniqProtID\n")
    for indx, (prot_seq, iso_name_li) in enumerate(uniq_prot_dict.items()):
        f_out_uniq_fa.write(">{0}\n{1}\n".format(indx, prot_seq))
        for iso_name in iso_name_li:
            f_out_seq2uniq.write("{0}\t{1}\n".format(iso_name, indx))
    
    f_out_pass_fa.close()
    f_out_fail_fa.close()
    f_out_seq2uniq.close()
    f_out_uniq_fa.close()
            

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
