import argparse
import random


def get_options():
    parser = argparse.ArgumentParser(
        description="Randomly transfer segments of tokens from genomes in a donor file "
                    "into the paired genomes of a recipient file (one direction only).",
        prog='python transfer_segments.py')
    parser.add_argument('--donor', required=True,
                        help='Donor file: one genome per line, "name<TAB>space-separated tokens"')
    parser.add_argument('--recipient', required=True,
                        help='Recipient file, same format. Genomes are paired with donor genomes by line order')
    parser.add_argument('--out', default='recipient_transferred.txt',
                        help='Output recipient file (same order as input). Default = "recipient_transferred.txt"')
    parser.add_argument('--summary', default='transfer_summary.tsv',
                        help='Output summary file. Default = "transfer_summary.tsv"')
    parser.add_argument('--n_transfers', type=int, default=1,
                        help='Number of transfers per genome pair. Default = 1')
    parser.add_argument('--min_length', type=int, default=1,
                        help='Minimum transfer length in tokens. Default = 1')
    parser.add_argument('--max_length', type=int, default=50,
                        help='Maximum transfer length in tokens. Default = 50')
    parser.add_argument('--seed', type=int, default=None, help='Random seed')
    return parser.parse_args()


def read_genomes(path):
    genomes = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            name, _, rest = line.partition('\t')
            genomes.append((name, rest.split()))
    return genomes


def contigs(tokens, sep="_"):
    """Return (start, end) half-open spans of tokens between contig break tokens."""
    spans = []
    start = None
    for i, tok in enumerate(tokens):
        if tok == sep:
            if start is not None:
                spans.append((start, i))
                start = None
        elif start is None:
            start = i
    if start is not None:
        spans.append((start, len(tokens)))
    return spans


def main():
    args = get_options()
    if args.min_length < 1 or args.max_length < args.min_length:
        raise ValueError("Require 1 <= min_length <= max_length")
    rng = random.Random(args.seed)

    donors = read_genomes(args.donor)
    recipients = read_genomes(args.recipient)
    if len(donors) != len(recipients):
        print(f"Warning: {len(donors)} donor and {len(recipients)} recipient genomes; "
              f"only the first {min(len(donors), len(recipients))} pairs are modified.")

    with open(args.out, 'w') as out, open(args.summary, 'w') as summ:
        summ.write("transfer\tdonor_genome\trecipient_genome\tlength\tdonor_start\tdonor_end\t"
                   "recipient_start\trecipient_end\tinserted_tokens\treplaced_tokens\n")
        for i, (rname, rtokens) in enumerate(recipients):
            if i < len(donors):
                dname, dtokens = donors[i]
                # Segments must lie entirely within a contig (no "_" break tokens),
                # both in the donor and in the recipient, and must not overlap.
                dcontigs = contigs(dtokens)
                rcontigs = contigs(rtokens)
                edits = []
                used = []
                for t in range(args.n_transfers):
                    for _ in range(100):
                        L = rng.randint(args.min_length, args.max_length)
                        dstarts = [s for cs, ce in dcontigs for s in range(cs, ce - L + 1)]
                        rstarts = [s for cs, ce in rcontigs for s in range(cs, ce - L + 1)
                                   if all(s + L <= us or s >= ue for us, ue in used)]
                        if dstarts and rstarts:
                            ds = rng.choice(dstarts)
                            rs = rng.choice(rstarts)
                            used.append((rs, rs + L))
                            edits.append((rs, L, ds))
                            break
                new = list(rtokens)
                records = []
                for rs, L, ds in sorted(edits, reverse=True):
                    ins = dtokens[ds:ds + L]
                    rep = new[rs:rs + L]
                    new[rs:rs + L] = ins
                    records.append((rs, L, ds, ins, rep))
                for n, (rs, L, ds, ins, rep) in enumerate(sorted(records), 1):
                    summ.write(f"{n}\t{dname}\t{rname}\t{L}\t{ds}\t{ds + L}\t{rs}\t{rs + L}\t"
                               f"{' '.join(ins)}\t{' '.join(rep)}\n")
                rtokens = new
            out.write(rname + '\t' + ' '.join(rtokens) + '\n')


if __name__ == "__main__":
    main()
