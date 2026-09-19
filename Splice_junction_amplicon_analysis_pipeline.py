#!/usr/bin/env python3
"""Splice junction amplicon analysis pipeline, filter, align, classify, plot.

Everything the run needs comes from two places, an annotated reference of the
full amplicon and one or more FASTQ files.

  1. Filter. A read is kept when it carries both ends of the reference, the
     first and last --anchor bases, in either orientation and within
     --mismatches of a perfect match. Kept reads are written out in reference
     orientation. Reads missing either end are discarded.
  2. Align. Kept reads are aligned to the reference with minimap2 and the
     alignments are coordinate sorted and indexed with samtools.
  3. Classify. The variable interval is taken from the reference annotation,
     the feature whose type or name matches --feature. A read is spliced when
     its alignment deletes at least --deletion-fraction of that interval,
     unspliced when the alignment spans the interval without such a deletion,
     and other in every remaining case.
  4. Plot. One plot per --reads group, a stacked bar per library within that
     group, spliced, unspliced and other as a percentage of the aligned reads.
     Libraries from different source folders never share a plot.

The reference may be SnapGene (.dna) or GenBank (.gb, .gbk), which carry the
annotation, or FASTA (.fa, .fasta) together with --interval START-END.

Example
  python3 Splice_junction_amplicon_analysis_pipeline.py \
      --reference "../RT-qPCR assays/mScarlet_REF.dna" \
      --reads LPXHYS_RTassay1_cDNA 7SY36N_gel_purified_bands \
      --outdir splice_junction_out \\
      --labels 1 2 3 heterodimer unspliced spliced \\
      --xlabel "Biological replicate" ""
"""
import argparse, csv, glob, json, os, re, struct, subprocess, sys, xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------- reference

def read_snapgene(path):
    """Return (sequence, [(name, type, start, end)]) with 1-based inclusive coords."""
    raw = open(path, "rb").read()
    seq, feats, i = "", [], 0
    while i < len(raw):
        kind = raw[i]
        size = struct.unpack(">I", raw[i + 1:i + 5])[0]
        body = raw[i + 5:i + 5 + size]
        if kind == 0:
            seq = body[1:].decode("latin1").upper()
        elif kind == 10:
            root = ET.fromstring(body.decode("utf8", "replace"))
            for f in root.iter("Feature"):
                spans = [s.get("range") for s in f.findall("Segment") if s.get("range")]
                if not spans:
                    continue
                lo = min(int(r.split("-")[0]) for r in spans)
                hi = max(int(r.split("-")[1]) for r in spans)
                feats.append((f.get("name") or "", f.get("type") or "", lo, hi))
        i += 5 + size
    return seq, feats


def read_genbank(path):
    text = open(path).read()
    feats = []
    for m in re.finditer(r"^ {5}(\S+)\s+(?:complement\()?<?(\d+)\.\.>?(\d+)", text, re.M):
        kind, lo, hi = m.group(1), int(m.group(2)), int(m.group(3))
        tail = text[m.end():m.end() + 400]
        lab = re.search(r'/(?:label|gene|product|note)="([^"]+)"', tail)
        feats.append((lab.group(1) if lab else kind, kind, lo, hi))
    origin = text.split("ORIGIN", 1)[1] if "ORIGIN" in text else ""
    seq = re.sub(r"[^A-Za-z]", "", origin).upper()
    return seq, feats


def read_fasta(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">"))
    return seq.upper(), []


def load_reference(path, feature, interval):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".dna":
        seq, feats = read_snapgene(path)
    elif ext in (".gb", ".gbk", ".genbank"):
        seq, feats = read_genbank(path)
    elif ext in (".fa", ".fasta", ".fna"):
        seq, feats = read_fasta(path)
    else:
        sys.exit(f"unrecognized reference format, {path}")
    if not seq:
        sys.exit(f"no sequence found in {path}")
    if interval:
        lo, hi = (int(x) for x in interval.replace("..", "-").split("-"))
    else:
        hits = [f for f in feats if feature.lower() in f[1].lower()] or \
               [f for f in feats if feature.lower() in f[0].lower()]
        if not hits:
            sys.exit(f"no feature matching '{feature}' in {path}, pass --interval instead")
        if len({(f[2], f[3]) for f in hits}) > 1:
            sys.exit("several features match --feature, narrow it or pass --interval")
        lo, hi = hits[0][2], hits[0][3]
        print(f"variable interval from reference annotation, {hits[0][0]} at {lo}-{hi}")
    return seq, lo, hi


# ------------------------------------------------------------------- filter

def within(hay, needle, mismatches, start=0):
    """Leftmost index at or after start where needle matches hay with <= mismatches."""
    n, h = len(needle), len(hay)
    hit = hay.find(needle, start)
    if hit >= 0:
        return hit
    if mismatches <= 0:
        return -1
    for i in range(start, h - n + 1):
        bad = 0
        for a, b in zip(hay[i:i + n], needle):
            if a != b:
                bad += 1
                if bad > mismatches:
                    break
        else:
            return i
    return -1


def oriented(seq, qual, head, tail, mismatches):
    """Return (seq, qual) in reference orientation when both anchors are present."""
    for s, q in ((seq, qual), (rc(seq), qual[::-1])):
        i = within(s, head, mismatches)
        if i < 0:
            continue
        j = within(s, tail, mismatches, i + len(head))
        if j < 0:
            continue
        return s, q
    return None


def filter_fastq(path, head, tail, mismatches, out_path):
    total = kept = 0
    with open(path) as fh, open(out_path, "w") as out:
        while True:
            name = fh.readline()
            if not name:
                break
            seq, plus, qual = fh.readline().strip(), fh.readline(), fh.readline().strip()
            total += 1
            got = oriented(seq.upper(), qual, head, tail, mismatches)
            if got:
                kept += 1
                out.write(f"{name.rstrip()}\n{got[0]}\n+\n{got[1]}\n")
    return total, kept


# -------------------------------------------------------------------- align

def align(fastq, ref_fa, bam, threads):
    with open(bam, "wb") as out:
        mm = subprocess.Popen(["minimap2", "-ax", "map-ont", "--MD", "-t", str(threads),
                               ref_fa, fastq], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        st = subprocess.Popen(["samtools", "sort", "-o", "-"], stdin=mm.stdout, stdout=out,
                              stderr=subprocess.DEVNULL)
        mm.stdout.close()
        st.communicate()
    subprocess.run(["samtools", "index", bam], check=True)


# ----------------------------------------------------------------- classify

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def classify(bam, lo, hi, fraction):
    """lo, hi are 1-based inclusive reference coordinates of the variable interval."""
    need = fraction * (hi - lo + 1)
    calls, per_read = {"spliced": 0, "unspliced": 0, "other": 0}, {}
    out = subprocess.run(["samtools", "view", "-F", "0x904", bam],
                         capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        f = line.split("\t")
        read, pos, cig = f[0], int(f[3]), f[5]
        ref = pos
        gap = 0
        for n, op in CIGAR.findall(cig):
            n = int(n)
            if op in "DN":
                gap += max(0, min(ref + n - 1, hi) - max(ref, lo) + 1)
                ref += n
            elif op in "M=X":
                ref += n
        end = ref - 1
        if gap >= need:
            call = "spliced"
        elif pos <= lo and end >= hi:
            call = "unspliced"
        else:
            call = "other"
        calls[call] += 1
        per_read[read] = (call, gap)
    return calls, per_read



# --------------------------------------------------------------------- plot

CLASS_COLOUR = {"spliced": "#56adde", "unspliced": "#0f6fb0", "other": "#9e9e9e"}
LABEL_COLOUR = "#ffffff"


def plot_composition(rows, labels, out_base, xlabel, geom):
    """geom carries the panel geometry, all in points, so a panel can be dropped
    straight into an existing figure with the same axes box as the panel it
    replaces."""
    fs = geom["fontsize"]
    plt.rcParams.update({"font.size": fs, "font.family": "sans-serif",
                         "font.sans-serif": ["Arial"], "axes.linewidth": 1.4,
                         "xtick.labelsize": fs, "ytick.labelsize": fs,
                         "xtick.direction": "out", "ytick.direction": "out",
                         "xtick.major.size": 3.5, "ytick.major.size": 3.5,
                         "xtick.major.width": 1.4, "ytick.major.width": 1.4,
                         "figure.facecolor": "white", "axes.facecolor": "white",
                         "savefig.facecolor": "white", "legend.frameon": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    n = len(rows)
    if geom["size"]:
        w_pt, h_pt = geom["size"]
        left, right, top, bottom = geom["axes"]
        fig = plt.figure(figsize=(w_pt / 72, h_pt / 72))
        ax = fig.add_axes([left / w_pt, (h_pt - bottom) / h_pt,
                           (right - left) / w_pt, (bottom - top) / h_pt])
        ax.set_ylim(0, geom["ytop"])
        half = geom["xpad"]
        ax.set_xlim(-half, n - 1 + half)
        width = geom["barwidth"]
    else:
        fig = plt.figure(figsize=(1.15 * n + 3.6, 4.6))
        ax = fig.add_axes([0.62 / (1.15 * n + 3.6) + 0.06, 0.26,
                           1.05 * n / (1.15 * n + 3.6), 0.62])
        ax.set_ylim(0, 100)
        ax.set_xlim(-0.62, n - 1 + 0.62)
        width = 0.60

    x = list(range(n))
    bottom_vals = [0.0] * n
    for cls in ("spliced", "unspliced", "other"):
        vals = [r[cls + "_pct"] for r in rows]
        ax.bar(x, vals, width, bottom=bottom_vals, color=CLASS_COLOUR[cls],
               edgecolor="white", lw=1.0, label=cls.capitalize())
        for i, (b, v) in enumerate(zip(bottom_vals, vals)):
            if cls == "other":
                ax.text(i, 101.5, f"{v:.0f}", ha="center", va="bottom",
                        fontsize=fs, color="#000000")
            elif v > 6:
                ax.text(i, b + v / 2, f"{v:.0f}", ha="center", va="center",
                        fontsize=fs, color=LABEL_COLOUR)
        bottom_vals = [b + v for b, v in zip(bottom_vals, vals)]

    ax.set_xticks(x)
    rot = geom["rotate"]
    if rot is None:
        rot = 30 if any(len(str(l)) > 4 for l in labels) else 0
    ax.set_xticklabels(labels, rotation=rot, ha="right" if rot else "center")
    ax.set_yticks([0, 50, 100])
    ax.set_ylabel("Percent of reads")
    ax.set_clip_on(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if xlabel:
        ax.set_xlabel(xlabel)
    if geom["legend"]:
        h, l = ax.get_legend_handles_labels()
        ax.legend(h, l, loc="upper left", bbox_to_anchor=(1.04, 1.0), handlelength=1.1,
                  handleheight=1.1, borderpad=0, labelspacing=0.9, fontsize=fs)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(out_base + "." + ext, dpi=300)
    plt.close(fig)


# --------------------------------------------------------------------- main

def read_groups(paths):
    """One group per --reads argument, so separate sources stay on separate plots."""
    groups = []
    for p in paths:
        if os.path.isdir(p):
            fq = sorted(glob.glob(os.path.join(p, "*.fastq"))) + \
                 sorted(glob.glob(os.path.join(p, "*.fq")))
            groups.append((os.path.basename(os.path.normpath(p)), fq))
        else:
            groups.append((os.path.splitext(os.path.basename(p))[0], [p]))
    return groups


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--reads", required=True, nargs="+",
                    help="FASTQ files, or folders holding them")
    ap.add_argument("--outdir", default="amplicon_pipeline_out")
    ap.add_argument("--feature", default="intron",
                    help="reference feature giving the variable interval")
    ap.add_argument("--interval", help="START-END on the reference, for an unannotated FASTA")
    ap.add_argument("--anchor", type=int, default=20,
                    help="bases at each end of the reference used as the primer anchors")
    ap.add_argument("--mismatches", type=int, default=2)
    ap.add_argument("--deletion-fraction", type=float, default=0.5)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--labels", nargs="+",
                    help="short x tick label per library, in the order the reads are given")
    ap.add_argument("--xlabel", nargs="+", default=None,
                    help="x axis label, one per --reads group")
    ap.add_argument("--panel-size", nargs=2, type=float, metavar=("W", "H"),
                    help="panel width and height in points, to match a figure panel")
    ap.add_argument("--panel-axes", nargs=4, type=float,
                    metavar=("LEFT", "RIGHT", "TOP", "BOTTOM"),
                    help="axes box in points, measured from the top left of the panel")
    ap.add_argument("--fontsize", type=float, default=20)
    ap.add_argument("--rotate", type=float, default=None,
                    help="x tick label rotation, default picks 0 or 30 by label length")
    ap.add_argument("--bar-width", type=float, default=0.645)
    ap.add_argument("--x-pad", type=float, default=0.6774,
                    help="x axis padding either side of the first and last bar")
    ap.add_argument("--y-top", type=float, default=112.0)
    ap.add_argument("--no-legend", action="store_true")
    a = ap.parse_args()

    seq, lo, hi = load_reference(a.reference, a.feature, a.interval)
    head, tail = seq[:a.anchor], seq[-a.anchor:]
    name = os.path.splitext(os.path.basename(a.reference))[0]

    for sub in ("reference", "filtered", "bam"):
        os.makedirs(os.path.join(a.outdir, sub), exist_ok=True)
    ref_fa = os.path.join(a.outdir, "reference", name + ".fa")
    with open(ref_fa, "w") as fh:
        fh.write(f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")

    print(f"reference {name}, {len(seq)} nt, variable interval {lo}-{hi} ({hi-lo+1} nt)")
    print(f"anchors {head} and {tail}, up to {a.mismatches} mismatches each\n")

    if bool(a.panel_size) != bool(a.panel_axes):
        sys.exit("--panel-size and --panel-axes go together")
    geom = {"size": a.panel_size, "axes": a.panel_axes, "fontsize": a.fontsize,
            "rotate": a.rotate, "barwidth": a.bar_width, "xpad": a.x_pad,
            "ytop": a.y_top, "legend": not a.no_legend}

    groups = read_groups(a.reads)
    if a.labels and len(a.labels) != sum(len(g[1]) for g in groups):
        sys.exit(f"{len(a.labels)} labels given for "
                 f"{sum(len(g[1]) for g in groups)} libraries")
    if a.xlabel and len(a.xlabel) != len(groups):
        sys.exit(f"{len(a.xlabel)} x axis labels given for {len(groups)} read groups")

    rows, classes, group_rows, taken = [], {}, [], 0
    for gname, fastqs in groups:
      these = []
      for fq in fastqs:
        lib = os.path.splitext(os.path.basename(fq))[0]
        kept_fq = os.path.join(a.outdir, "filtered", lib + ".fastq")
        total, kept = filter_fastq(fq, head, tail, a.mismatches, kept_fq)
        bam = os.path.join(a.outdir, "bam", lib + ".bam")
        align(kept_fq, ref_fa, bam, a.threads)
        calls, per_read = classify(bam, lo, hi, a.deletion_fraction)
        classes[lib] = per_read
        aligned = sum(calls.values())
        pct = lambda n: round(100 * n / aligned, 1) if aligned else 0.0
        row = {"group": gname, "library": lib, "reads": total, "passed_filter": kept,
               "aligned": aligned, "spliced": calls["spliced"],
               "unspliced": calls["unspliced"], "other": calls["other"],
               "spliced_pct": pct(calls["spliced"]),
               "unspliced_pct": pct(calls["unspliced"]),
               "other_pct": pct(calls["other"])}
        rows.append(row)
        these.append(row)
        print(f"{lib:52s} reads {total:6d}  filter {kept:6d}  "
              f"spliced {pct(calls['spliced']):5.1f}  unspliced {pct(calls['unspliced']):5.1f}  "
              f"other {pct(calls['other']):5.1f}")
      group_rows.append((gname, these))

    with open(os.path.join(a.outdir, "classification.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(a.outdir, "read_classes.json"), "w") as fh:
        json.dump({k: {r: v[0] for r, v in d.items()} for k, d in classes.items()}, fh)
    print(f"\nwrote {os.path.join(a.outdir, 'classification.csv')}")
    for i, (gname, these) in enumerate(group_rows):
        if a.labels:
            labels = a.labels[taken:taken + len(these)]
            taken += len(these)
        else:
            labels = [r["library"] for r in these]
        out_base = os.path.join(a.outdir, "composition_" + gname)
        plot_composition(these, labels, out_base, a.xlabel[i] if a.xlabel else "", geom)
        print(f"wrote {out_base}.pdf and {out_base}.png")


if __name__ == "__main__":
    main()
