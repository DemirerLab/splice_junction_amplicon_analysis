# Splice Junction Amplicon Analysis

Filters, aligns and classifies amplicon sequencing reads by the presence or absence of a variable interval, such as a spliced intron, and plots the composition of each library.

Demirer Lab, Caltech Division of Chemistry and Chemical Engineering (CCE)

---

## What it does

| Step | Description |
|------|-------------|
| 1 | **Filter.** A read is kept when it carries both ends of the reference, the first and last `--anchor` bases, in either orientation and within `--mismatches` of a perfect match. Since the reference is the full amplicon, its ends are the primer sites, so no separate primer file is needed. Kept reads are written out in reference orientation. |
| 2 | **Align.** Kept reads are aligned to the reference with minimap2 and the alignments are coordinate sorted and indexed with samtools. |
| 3 | **Classify.** The variable interval is read from the reference annotation. A read is spliced when its alignment deletes at least `--deletion-fraction` of that interval, unspliced when the alignment spans the interval without such a deletion, and other in every remaining case. |
| 4 | **Plot.** One stacked bar chart per read group, spliced, unspliced and other as a percentage of the aligned reads. Libraries from different source folders stay on separate plots. |

The deletion fraction is what groups minor splice variants with the main spliced product. At the default of 0.5, a deletion anywhere from half the interval to its full length counts as spliced.

---

## Dependencies

```bash
pip install matplotlib
```

Python ≥ 3.8, plus [minimap2](https://github.com/lh3/minimap2) and [samtools](https://www.htslib.org) on the PATH.

Developed against minimap2 2.31 and samtools 1.24.

---

## Usage

```bash
python Splice_junction_amplicon_analysis_pipeline.py \
    --reference reference.dna \
    --reads reads_folder \
    --outdir results
```

The reference may be SnapGene (`.dna`) or GenBank (`.gb`, `.gbk`), which carry the annotation, or FASTA (`.fa`, `.fasta`) together with `--interval START-END`. Each `--reads` argument, a folder or a single FASTQ, becomes its own group and its own plot.

### Options

| Option | Default | Description |
|---|---|---|
| `--feature` | `intron` | Reference feature giving the variable interval |
| `--interval` | | `START-END` on the reference, for an unannotated FASTA |
| `--anchor` | `20` | Bases at each end of the reference used as the primer anchors |
| `--mismatches` | `2` | Mismatches allowed per anchor |
| `--deletion-fraction` | `0.5` | Fraction of the interval a deletion must cover to count as spliced |
| `--threads` | `8` | Threads passed to minimap2 |
| `--labels` | | One x tick label per library, in the order the reads are given |
| `--xlabel` | | One x axis label per read group |
| `--no-legend` | | Omit the colour key, for a panel that shares one |
| `--panel-size` | | Panel width and height in points, to match an existing figure panel |
| `--panel-axes` | | Axes box in points, measured from the top left of the panel |
| `--bar-width`, `--x-pad`, `--y-top`, `--fontsize`, `--rotate` | | Panel geometry and type size |

Passing `--panel-size` and `--panel-axes` makes the plot at an exact size with an exact axes box, so it can be dropped straight into a figure in place of an existing panel.

---

## Outputs

| Path | Contents |
|---|---|
| `reference/` | The reference written out as FASTA |
| `filtered/` | Reads that carried both anchors, in reference orientation |
| `bam/` | Indexed alignments of the filtered reads |
| `classification.csv` | Per library, reads, reads passing the filter, reads aligned, and the counts and percentages in each class |
| `read_classes.json` | Every read with the class it was called |
| `composition_<group>.pdf/.png/.svg` | The stacked bar chart for each read group |

---

## License

MIT — see [LICENSE](LICENSE)
