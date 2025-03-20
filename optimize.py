import marimo

__generated_with = "0.11.22"
app = marimo.App(width="medium")


@app.cell
def _():
    from biotite.sequence.io.fasta import FastaFile

    from dnachisel import (
        AvoidPattern, 
        CodonOptimize,
        DnaOptimizationProblem,
        EnforceGCContent,
        EnforceTranslation,
    )

    from io import StringIO

    import marimo as mo

    import polars as pl

    from typing import Iterator
    return (
        AvoidPattern,
        CodonOptimize,
        DnaOptimizationProblem,
        EnforceGCContent,
        EnforceTranslation,
        FastaFile,
        Iterator,
        StringIO,
        mo,
        pl,
    )


@app.cell
def _(Iterator):
    aa_to_most_frequent_codons = dict(
        [
            ("*", "TGA"),
            ("A", "GCC"),
            ("C", "TGC"),
            ("D", "GAC"),
            ("E", "GAG"),
            ("F", "TTC"),
            ("G", "GGC"),
            ("H", "CAC"),
            ("I", "ATC"),
            ("K", "AAG"),
            ("L", "CTG"),
            ("M", "ATG"),
            ("N", "AAC"),
            ("P", "CCC"),
            ("Q", "CAG"),
            ("R", "AGA"),
            ("S", "AGC"),
            ("T", "ACC"),
            ("V", "GTG"),
            ("W", "TGG"),
            ("Y", "TAC"),
        ],
    )
    def backtranslate_orf(seq: str) -> Iterator[str]: 
        for aa in seq:
            if not aa in aa_to_most_frequent_codons:
                raise ValueError("Sequence contains invalid amino acids.")
            yield aa_to_most_frequent_codons[aa]
    return aa_to_most_frequent_codons, backtranslate_orf


@app.cell
def _(mo):
    sequences_input = mo.ui.text_area(
        full_width = True,
    )
    mo.vstack(
        [
            mo.md("### Sequences to optimize [FASTA format]"),
            sequences_input,
        ],
    )
    return (sequences_input,)


@app.cell
def _(mo):
    site_exclusion = mo.ui.multiselect(
        [
            "XbaI_site",
            "BamHI_site",
            "NotI_site",
            "NheI_site",
            "HindIII_site",
        ],
        full_width = True,
    )
    species = mo.ui.multiselect(
        [
            "b_subtilis",
            "c_elegans",
            "d_melanogaster",
            "e_coli",
            "g_gallus",
            "h_sapiens",
            "m_musculus",
            "m_musculus_domesticus",
            "s_cerevisiae",
        ],
        value = ["h_sapiens", "m_musculus"],
        full_width = True,
    )
    mo.hstack(
        [
            mo.vstack(
                [
                    mo.md("### Sites to exclude"),
                    site_exclusion,
                ],
            ),
            mo.vstack(
                [
                    mo.md("### Species"),
                    species,
                ],
            ),
        ],
        widths = [0.5, 0.5],
    )
    return site_exclusion, species


@app.cell
def _(
    AvoidPattern,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceGCContent,
    EnforceTranslation,
    FastaFile,
    StringIO,
    backtranslate_orf,
    mo,
    sequences_input,
    site_exclusion,
    species,
):
    sequences = FastaFile.read_iter(StringIO(sequences_input.value))
    site_constraints = None
    optimized_seqs = {}
    for header, seq in sequences:
        problem = DnaOptimizationProblem(
            sequence = "".join(backtranslate_orf(seq)),
            constraints = [
                *(
                    AvoidPattern(restriction_site, strand = 0)
                    for restriction_site in site_exclusion.value
                ),
                EnforceGCContent(mini = 0.2, maxi = 0.5, window = 100),
                EnforceTranslation(),
            ],
            objectives = [
                CodonOptimize(species = s)
                for s in species.value
            ]
        )
        problem.max_random_iters = 5000
        problem.resolve_constraints()
        problem.optimize()
        optimized_seqs[header] = problem.sequence

    mo.hstack(
        [
            mo.vstack(
                [
                    mo.ui.text(s, disabled = True)
                    for s in site_exclusion.value
                ],
            ),
            mo.vstack(
                [
                    mo.ui.text(s, disabled = True)
                    for s in species.value
                ],
            ),
        ],
        widths = [0.5, 0.5],
    )
    return header, optimized_seqs, problem, seq, sequences, site_constraints


@app.cell
def _(mo, pl):
    components = pl.read_csv(
        "component_db.csv",
        has_header = False,
        schema = {
            "name": pl.String,
            "value": pl.String,
            "length": pl.Int32,
        },
    )

    five_prime_flank = mo.ui.dropdown(
        components["name"].unique(),
        searchable = True,
        full_width = True,
    )
    three_prime_flank = mo.ui.dropdown(
        components["name"].unique(),
        searchable = True,
        full_width = True,
    )


    mo.hstack(
        [
            mo.vstack([mo.md("5' end"), five_prime_flank]),
            mo.vstack([mo.md("3' end"), three_prime_flank]),
        ],
        widths = [0.3, 0.3],
    )
    return components, five_prime_flank, three_prime_flank


@app.cell
def _(components, five_prime_flank, mo, optimized_seqs, pl, three_prime_flank):
    if five_prime_flank.value is None:
        five_prime_flank_value = ""
    else:
        five_prime_flank_value = components.filter(
            pl.col("name") == five_prime_flank.value
        )["value"].item()

    if three_prime_flank.value is None:
        three_prime_flank_value = ""
    else: 
        three_prime_flank_value = components.filter(
            pl.col("name") == three_prime_flank.value
        )["value"].item()

    final_seqs = {
        header: f"{five_prime_flank_value}{seq}{three_prime_flank_value}"
        for header, seq in optimized_seqs.items()
    }
    mo.vstack(
        [
            mo.md("### Optimized sequences"),
            mo.ui.table(
                pl.DataFrame({"name": final_seqs.keys(), "seq": final_seqs.values()})
            ),
        ],
    )
    return final_seqs, five_prime_flank_value, three_prime_flank_value


if __name__ == "__main__":
    app.run()
