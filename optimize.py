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

    from pathlib import Path

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
        Path,
        StringIO,
        mo,
        pl,
    )


@app.cell
def _(Iterator):
    def backtranslate_orf(seq: str) -> Iterator[str]: 
        aa_to_codon = dict(
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
        for aa in seq:
            if not aa in aa_to_codon:
                raise ValueError("Sequence contains invalid amino acids.")
            yield aa_to_codon[aa]
    return (backtranslate_orf,)


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
def _(Path, pl):
    CODON_TABLES_PATH = Path("./codon_usage_tables").resolve()
    def get_codon_table_species(species: str, convert_U_to_T: bool = True) -> dict[str, dict[str, float]]:
        codon_table = pl.read_csv(CODON_TABLES_PATH / f"{species}.csv")
        return {
            aa: {
                codon.replace("U", "T"): freq
                for codon, freq in codon_table.filter(
                    pl.col("amino_acid") == aa
                ).select(
                    pl.col("codon", "relative_frequency"),
                ).iter_rows()
            }
            for aa in codon_table["amino_acid"].unique()
        }

    def get_available_species() -> set[str]:
        return {
            s.stem for s in CODON_TABLES_PATH.glob("*.csv")
        }
    return CODON_TABLES_PATH, get_available_species, get_codon_table_species


@app.cell
def _(get_available_species, mo):
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
        get_available_species(),
        value = ["h_sapiens_9606"],
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
    get_codon_table_species,
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
                CodonOptimize(codon_usage_table = get_codon_table_species(s))
                for s in species.value
            ],
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


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
