import os

import limix
from limix._cli.pipeline import Pipeline
from limix._cli.preprocess import where_filter


def test_pipeline():

    filenames = [
        "chrom22_subsample20_maf0.10.bed",
        "chrom22_subsample20_maf0.10.fam",
        "chrom22_subsample20_maf0.10.bim",
    ]
    shapes = [
        (274, 5647),
        (274, 0),
        (274, 49008),
        (4, 49008),
        (274, 5647),
        (274, 0),
        (274, 49008),
    ]
    exprs = [
        "(16050612 <= pos) & (pos < 21050612)",
        "(16050612 <= pos) & (pos < 16050612)",
        "",
        "sample.isin(['HG00111', 'HG00112', 'NA20775', 'NA20804'])",
        "(chrom == '22') & (16050612 <= pos) & (pos < 21050612)",
        "a0 == a1",
        "sample != None",
    ]
    with limix.file_example(filenames) as filepaths:
        folder = os.path.dirname(filepaths[0])
        filepath = os.path.join(folder, "chrom22_subsample20_maf0.10")

        for shape, expr in zip(shapes, exprs):
            G = limix.io.fetch("genotype", f"{filepath}", verbose=False)
            data = {"G": G}
            pipeline = Pipeline(data)
            pipeline.append(where_filter, expr=expr, target="genotype")
            data = pipeline.run(verbose=False)
            assert data["G"].shape == shape
            assert data["G"].dims == ("sample", "candidate")
