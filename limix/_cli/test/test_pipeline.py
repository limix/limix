import os

from numpy import nan
from numpy.random import RandomState
from numpy.testing import assert_allclose

import limix
from limix._cli.pipeline import Pipeline
from limix._cli.preprocess import impute, normalize, where
from limix.io._fetch import fetch


def test_pipeline_where_filter():

    filenames = [
        "chrom22_subsample20_maf0.10.bed",
        "chrom22_subsample20_maf0.10.fam",
        "chrom22_subsample20_maf0.10.bim",
    ]
    shapes = [
        (274, 5647),
        (274, 0),
        (274, 49008),
        (274, 49008),
        (4, 49008),
        (274, 5647),
        (274, 0),
        (274, 49008),
    ]
    specs = [
        "genotype:(16050612 <= pos) & (pos < 21050612)",
        "genotype:(16050612 <= pos) & (pos < 16050612)",
        "genotype",
        "genotype:",
        "genotype:sample.isin(['HG00111', 'HG00112', 'NA20775', 'NA20804'])",
        "genotype:(chrom == '22') & (16050612 <= pos) & (pos < 21050612)",
        "genotype: a0 == a1",
        "genotype: sample != None",
    ]
    random = RandomState(0)
    with limix.file_example(filenames) as filepaths:
        folder = os.path.dirname(filepaths[0])
        filepath = os.path.join(folder, "chrom22_subsample20_maf0.10")

        for shape, spec in zip(shapes, specs):
            G = fetch("genotype", f"{filepath}", verbose=False)
            y = random.randn(G.shape[0])
            data = {"G": G, "y": y}
            pipeline = Pipeline(data)
            pipeline.append(where, "where", spec)
            data = pipeline.run(verbose=False)
            assert data["G"].shape == shape
            assert data["G"].dims == ("sample", "candidate")


def test_pipeline_normalize():

    with limix.file_example("expr.csv") as filepath:
        folder = os.path.dirname(filepath)
        filepath = os.path.join(folder, "expr.csv")

        y = fetch("trait", f"{filepath}::row=trait", verbose=False)
        data = {"y": y}
        pipeline = Pipeline(data)
        pipeline.append(normalize, "normalize", spec="trait:trait:gaussianize")
        data = pipeline.run(verbose=False)
        assert_allclose(
            data["y"].values[0, :3],
            [
                -2.684_091_161_327_623_7,
                -0.348_755_695_517_044_7,
                -1.429_964_275_028_744_2,
            ],
        )

        y = fetch("trait", f"{filepath}::row=trait", verbose=False)
        data = {"y": y}
        pipeline = Pipeline(data)
        pipeline.append(normalize, "normalize", spec="trait:sample:gaussianize")
        data = pipeline.run(verbose=False)
        assert_allclose(
            data["y"].values[0, :3],
            [-1.382_994_127_100_638, 0.0, -0.674_489_750_196_081_7],
        )


def test_pipeline_impute():

    with limix.file_example("expr_nan.csv") as filepath:
        folder = os.path.dirname(filepath)
        filepath = os.path.join(folder, "expr_nan.csv")

        y = fetch("trait", f"{filepath}::row=trait", verbose=False)
        data = {"y": y}
        pipeline = Pipeline(data)
        pipeline.append(impute, "impute", spec="trait:trait:mean")
        data = pipeline.run(verbose=False)
        assert_allclose(
            data["y"].values[0, :3],
            [0.086_810_659_809_597_32, -0.351_453_183_202_999_95, -1.319_967_938_81],
        )


def test_pipeline_normalize_nan():
    with limix.file_example("expr_nan.csv") as filepath:
        folder = os.path.dirname(filepath)
        filepath = os.path.join(folder, "expr_nan.csv")

        y = fetch("trait", f"{filepath}::row=trait", verbose=False)
        data = {"y": y}
        pipeline = Pipeline(data)
        pipeline.append(normalize, "normalize", spec="trait:sample:gaussianize")
        data = pipeline.run(verbose=False)
        assert_allclose(data["y"].values[0, :3], [nan, 0.0, -0.841_621_233_572_914_3])

        y = fetch("trait", f"{filepath}::row=trait", verbose=False)
        data = {"y": y}
        pipeline = Pipeline(data)
        pipeline.append(normalize, "normalize", spec="trait:trait:gaussianize")
        data = pipeline.run(verbose=False)
        assert_allclose(
            data["y"].values[0, :3],
            [nan, -0.345_222_629_722_377_3, -1.429_964_275_028_744_2],
        )

        y = fetch("trait", f"{filepath}::row=trait", verbose=False)
        data = {"y": y}
        pipeline = Pipeline(data)
        pipeline.append(normalize, "normalize", spec="trait::gaussianize")
        data = pipeline.run(verbose=False)
        assert_allclose(
            data["y"].values[0, :3],
            [nan, -0.345_222_629_722_377_3, -1.429_964_275_028_744_2],
        )

        y = fetch("trait", f"{filepath}::row=trait", verbose=False)
        data = {"y": y}
        pipeline = Pipeline(data)
        pipeline.append(normalize, "normalize", spec="trait")
        data = pipeline.run(verbose=False)
        assert_allclose(
            data["y"].values[0, :3],
            [nan, -0.345_222_629_722_377_3, -1.429_964_275_028_744_2],
        )


def test_pipeline_impute_and_normalize():

    with limix.file_example("expr_nan.csv") as filepath:
        folder = os.path.dirname(filepath)
        filepath = os.path.join(folder, "expr_nan.csv")

        y = fetch("trait", f"{filepath}::row=trait", verbose=False)
        data = {"y": y}
        pipeline = Pipeline(data)
        pipeline.append(impute, "impute", spec="trait:trait:mean")
        pipeline.append(normalize, "normalize", spec="trait:trait:gaussianize")
        data = pipeline.run(verbose=False)
        assert_allclose(
            data["y"].values[0, :3],
            [
                -0.013_672_943_874_439_595,
                -0.348_755_695_517_044_7,
                -1.429_964_275_028_744_2,
            ],
        )

        y = fetch("trait", f"{filepath}::row=trait", verbose=False)
        data = {"y": y}
        pipeline = Pipeline(data)
        pipeline.append(impute, "impute", spec="trait::")
        pipeline.append(normalize, "normalize", spec="trait:")
        data = pipeline.run(verbose=False)
        assert_allclose(
            data["y"].values[0, :3],
            [
                -0.013_672_943_874_439_595,
                -0.348_755_695_517_044_7,
                -1.429_964_275_028_744_2,
            ],
        )
