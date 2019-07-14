from loguru import logger
from typing import List


class QTLTask:
    def run(self):
        raise NotImplementedError

    def _write_toml(self, outdir, params):
        import qtoml

        with open(outdir / "task.toml", "w") as f:
            qtoml.dump(params, f)


class QTLSTTask(QTLTask):
    def __init__(self, trait: str):
        self._trait = trait

    def run(self, input, task_id: int, verbose: bool):
        import limix
        from limix._data import conform_dataset

        logger.info(f"Single-trait run for {self._trait}.")
        y = input.traits.sel(trait=self._trait)
        covariates = input.covariates
        kinship = input.kinship
        genotype = input.genotype

        data = {"y": y, "M": covariates, "K": kinship, "G": genotype}
        data = conform_dataset(**data)
        data = {k: v for k, v in data.items() if v is not None}

        if "K" not in data:
            data["K"] = None

        results = limix.qtl.scan(
            data["G"],
            data["y"],
            lik="normal",
            K=data["K"],
            M=data["M"],
            verbose=verbose,
        )

        outdir = input.outdir / f"task{task_id}"
        logger.info(f"Output directory: {outdir}")
        if not outdir.exists():
            outdir.mkdir()

        self._write_results(results, outdir, verbose)
        self._write_toml(outdir)

    def _write_results(self, results, outdir, verbose: bool):
        from limix._display import session_line

        outdir_str = str(outdir.absolute())
        msg = f"Saving results to folder <{outdir_str}>... "
        with session_line(msg, disable=not verbose):
            results.to_csv(
                outdir / "h0_effsizes.csv",
                outdir / "h0_variances.csv",
                outdir / "h2_effsizes.csv",
                outdir / "stats.csv",
            )

    def _write_toml(self, outdir):
        super()._write_toml(outdir, {"method": "st", "trait": self._trait})


class QTLMTTask(QTLTask):
    def __init__(self, traits: List[str]):
        self._traits = traits
