#!/usr/bin/env python

import os
import sys
from pathlib import Path
from typing import Optional

from tadashi.apps import App
from tadashi.translators import Polly, Translator

ml4tadashi = os.path.dirname(__file__)
ml4tadashi = os.path.dirname(ml4tadashi)
ml4tadashi = os.path.dirname(ml4tadashi)
ml4tadashi = os.path.dirname(ml4tadashi)
ml4tadashi = os.path.dirname(ml4tadashi)
ml4tadashi = os.path.dirname(ml4tadashi)
sys.path.append(ml4tadashi)
ml4tadashi = os.path.join(ml4tadashi, "ML4TADASHI")
sys.path.append(ml4tadashi)

from ML4TADASHI import run


class NicamViRhowSolver(App):
    def __init__(
        self,
        source: str | Path = Path("postK_nicam_vi_rhow_solver.f90"),
        translator: Optional[Translator] = None,
        compiler_options: Optional[list[str]] = None,
        ephemeral: bool = False,
        populate_scops: bool = True,
    ):
        super().__init__(
            source=source,
            translator=translator,
            compiler_options=compiler_options,
            ephemeral=ephemeral,
            populate_scops=populate_scops,
        )

    def codegen_init_args(self):
        return []

    def app_required_options(self) -> list[str]:
        return ["-Icommon/include"]

    def compile_cmd(self, suffix: str) -> list[str]:
        self.source.touch()
        cmd = ["make", f"{str(self.output_binary)}.x", f"APP={str(self.output_binary)}"]
        return cmd

    def run_cmd(self):
        return [f"./{str(self.output_binary)}.x"]

    def extract_runtime(self, proc):
        stdout = proc.stdout.decode()
        for line in stdout.split("\n"):
            if "WALLTIME:" in line:
                rv = float(line.replace("s", "").split()[1])
                return rv
        raise RuntimeError()


def main():
    kwargs = {
            "translator": "Polly",
            "translator_params": "flang",
            }
    # app = NicamViRhowSolver.mkapp(kwargs)
    run(NicamViRhowSolver, kwargs)


if __name__ == "__main__":
    main()


