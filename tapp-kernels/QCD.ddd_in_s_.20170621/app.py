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


class QcdDdd(App):
    def __init__(
        self,
        source: str | Path = Path("ddd_in_s.cc"),
        translator: Optional[Translator] = None,
        compiler_options: Optional[list[str]] = None,
        ephemeral: bool = False,
        populate_scops: bool = True,
        target: str = "DDD_IN_S",
    ):
        self.target = target
        super().__init__(
            source=source,
            translator=translator,
            compiler_options=compiler_options,
            ephemeral=ephemeral,
            populate_scops=populate_scops,
        )

    def codegen_init_args(self):
        return {"target": self.target}

    def app_required_options(self) -> list[str]:
        return [
            "-Icommon/include",
            "-fopenmp",
            "-DRDC",
            "-DVLENS=16",
            "-DEML_LIB",
            "-DPREFETCH",
            "-D_CHECK_SIM",
            "-DTARGET_JINV",
        ]

    def compile_cmd(self, suffix: str) -> list[str]:
        cmd = [
            "make",
            f"{str(self.output_binary)}",
            f"PROG={str(self.output_binary)}",
            f"{self.target}={str(self.source)}",
        ]
        return cmd

    def run_cmd(self):
        return [f"./{str(self.output_binary)}"]

    def extract_runtime(self, proc):
        print("{proc.stdout.decode()=")
        return 0.3


def main():
    kwargs = {"translator": "Polly"}
    run(QcdDdd, kwargs)


if __name__ == "__main__":
    main()
