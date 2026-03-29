#!/usr/bin/env python

import os
from pathlib import Path
from typing import Optional

from tadashi.apps import Simple
from tadashi.translators import Polly, Translator

os.environ["LOG_LEVEL"] = "DEBUG"


class Advanture0(Simple):
    def __init__(
        self,
        source: str | Path = Path("adventure_kernel_region0_tune4_arm_pad-acle.c"),
        translator: Optional[Translator] = None,
        compiler_options: Optional[list[str]] = None,
        ephemeral: bool = False,
        populate_scops: bool = True,
        *,
        runtime_prefix: str = "WALLTIME: ",
    ):
        self.runtime_prefix = runtime_prefix
        super().__init__(
            source=source,
            translator=translator,
            compiler_options=compiler_options,
            ephemeral=ephemeral,
            populate_scops=populate_scops,
        )

    def app_required_options(self) -> list[str]:
        return ["-Icommon/include", "-fopenmp"]


def main():
    app = Advanture0(translator=Polly())
    print(f"{len(app.scops)=}")


if __name__ == "__main__":
    main()
