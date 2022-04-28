import datetime
import logging
import os
from dataclasses import dataclass
from typing import Optional

import pkg_resources


PIPOLB_HMM_PROFILE = pkg_resources.resource_filename(
        'explore_pipolin',
        'data/pipolb_expanded_definitive.hmm',
    )

REF_ATT = pkg_resources.resource_filename('explore_pipolin', 'data/attL.fa')

PROTEINS = pkg_resources.resource_filename('explore_pipolin', '/data/HHpred_proteins.faa')

_BORDER_INFLATE = 0
_NO_BORDER_INFLATE = 30_000


@dataclass()
class GlobalSettings:
    out_dir: str
    pipolb_hmm_profile: str
    ref_att: str
    percent_identity: int
    border_inflate = _BORDER_INFLATE
    no_border_inflate: int
    proteins: str
    prokka_cpus: int
    skip_colours: bool

    @staticmethod
    def create_instance(
            out_dir_prefix,
            user_defined_out_dir,
            user_defined_profile,
            user_defined_att,
            percent_identity,
            max_inflate,
            user_defined_proteins,
            prokka_cpus,
            skip_colours,
    ):

        profile = user_defined_profile if user_defined_profile is not None else PIPOLB_HMM_PROFILE
        att = user_defined_att if user_defined_att is not None else REF_ATT
        proteins = user_defined_proteins if user_defined_proteins is not None else PROTEINS

        return GlobalSettings(
            out_dir=_get_out_dir(out_dir_prefix, user_defined_out_dir),
            pipolb_hmm_profile=profile,
            ref_att=att,
            percent_identity=percent_identity,
            no_border_inflate=max_inflate,
            proteins=proteins,
            prokka_cpus=prokka_cpus,
            skip_colours=skip_colours
        )


_GLOBAL_SETTINGS_INSTANCE: Optional[GlobalSettings] = None


def set_instance(settings: GlobalSettings):
    global _GLOBAL_SETTINGS_INSTANCE
    _GLOBAL_SETTINGS_INSTANCE = settings


def get_instance() -> GlobalSettings:
    if _GLOBAL_SETTINGS_INSTANCE is None:
        raise AssertionError('Settings were not set!')
    return _GLOBAL_SETTINGS_INSTANCE


_DEFAULT_OUT_DIR_PREFIX = 'results'
_SUFFIX = datetime.datetime.now().strftime('_%H%M%S')


def _get_out_dir(out_dir_prefix: Optional[str], out_dir: Optional[str]) -> str:
    if out_dir_prefix and out_dir:
        logging.fatal('Options --out-dir-prefix and --out-dir are mutually exclusive!')
        exit(1)
    elif out_dir_prefix:
        return os.path.join(os.getcwd(), out_dir_prefix + _SUFFIX)
    elif out_dir:
        return out_dir
    else:
        return os.path.join(os.getcwd(), _DEFAULT_OUT_DIR_PREFIX + _SUFFIX)
