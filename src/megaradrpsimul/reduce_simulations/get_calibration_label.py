#
# Copyright 2025-2026 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Determine the calibration tree from the INSCONF keyword
    of the simulated raw MEGARA images."""

from astropy.io import fits
from pathlib import Path

import logging
logger = logging.getLogger(__name__)

def get_calibration_label(data_dir, 
                         calibrations_dir):
    """
    Determine the calibration tree from the INSCONF keyword
    of the simulated raw MEGARA images.

    Parameters
    ----------
    data_dir : Path instance
        Path to the directory containing the simulated raw MEGARA images.
    calibrations_dir : Path instance
        Path to the directory containing the calibration trees.
    """

    raw_files = sorted(data_dir.glob("0*.fits"))

    if not raw_files:
        raise FileNotFoundError(
            f"No simulated raw MEGARA images found in {data_dir}"
        )

    insconf_values = set()
    for raw_file in raw_files:
        header = fits.getheader(raw_file, 0)
        insconf = header.get("INSCONF")

        if not insconf:
            raise KeyError(
                f"INSCONF keyword not found in {raw_file}"
            )

        insconf_values.add(str(insconf).strip())

    if len(insconf_values) != 1:
        raise ValueError(
            "Different INSCONF values were found in the simulated images: "
            f"{sorted(insconf_values)}"
        )

    insconf = insconf_values.pop()
    calibration_label_dir = calibrations_dir / insconf

    if not calibration_label_dir.is_dir():
        raise FileNotFoundError(
            f"Calibration tree for INSCONF={insconf} was not found: "
            f"{calibration_label_dir}"
        )

    logger.info("MEGARA instrument configuration: %s", insconf)
    logger.debug("Calibration root: %s", calibration_label_dir)

    return calibration_label_dir