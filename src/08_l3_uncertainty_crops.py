'''
Extracts the uncertainty in the different levels for the selected field parcel polygons
and relates to the different crop types available. This way it is also possible to plot
uncertainty over time per crop type.
'''

import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime
from typing import Optional
from typing import List
from pathlib import Path

plt.style.use('ggplot')


def extract_uncertainty_crops(
        analysis_dir: Path,
        parameter_name: str,
        out_dir: Path
    ):
    """
    Extracts the uncertainty information for the field parcel polygons.

    :param analysis_dir:
        directory where the uncertainty results in the L1C, L2A and L3 products
        are stored
    :param vi_name:
        name of the vegetation index/ parameter to analyze
    :param out_dir:
        directory where to save the results of the analysis to
    """

    # find results for the scenes available 



if __name__ == '__main__':
    
    analysis_dir = Path(
        '../S2A_MSIL2A_Analysis/150_scenarios'
    )

    output_dir = analysis_dir.joinpath('uncertainty_crops')
    if not output_dir.exists(): output_dir.mkdir()

    
        