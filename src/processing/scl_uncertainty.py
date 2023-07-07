'''
Uncertainty in the scene classification layer (SCL).
'''

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import mode
from uncertainties import ufloat
from pathlib import Path
from sklearn.metrics import confusion_matrix
from pretty_print_confusion import pp_matrix

from eodal.core.sensors import Sentinel2
from eodal.utils.constants.sentinel2 import SCL_Classes

from utils.logger import get_logger

logger = get_logger('SCL_Uncertainty')
scl_classes = SCL_Classes().values()

def scl_uncertainty(
        scenario_dir: Path,
        output_dir: Path,
        aoi: Path
    ) -> None:
    """
    Classification uncertainty of the scene classification layer (SCL)
    of the L2A Sentinel-2 product

    :param scenario_dir:
        directory with L1C/L2A scenario runs
    :param output_dir:
        directory where to save the results to
    :param aoi:
        vector file defining the area of interest to read
    """
    # find scenes for which scenarios are available
    scenes = glob.glob(scenario_dir.joinpath('S2*_MSIL1C*').as_posix())

    # loop over scenes and their scenarios
    area_stats_scenes_list = []
    class_assignment_confidence_list = []
    for idx, scene in enumerate(scenes):

        logger.info(f'Working on scene {scene} ({idx+1}/{len(scenes)})')
        # find L2A scenes
        scenarios = glob.glob(Path(scene).joinpath('*/S2*_MSIL2A*.SAFE').as_posix())

        # loop over scenarios of the current scene
        area_stats_list = []
        scl_array_list = []
        n_scenarios = len(scenarios)
        for jdx, scenario in enumerate(scenarios):

            logger.info(f'Processing Scene {jdx+1}/{len(scenarios)} ({scenario})')
            # define input and outputs for the LAI model
            scenario = Path(scenario)
            s2_ds = Sentinel2().from_safe(
                in_dir=scenario,
                band_selection=['SCL'],
                vector_features=aoi
            )
            # get scl layer
            scl = s2_ds.get_values()
            # get relative number of pixels per class
            scl_vals, scl_counts = np.unique(scl, return_counts=True)
            scl_stats = dict.fromkeys(scl_classes.values(), 0)
            for scl_class in scl_classes:
                if (scl_vals == scl_class).any():
                    class_idx = np.where(scl_vals == scl_class)
                    # calculate percentage of pixels belonging to that SCL class
                    class_count = scl_counts[class_idx][0] / scl.size * 100
                    scl_stats[scl_classes[scl_class]] = class_count
            area_stats_list.append(scl_stats)
            scl_array_list.append(scl)

        # calculate the uncertainties. For the number of pixels per class check how this
        # number differs among the scenario runs:
        area_stats_raw = pd.DataFrame(area_stats_list)
        mean_pixel_num_per_class = area_stats_raw.mean()
        std_pixel_num_per_class = area_stats_raw.std()
        # combine mean and std into a uncertainties object
        area_stats_unc = {}
        area_stats_unc['date'] = pd.to_datetime(Path(scene).name.split('_')[2][0:8])
        for scl_class in scl_classes.values():
            mu = mean_pixel_num_per_class.at[scl_class]
            sigma = std_pixel_num_per_class.at[scl_class]
            # convert into relative uncertainties
            sigma_rel = sigma / mu * 100.
            area_stats_unc[scl_class] = ufloat(nominal_value=mu, std_dev=sigma_rel)
        area_stats_scenes_list.append(area_stats_unc)

        # For the class assignment confidence use scipy.stats mode and check how
        # the class assignment confidence varies within a class (always looking at
        # the majority vote of a pixel
        data_arr = np.vstack(scl_array_list)
        scl_analysis = mode(data_arr, axis=0)

        # majority vote
        majority = scl_analysis.mode[:,:]
        # confidence of the majority vote
        confidence = scl_analysis.count[:,:] / n_scenarios * 100.

        # create confusion matrix
        # the majority vote serves as the "truth" against which the scenario outcomes are evaluated
        y_test = np.empty(shape=(len(scl_array_list), majority.shape[0], majority.shape[1]))
        y_pred = np.empty_like(y_test)
        for ii in range(len(scl_array_list)):
            y_test[ii,:,:] = majority
            y_pred[ii,:,:] = scl_array_list[ii][0,:,:]
        y_test = y_test.flatten()
        y_pred = y_pred.flatten()
        conf_matrix = confusion_matrix(y_true=y_test, y_pred=y_pred)
        # label classes correctly based on their occurence in the scene
        unique_true_scl_classes = np.unique(majority)
        unique_pred_scl_classes = np.unique(y_pred).astype(int)
        # check how many SCL classes are in the confusion matrix
        if len(unique_true_scl_classes) < len(unique_pred_scl_classes):
            unique_scl_classes = unique_pred_scl_classes
        elif len(unique_true_scl_classes) > len(unique_pred_scl_classes):
            unique_scl_classes = unique_true_scl_classes
        else:
            unique_scl_classes = unique_true_scl_classes
        full_conf_matrix = np.zeros(shape=(12,12), dtype=int)
        for ii, unique_scl_class in enumerate(unique_scl_classes):
            for jj, _unique_scl_class in enumerate(unique_scl_classes):
                full_conf_matrix[unique_scl_class,_unique_scl_class] = conf_matrix[ii,jj]
        # convert to pandas data-frame and assign scl class labels
        scene_date = datetime.strptime((Path(scene).name.split('_')[2][0:8]), '%Y%m%d').date()
        conf_df = pd.DataFrame(
            full_conf_matrix,
            columns=list(scl_classes.values()),
            index=list(scl_classes.values())
        )
        fname_conf = output_dir.joinpath(f'{scene_date}_SCL_conf.png')
        f = pp_matrix(conf_df, cmap='Accent_r', figsize=(12,12))
        f.savefig(fname_conf, bbox_inches='tight')
        plt.close(f)

        fname_csv = output_dir.joinpath(f'{scene_date}_SCL_conf.csv')
        conf_df.to_csv(fname_csv)

        # go through the majority votes and check the variability of the majority
        # vote confidence
        majority_confidence_variability = {}
        majority_confidence_variability['date'] = pd.to_datetime(Path(scene).name.split('_')[2][0:8])
        majority_classes = np.unique(majority)
        for scl_class in scl_classes:
            if (scl_class == majority_classes).any():
                class_idx = np.where(majority_classes == scl_class)
                # check variability of confidence values of that class
                confidence_variability = np.std(confidence[majority == scl_class])
                confidence_mean = np.mean(confidence[majority == scl_class])
                confidence_unc = ufloat(
                    nominal_value=confidence_mean,
                    std_dev=confidence_variability
                )
                majority_confidence_variability[scl_classes[scl_class]] = confidence_unc
        class_assignment_confidence_list.append(majority_confidence_variability)

    # combine results from all scenes
    area_stats_complete = pd.DataFrame(area_stats_scenes_list)
    class_assignment_confidence_complete = pd.DataFrame(class_assignment_confidence_list)

    # save results to CSV
    fname_area_stats = output_dir.joinpath('SCL_relative-number-of-pixels-per-class_abs-uncertainty.csv')
    area_stats_complete.to_csv(fname_area_stats, index=False)

    fname_classignment_confidence = output_dir.joinpath('SCL_class-assignment-confidence-abs-variability.csv')
    class_assignment_confidence_complete.to_csv(fname_classignment_confidence, index=False)

if __name__ == '__main__':

    aoi = Path('../../shp/AOI_Esch_EPSG32632.shp')

    scenario_dir = Path('../../S2_MSIL1C_RUT-Scenarios')
    output_dir = scenario_dir.joinpath('SCL_Uncertainty')
    output_dir.mkdir(exist_ok=True)

    scl_uncertainty(
        scenario_dir=scenario_dir,
        output_dir=output_dir,
        aoi=aoi
    )

    # combine results from batches into two single DataFrames and save them to csv and latex
    area_stats_list = []
    class_confidence_list = []
    
    fname_area_stats = Path(
        '../../S2_MSIL1C_RUT-Scenarios/SCL_Uncertainty/SCL_relative-number-of-pixels-per-class_abs-uncertainty.csv'
    )
    area_stats = pd.read_csv(fname_area_stats)
    fname_class_confidence = Path(
        '../../S2_MSIL1C_RUT-Scenarios/SCL_Uncertainty/SCL_class-assignment-confidence-abs-variability.csv'
    )
    class_confidence = pd.read_csv(fname_class_confidence)

    fname_area_stats = Path('../../S2_MSIL2A_Analysis').joinpath(
        'SCL_relative-number-of-pixels-per-class_abs-uncertainty.csv'
    )
    area_stats = area_stats.sort_values(by='date')
    area_stats.to_csv(fname_area_stats, index=False)

    fname_class_confidence = Path('../../S2_MSIL2A_Analysis').joinpath(
        'SCL_class-assignment-confidence-abs-variability.csv'
    )
    class_confidence = class_confidence.sort_values(by='date')
    class_confidence.to_csv(fname_class_confidence, index=False)

    # output latex table (optional)
    # df = pd.read_csv(fname_area_stats)
    # sel_cols = ['date', 'cloud_shadows', 'vegetation', 'non_vegetated', 'water', 'unclassified',
    #    'cloud_medium_probability', 'cloud_high_probability', 'thin_cirrus']
    # latex_expr = df[sel_cols].to_latex(index=False)
