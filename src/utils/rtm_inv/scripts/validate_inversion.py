'''
Created on Feb 20, 2022

@author: graflu
'''

import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from rtm_inv.inversion import inv_img, retrieve_traits
from scipy.stats import linregress
from sklearn.metrics import mean_squared_error

plt.style.use('ggplot')

if __name__ == '__main__':

    s2_val_data = Path('/mnt/ides/Lukas/02_Research/PhenomEn/01_Data/04_inSitu_Data')
    search_expr = '*_WW_greenLAI_S2spectra.csv'
    lut_path = Path('/mnt/ides/Lukas/software/scripts_paper_uncertainty/S2_ProSAIL_LUTs/S2A_MSIL2A_20190530T103031_N9999_R108_T32TMT_20211122T115309-prosail_lut.pkl')
    traits = ['lai']

    # read S2 spectra and in-situ data
    df_list = []
    for s2_val_dataset in s2_val_data.glob(search_expr):
        df_list.append(pd.read_csv(s2_val_dataset))

    s2_val = pd.concat(df_list)
    s2_val['date'] = pd.to_datetime(s2_val['date'])
    s2_val['year'] = s2_val['date'].dt.year
    s2_band_names = s2_val.columns[s2_val.columns.str.startswith('B')]
    s2_spectra = s2_val[s2_band_names].values.T * 0.0001
    s2_spectra_shape = (len(s2_band_names), s2_val.shape[0], 1)
    s2_spectra = s2_spectra.reshape(s2_spectra_shape)
    mask = np.zeros(shape=s2_spectra_shape[1::], dtype='uint8').astype('bool')

    # load LUT and run inversion on S2 spectra available
    s2_lut = pd.read_pickle(lut_path)
    s2_lut_spectra = s2_lut[['B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B8A', 'B11', 'B12']].values

    n_solutions = [1, 5, 10, 20, 100, 200, 1000, 2000, 2500, 5000]
    cost_functions = ['rmse', 'mae', 'contrast_function']

    inv_configs = itertools.product(n_solutions, cost_functions)

    res = []
    for inv_config in inv_configs:
        n_solution = inv_config[0]
        cost_function = inv_config[1]
        lut_idxs = inv_img(
            lut=s2_lut_spectra,
            img=s2_spectra,
            mask=mask,
            cost_function=cost_function,
            n_solutions=n_solution
        )
        trait_vals = retrieve_traits(
            lut=s2_lut,
            lut_idxs=lut_idxs,
            traits=traits
        )
        # compute error statistics
        rmse = mean_squared_error(
            s2_val['greenLAI'].values, trait_vals[0,:,0], squared=False
        )
        nrmse = rmse / s2_val['greenLAI'].mean() * 100
        slope, intercept, r_value, p_value, std_err = linregress(
            s2_val['greenLAI'].values,
            trait_vals[0,:,0]
        )
        textstr = f'RMSE = {np.round(rmse,2)}' + r'$m^2$/$m^2$'
        textstr += f'\nnRMSE = {np.round(nrmse,2)}%\n'
        textstr += r'$R^2$ = ' + f'{np.round(r_value**2,2)}'
    
        # compare against in-situ data
        f = plt.figure(figsize=(6,6))
        ax = f.add_subplot(111)
        ax.scatter(x=s2_val['greenLAI'], y=trait_vals[0,:,0])
        ax.set_xlabel(r'In-situ measured green LAI [$m^2$/$m^2$]', fontsize=16)
        ax.set_ylabel(r'ProSAIL estimated green LAI [$m^2$/$m^2$]', fontsize=16)
        ax.set_title(
            f'Winter Wheat (N={s2_val.shape[0]})\nCost Function: {cost_function}, '\
            f'#Solutions: {n_solution}'
        )
        ax.set_xlim(0,7)
        ax.set_ylim(0,7)
        ax.plot([0,7], [0,7], 'k-', label='1:1 Line')
        ax.legend()
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.05, 0.9, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
    
        f.savefig(
            s2_val_data.joinpath(f'ProSAIL-insitu_greenLAI20172018_{cost_function}-N{n_solution}.png'),
            bbox_inches='tight'
        )
        plt.close(f)

        res_dict = {
            'cost_function': cost_function,
            'n_solutions': n_solution,
            'rmse': rmse,
            'nrmse': nrmse,
            'r_squared': r_value**2,
            'slope': slope,
            'intercept': intercept,
            'p_value': p_value
        }
        res.append(res_dict)

    error_stats = pd.DataFrame(res)
    error_stats.to_csv(s2_val_data.joinpath('error_stats.csv'), index=False)
    

