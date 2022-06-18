'''
Adapter to RTMs. RTMs currently implemented

    - ProSAIL (4SAIL with either Prospect-5 or Prospect-D as leaf model)
    - SPART (BSM, 4SAIL, SMAC and Prospect-5 or Prospect-PRO)
'''

import numpy as np
import prosail

from pathlib import Path
from spectral import BandResampler
from typing import Optional

from rtm_inv.lookup_table import LookupTable
from rtm_inv._sensors import Sensors

class SPARTParameters:
    """
    class defining which leaf, canopy, soil and atmosphere parameters
    are required to run SPART simulations.

    This class helps mapping the entries of the CSV with the SPART model
    parameterization to the actual SPART function call so that user do not
    have to take care about the order of parameters.
    """

    __name__ = 'SPART'

    # define the entries required to run the SPART submodels
    SMB = ['B', 'lat', 'lon', 'SMp', 'SMC', 'film']  # soil model
    prospect_5d = ['Cab', 'Cca', 'Cw', 'Cdm', 'Cs', 'Cant', 'N', 'PROT', 'CBC']  # leaf model
    sailh = ['LAI', 'LIDFa', 'LIDFb', 'q']  # canopy model
    SMAC = ['aot550', 'uo3', 'uh2o', 'Pa']  # atmosphere model
    angles = ['sol_angle', 'obs_angle', 'rel_angle']  # sun and observer geometry


class ProSAILParameters:
    """
    class defining which leaf and canopy and soil parameters are required to run
    ProSAIL simulations
    """

    __name__ = 'prosail'

    prospect5 = ['n', 'cab', 'car', 'cbrown', 'cw', 'cm']
    prospectD = ['n', 'cab', 'car', 'cbrown', 'cw', 'cm', 'ant']
    fourSAIL= [
        'lai', 'lidfa', 'lidfb', 'psoil', 'rsoil', \
        'hspot', 'tts', 'tto', 'phi', 'typelidf', \
        'rsoil0', 'soil_spectrum1', 'soil_spectrum2', 'alpha'
    ]

class RTM:
    """
    Class for simulating synthetic vegetation spectra
    """
    def __init__(
            self,
            lut: LookupTable,
            rtm: str,
            n_step: Optional[int] = 500
        ):
        """
        Class constructor

        :param lut:
            lookup-table with vegetation traits and parameters for
            which to simulate spectra
        :param rtm:
            name of the RTM to run ("prosail", "SPART")
        :param n_step:
            step at which to write output to logs when created spectra
        """
        if lut.samples.empty:
            raise ValueError('LUT must not be empty')
        if rtm not in ['prosail', 'SPART']:
            raise ValueError('Unknown RTM name')
        if n_step <= 0:
            raise ValueError('Steps must be > 0')

        self._lut = lut
        self._rtm = rtm
        self._nstep = n_step

    def _run_prosail(self, sensor: str, **kwargs) -> None:
        """
        Runs the ProSAIL RTM
        """
        # check if Prospect version
        if set(ProSAILParameters.prospect5).issubset(set(self._lut.samples.columns)):
            prospect_version = '5'
        elif set(self._lut.samples.columns).issubset(ProSAILParameters.prospectD):
            prospect_version = 'D'
        else:
            raise ValueError('Cannot determine Prospect Version')

        # get sensor
        traits = self._lut.samples.columns
        try:
            sensor = eval(f'Sensors.{sensor}()')
        except Exception as e:
            raise Exception(f'No such sensor: {sensor}: {e}')

        # get band names
        sensor_bands = sensor.band_names
        self._lut.samples[sensor_bands] = np.nan
        # get central wavelengths and band width per band
        centers_sensor, fwhm_sensor = sensor.central_wvls, sensor.band_widths
        # convert band withs to FWHM (full-width-half-maximum)
        fwhm_sensor = [x*0.5 for x in fwhm_sensor]
        # define centers and bandwidth of ProSAIL output
        centers_prosail = np.arange(400,2501,1)
        fwhm_prosail = np.ones(centers_prosail.size)
        # initialize spectral sampler object to perform the spectral
        # convolution from 1nm ProSAIL output to Sentinel-2 spectral
        # resolution using a Gaussian spectral response function
        resampler = BandResampler(
            centers1=centers_prosail,
            centers2=centers_sensor,
            fwhm1=fwhm_prosail,
            fwhm2=fwhm_sensor
        )

        # iterate through LUT and run ProSAIL
        spectrum = None
        traits = self._lut.samples[traits].copy()
        for idx, record in traits.iterrows():
            record_inp = record.to_dict()
            record_inp.update({
                'prospect_version': prospect_version
            })
            record.update(kwargs)
            # run ProSAIL
            try:
                spectrum = prosail.run_prosail(**record_inp)
            except Exception as e:
                print(e)
            if (idx+1)%self._nstep == 0:
                print(f'Simulated spectrum {idx+1}/{self._lut.samples.shape[0]}')

            # resample to spectral resolution of sensor
            sensor_spectrum = resampler(spectrum)
            self._lut.samples.at[idx,sensor_bands] = sensor_spectrum

    def simulate_spectra(self, sensor: str, **kwargs):
        """
        Simulation of spectra for all entries in the lookup-table

        :paran sensor:
            name of the sensor for which to generate spectra
        """
        # call different RTMs
        if self._rtm == 'prosail':
            self._run_prosail(sensor=sensor, **kwargs)

        return self._lut.samples

if __name__ == '__main__':

    import matplotlib.pyplot as plt

    fpath_rtm_params = Path('../../parameters/prosail_s2.csv')

    lut = LookupTable(params_csv=fpath_rtm_params)
    num_samples = 50000
    method = 'LHS'
    lut.generate_samples(num_samples, method)

    rtm = RTM(lut=lut, rtm='prosail')
    s2_rtm_simulations = rtm.simulate_spectra(sensor='Sentinel2A')

    # save to CSV
    s2_rtm_simulations.to_csv('../../parameters/s2_prosail_demo.csv', index=False)

    s2a = Sensors.Sentinel2A
    plt.plot(
        s2a.central_wvls,
        s2_rtm_simulations[s2a.band_names].values.T
    )
    plt.ylabel('Surface Reflectance Factor [-]')
    plt.xlabel('Wavelength [nm]')
    plt.show()
