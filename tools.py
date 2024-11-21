#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------------------
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# The text of the license is located in the root directory of the source package.
# -------------------------------------------------------------------------------------

__author__ = "Tristan Hoellinger"
__version__ = "0.1"
__date__ = "2024"
__license__ = "GPLv3"

"""
"""

import numpy as np


def generate_sim_params(params_dict, ICs, workdir, outdir, file_ext=None, force=False):
    """Write the parameter file.

    Parameters
    ----------
    params_dict : dict
        Dictionary containing the parameters for the simulation.
    ICs : str
        Path to the initial conditions.
    workdir : str
        Directory where to store the parameter file.
    outdir : str
        Directory where to store the simulation outputs.
    file_ext : str, optional
        Prefix for the output files.

    Returns
    -------
    sbmy_path : str
        Path to the parameter file generated.

    """
    from os.path import isfile
    from pysbmy import param_file
    from pysbmy.timestepping import StandardTimeStepping, BullFrogTimeStepping
    from params import cosmo, cosmo_small_to_full_dict
    from pysbmy.timestepping import StandardTimeStepping, BullFrogTimeStepping
    from params import cosmo, cosmo_small_to_full_dict

    method = params_dict["method"]
    path = workdir + file_ext + "_" if file_ext else workdir
    simpath = outdir + file_ext + "_" if file_ext else outdir
    sbmy_path = path + "example_" + method + ".sbmy"

    # Parameters shared by all methods for this run
    Particles = params_dict["Np"]
    Mesh = params_dict["N"]
    BoxSize = params_dict["L"]
    corner0 = params_dict["corner0"]
    corner1 = params_dict["corner1"]
    corner2 = params_dict["corner2"]
    h = params_dict["h"]
    Omega_m = params_dict["Omega_m"]
    Omega_b = params_dict["Omega_b"]
    n_s = params_dict["n_s"]
    sigma8 = params_dict["sigma8"]

    # Generate the time-stepping distribution
    if method!="lpt":
        ts_filename = path + "ts_" + method + ".h5"
        print(">> Generating time-stepping distribution...")
        if not isfile(ts_filename) or force:
            TimeStepDistribution = params_dict["TimeStepDistribution"]
            ai = params_dict["ai"]
            af = params_dict["af"]
            nsteps = params_dict["nsteps"]
            forces = np.full(nsteps, True)
            forces = np.full(nsteps, True)
            snapshots = np.full((nsteps), False)
            if TimeStepDistribution != 3:
                if method!="bullfrog":
                    TS = StandardTimeStepping(ai, af, snapshots, TimeStepDistribution)
                else:
                    TS = StandardTimeStepping(ai, af, np.full((nsteps), True), TimeStepDistribution)
                    TS.snapshots*=False
            else:
                snapshots = np.full((nsteps), True)
                TS = BullFrogTimeStepping(ai, af, cosmo_small_to_full_dict(cosmo), snapshots, forces=forces)
            TS.write(ts_filename)
            TS.plot(savepath=path + "ts_" + method + ".png")
        else:
            print(">> Using existing time-stepping distribution.")

    # Write the parameter file
    print(">> Generating parameter file...")
    if params_dict["method"] == "lpt":
        S = param_file(
            OutputRngStateLPT=simpath + "dummy.rng",
            Particles=Particles,
            Mesh=Mesh,
            BoxSize=BoxSize,
            corner0=corner0,
            corner1=corner1,
            corner2=corner2,
            ICsMode=params_dict["ICsMode"],
            InputWhiteNoise=params_dict["InputWhiteNoise"],  # None or str
            InputPowerSpectrum=params_dict["InputPowerSpectrum"],
            WriteInitialConditions=1,
            OutputInitialConditions=ICs,
            OutputLPTSnapshot=simpath + "lpt_particles.gadget3",
            OutputLPTDensity=simpath + "lpt_density.h5",
            h=h,
            Omega_m=Omega_m,
            Omega_b=Omega_b,
            n_s=n_s,
            sigma8=sigma8,
            Omega_q=1.0 - Omega_m,
            Omega_k=0.0,
            w0_fld=-1.0,
            wa_fld=0.0,
        )
    if params_dict["method"] == "pm":
        S = param_file(
            Particles=Particles,
            Mesh=Mesh,
            BoxSize=BoxSize,
            corner0=corner0,
            corner1=corner1,
            corner2=corner2,
            ICsMode=2,
            InputInitialConditions=ICs,
            RedshiftLPT=params_dict["RedshiftLPT"],
            WriteLPTSnapshot=0,
            WriteLPTDensity=0,
            ModulePMCOLA=1,
            EvolutionMode=1,
            ParticleMesh=params_dict["Npm"],
            TimeStepDistribution=ts_filename,
            RedshiftFCs=params_dict["RedshiftFCs"],
            WriteFinalDensity=1,
            OutputFinalDensity=simpath + "final_density_pm.h5",
            h=h,
            Omega_m=Omega_m,
            Omega_b=Omega_b,
            n_s=n_s,
            sigma8=sigma8,
            Omega_q=1.0 - Omega_m,
            Omega_k=0.0,
            w0_fld=-1.0,
            wa_fld=0.0,
        )
    elif params_dict["method"] == "cola":
        S = param_file(
            Particles=Particles,
            Mesh=Mesh,
            BoxSize=BoxSize,
            corner0=corner0,
            corner1=corner1,
            corner2=corner2,
            ICsMode=2,
            InputInitialConditions=ICs,
            RedshiftLPT=params_dict["RedshiftLPT"],
            WriteLPTSnapshot=0,
            WriteLPTDensity=0,
            ModulePMCOLA=1,
            EvolutionMode=2,
            ParticleMesh=params_dict["Npm"],
            TimeStepDistribution=ts_filename,
            RedshiftFCs=params_dict["RedshiftFCs"],
            WriteDensities=1,
            OutputDensitiesBase=simpath + "density_cola_",
            WriteSnapshots=1,
            OutputSnapshotsBase=simpath + "particles_cola_",
            WriteFinalSnapshot=0,
            WriteFinalDensity=1,
            OutputFinalDensity=simpath + "final_density_cola.h5",
            h=h,
            Omega_m=Omega_m,
            Omega_b=Omega_b,
            n_s=n_s,
            sigma8=sigma8,
            Omega_q=1.0 - Omega_m,
            Omega_k=0.0,
            w0_fld=-1.0,
            wa_fld=0.0,
        )
    elif params_dict["method"] == "bullfrog":
        S = param_file(
            Particles=Particles,
            Mesh=Mesh,
            BoxSize=BoxSize,
            corner0=corner0,
            corner1=corner1,
            corner2=corner2,
            ICsMode=2,
            InputInitialConditions=ICs,
            RedshiftLPT=params_dict["RedshiftLPT"],
            WriteLPTSnapshot=0,
            WriteLPTDensity=0,
            ModulePMCOLA=1,
            EvolutionMode=4,
            ParticleMesh=params_dict["Npm"],
            TimeStepDistribution=ts_filename,
            RedshiftFCs=params_dict["RedshiftFCs"],
            WriteDensities=1,
            OutputDensitiesBase=simpath + "density_bullfrog_",
            WriteSnapshots=1,
            OutputSnapshotsBase=simpath + "particles_bullfrog_",
            WriteFinalSnapshot=0,
            WriteFinalDensity=1,
            OutputFinalDensity=simpath + "final_density_bullfrog.h5",
            h=h,
            Omega_m=Omega_m,
            Omega_b=Omega_b,
            n_s=n_s,
            sigma8=sigma8,
            Omega_q=1.0 - Omega_m,
            Omega_k=0.0,
            w0_fld=-1.0,
            wa_fld=0.0,
        )

    if not isfile(sbmy_path) or force:
        S.write(sbmy_path)
        print(">> Parameter file written to {}.".format(sbmy_path))
    else:
        print(">> Parameter file already exists.")

    return sbmy_path


def read_field(*args):
    from io import BytesIO
    from low_level import stdout_redirector
    from pysbmy.field import read_field as _read_field

    with BytesIO() as f:
        with stdout_redirector(f):
            return _read_field(*args)


def generate_white_noise_Field(
    seedphase,
    fname_whitenoise,
    seedname_whitenoise,
    L,
    N,
    corner,
    force_phase=False,
):
    """
    Generates a white noise realization in physical space using SeedSequence and
    default_rng from numpy.random.

    Parameters
    ----------
    seedphase : int or list of int
        user-provided seed to generate the “initial” white noise realization
    fname_whitenoise : str
        name of the output white noise file
    seedname_whitenoise : str
        name of the output white noise seed file
    force_phase : bool, optional, default=False
        force recomputation of the white noise
    L : float
        the size in Mpc of the box
    N : int
        the number of grid points
    corner : float
        the position of the corner in Mpc
    """
    from os.path import exists

    if not exists(fname_whitenoise) or force_phase:
        from gc import collect
        from numpy import random, save
        from pysbmy.field import BaseField

        rng = random.default_rng(seedphase)
        save(seedname_whitenoise, rng.bit_generator.state)
        with open(seedname_whitenoise + ".txt", "w") as f:
            f.write(str(rng.bit_generator.state))
        data = rng.standard_normal(size=N**3)
        wn = BaseField(L, L, L, corner, corner, corner, 1, N, N, N, data)
        del data
        wn.write(fname_whitenoise)
        del wn
        collect()
