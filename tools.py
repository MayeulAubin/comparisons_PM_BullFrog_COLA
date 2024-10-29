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
        Directory where the parameter file will be written.
    outdir : str
        Directory where the output files of the simulations will be
        written.
    file_ext : str, optional
        Prefix to add to the output file name.

    Returns
    -------
    sbmy_path : str
        Path to the created parameter file.

    """
    from os.path import isfile
    from pysbmy import param_file
    from pysbmy.timestepping import StandardTimeStepping, BullFrogTimeStepping
    from params import cosmo, cosmo_small_to_full_dict

    method = params_dict["method"]
    path = workdir + file_ext + "_" if file_ext else workdir
    simpath = outdir + file_ext + "_" if file_ext else outdir
    sbmy_path = path + "example_" + method + ".sbmy"

    # Parameters shared by all methods for this run
    Particles = params_dict["Np0"]
    Mesh = params_dict["N0"]
    BoxSize = params_dict["L0"]
    corner0 = params_dict["corner0"]
    corner1 = params_dict["corner1"]
    corner2 = params_dict["corner2"]
    h = params_dict["h"]
    Omega_m = params_dict["Omega_m"]
    Omega_b = params_dict["Omega_b"]
    n_s = params_dict["n_s"]
    sigma8 = params_dict["sigma8"]

    # Generate the time-stepping distribution
    if method in ["cola", "pm"]:
        ts_filename = path + "ts_" + method + ".h5"
        print(">> Generating time-stepping distribution...")
        if not isfile(ts_filename) or force:
            TimeStepDistribution = params_dict["TimeStepDistribution"]
            ai = params_dict["ai"]
            af = params_dict["af"]
            nsteps = params_dict["nsteps"]
            forces = np.full(nsteps, True)
            snapshots = np.full((nsteps), False)
            TS = StandardTimeStepping(ai, af, snapshots, TimeStepDistribution, forces=forces)
            TS.write(ts_filename)
        else:
            print(">> Using existing time-stepping distribution.")
    elif method == "bullfrog":
        ts_filename = path + "ts_bullfrog.h5"
        print(">> Generating time-stepping distribution...")
        if not isfile(ts_filename) or force:
            # TimeStepDistribution = params_dict["TimeStepDistribution"]
            ai = params_dict["ai"]
            af = params_dict["af"]
            nsteps = params_dict["nsteps"]
            forces = np.full(nsteps, True)
            snapshots = np.full((nsteps), False)
            TS = BullFrogTimeStepping(ai, af, cosmo_small_to_full_dict(cosmo), snapshots, forces=forces)
            TS.write(ts_filename)
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
            ParticleMesh=params_dict["Npm0"],
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
            ParticleMesh=params_dict["Npm0"],
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
            ParticleMesh=params_dict["Npm0"],
            TimeStepDistribution=ts_filename,
            RedshiftFCs=params_dict["RedshiftFCs"],
            WriteDensities=1,
            OutputDensitiesBase=simpath + "density_bullfrog_",
            WriteSnapshots=1,
            OutputSnapshotsBase=simpath + "particles_bullfrog_",
            WriteFinalSnapshot=1,
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


### This is a snippet to redirect C stdout to a file
### https://dzone.com/articles/redirecting-all-kinds-stdout

from contextlib import contextmanager
import ctypes
import io
import os, sys
import tempfile

libc = ctypes.CDLL(None)
c_stdout = ctypes.c_void_p.in_dll(libc, 'stdout')

@contextmanager
def stdout_redirector(stream):
    # The original fd stdout points to. Usually 1 on POSIX systems.
    original_stdout_fd = sys.stdout.fileno()

    def _redirect_stdout(to_fd):
        """Redirect stdout to the given file descriptor."""
        # Flush the C-level buffer stdout
        libc.fflush(c_stdout)
        # Flush and close sys.stdout - also closes the file descriptor (fd)
        sys.stdout.close()
        # Make original_stdout_fd point to the same file as to_fd
        os.dup2(to_fd, original_stdout_fd)
        # Create a new sys.stdout that points to the redirected fd
        sys.stdout = io.TextIOWrapper(os.fdopen(original_stdout_fd, 'wb'))

    # Save a copy of the original stdout fd in saved_stdout_fd
    saved_stdout_fd = os.dup(original_stdout_fd)
    try:
        # Create a temporary file and redirect stdout to it
        tfile = tempfile.TemporaryFile(mode='w+b')
        _redirect_stdout(tfile.fileno())
        # Yield to caller, then redirect stdout back to the saved fd
        yield
        _redirect_stdout(saved_stdout_fd)
        # Copy contents of temporary file to the given stream
        tfile.flush()
        tfile.seek(0, io.SEEK_SET)
        stream.write(tfile.read())
    finally:
        tfile.close()
        os.close(saved_stdout_fd)
