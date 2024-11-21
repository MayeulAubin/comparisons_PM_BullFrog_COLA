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

__author__ = "Mayeul Aubin, Tristan Hoellinger"
__version__ = "0.1"
__date__ = "2024"
__license__ = "GPLv3"

"""
"""

from os.path import isfile
from pathlib import Path
from pickle import dump
import numpy as np

from pysbmy.power import PowerSpectrum
from pysbmy.fft import FourierGrid
from pysbmy import pySbmy

from params import *
from parser import ArgumentParser, intNone, bool_sh

parser = ArgumentParser(
    description="Run and compare PM, COLA and BullFrog simulations.",
)
parser.add_argument(
    "--workdir",
    type=str,
    help="Directory where the parameter files will be written.",
)
parser.add_argument(
    "--simdir_root",
    type=str,
    help="Directory where the output files of the simulations will be written.",
)
parser.add_argument(
    "--name",
    type=str,
    help="Suffix to the working directory for this run.",
    default="std",
)
parser.add_argument(
    "-pm",
    "--nsteps_pm",
    type=int,
    nargs="*",
    default=[],
    help="Number of time steps for PM",
)
parser.add_argument(
    "-cola",
    "--nsteps_cola",
    type=int,
    nargs="*",
    default=[],
    help="Number of time steps for COLA",
)
parser.add_argument(
    "-bf",
    "--nsteps_BullFrog",
    type=int,
    nargs="*",
    default=[],
    help="Number of time steps for BullFrog",
)
parser.add_argument("-L", type=int, help="Length of the simulation box in Mpc/h.", default=2000)
parser.add_argument(
    "--size",
    type=int,
    help="Number of elements of the simulation grid per dimension.",
    default=512,
)
parser.add_argument(
    "--Np",
    type=intNone,
    help="Number of DM particles per dimension.",
    default=1024,
)
parser.add_argument(
    "--Npm",
    type=intNone,
    help="Number of elements per dimension of the particle-mesh grid.",
    default=1024,
)
parser.add_argument("--RedshiftLPT", type=float, help="Redshift of the LPT pre-simulation")
parser.add_argument(
    "--RedshiftFCs", type=float, default=0.0, help="Redshift of the final conditions (output)"
)
parser.add_argument(
    "--verbosity",
    type=int,
    help="Verbosity level. 0: low verbosity, 1: only important outputs, 2: all outputs",
    default=1,
)
parser.add_argument("--force", type=bool_sh, help="Force the computations.", default=False)
parser.add_argument(
    "--seedphase", type=int, default=None, help="User-provided seed for the initial phase."
)

parser.add_argument("-ts","--timestepping",type=int, default=None, help="Time stepping distribution: 0 linear in a, 1 log in a, 2 exp in a, 3 bullfrog KDK linear in D")


args = parser.parse_args()
workdir = args.workdir
simdir_root = args.simdir_root
wd = workdir + args.name + "/"
simdir = simdir_root + args.name + "/"
nsteps_pm_list = args.nsteps_pm
nsteps_cola_list = args.nsteps_cola
nsteps_bullfrog_list = args.nsteps_BullFrog
L = args.L
N = args.size
Np = args.Np
Npm = args.Npm
RedshiftLPT = args.RedshiftLPT
RedshiftFCs = args.RedshiftFCs
verbose = args.verbosity
force = args.force
seedphase = args.seedphase if args.seedphase is not None else BASELINE_SEEDPHASE

time_stepping_distribution = args.timestepping
default_time_stepping_distribution = (time_stepping_distribution is None)

corner = -L / 2.0

logdir = simdir + "logs/"
nsim_pm = len(nsteps_pm_list)
nsim_cola = len(nsteps_cola_list)
nsim_bullfrog = len(nsteps_bullfrog_list)

Path(wd).mkdir(parents=True, exist_ok=True)
# Path(simdir).mkdir(parents=True, exist_ok=True)
Path(logdir).mkdir(parents=True, exist_ok=True)

params = {
    "simdir": simdir,
    "nsteps_pm": nsteps_pm_list,
    "nsteps_cola": nsteps_cola_list,
    "nsteps_bullfrog": nsteps_bullfrog_list,
    "L": L,
    "N": N,
    "Np": Np,
    "Npm": Npm,
    "RedshiftLPT": RedshiftLPT,
    "RedshiftFCs": RedshiftFCs,
    "verbose": verbose,
    "force": force,
}

# Save parameters for later use e.g. for plotting
with open(wd + "params.pkl", "wb") as f:
    dump(params, f)

# Save parameters for logging
with open(wd + "params.txt", "w") as f:
    f.write("Parameters for this run:\n")
    for key, value in params.items():
        f.write(f"{key}: {value}\n")

if __name__ == "__main__":
    from tools import generate_sim_params, generate_white_noise_Field

    print("################################")
    print("## Setting up main parameters ##")
    print("################################")

    # Path to the initial conditions (generated later)
    ICs_path = simdir + "initial_density.h5"

    # Path to the input matter power spectrum (generated later)
    input_power_file = simdir + "input_power.h5"

    # path to the initial white noise (generated later)
    input_white_noise_file = simdir + "input_white_noise.h5"
    input_seed_phase_file = simdir + "seed"

    # Initial and final scale factors
    ai = z2a(RedshiftLPT)
    af = z2a(RedshiftFCs)

    # Common parameters for the simulations
    common_params = {
        "Np": Np,
        "N": N,
        "L": L,
        "corner0": corner,
        "corner1": corner,
        "corner2": corner,
        "h": cosmo["h"],
        "Omega_m": cosmo["Omega_m"],
        "Omega_b": cosmo["Omega_b"],
        "n_s": cosmo["n_s"],
        "sigma8": cosmo["sigma8"],
    }

    # Parameters for LPT simulations
    lpt_params = common_params.copy()
    lpt_params["method"] = "lpt"
    lpt_params["InputPowerSpectrum"] = input_power_file
    lpt_params["ICsMode"] = 1
    # 0 : the codes generates white noise, then initial conditions
    # 1 : external white noise specified, the code multiplies by the power spectrum
    # 2 : external initial conditions specified
    lpt_params["InputWhiteNoise"] = input_white_noise_file

    # Parameters specific to PM
    pm_params = common_params.copy()
    pm_params["method"] = "pm"
    pm_params["TimeStepDistribution"] = 0 if default_time_stepping_distribution else time_stepping_distribution
    pm_params["ai"] = ai
    pm_params["af"] = af
    pm_params["RedshiftLPT"] = RedshiftLPT
    pm_params["RedshiftFCs"] = RedshiftFCs
    pm_params["Npm"] = Npm

    # Parameters specific to COLA
    cola_params = common_params.copy()
    cola_params["method"] = "cola"
    cola_params["TimeStepDistribution"] = 0 if default_time_stepping_distribution else time_stepping_distribution
    cola_params["ai"] = ai
    cola_params["af"] = af
    cola_params["RedshiftLPT"] = RedshiftLPT
    cola_params["RedshiftFCs"] = RedshiftFCs
    cola_params["Npm"] = Npm

    # Parameters specific to BullFrog
    bullfrog_params = common_params.copy()
    bullfrog_params["method"] = "bullfrog"
    bullfrog_params["TimeStepDistribution"] = 3 if default_time_stepping_distribution else time_stepping_distribution
    bullfrog_params["ai"] = ai
    bullfrog_params["af"] = af
    bullfrog_params["RedshiftLPT"] = RedshiftLPT
    bullfrog_params["RedshiftFCs"] = RedshiftFCs
    bullfrog_params["Npm"] = Npm

    print("#############################")
    print("## Writing parameter files ##")
    print("#############################")

    print("> Generating simulations cards...")

    generate_sim_params(lpt_params, ICs_path, wd, simdir, None, force)

    for i in range(nsim_pm):
        print(f"PM nsteps = {nsteps_pm_list[i]}:")
        pm_params["nsteps"] = nsteps_pm_list[i]
        file_ext = f"nsteps{nsteps_pm_list[i]}"  # "pm" is already in the filename
        generate_sim_params(pm_params, ICs_path, wd, simdir, file_ext, force)

    for i in range(nsim_cola):
        print(f"COLA nsteps = {nsteps_cola_list[i]}:")
        cola_params["nsteps"] = nsteps_cola_list[i]
        file_ext = f"nsteps{nsteps_cola_list[i]}"  # "cola" is already in the filename
        generate_sim_params(cola_params, ICs_path, wd, simdir, file_ext, force)

    for i in range(nsim_bullfrog):
        print(f"BULLFROG nsteps = {nsteps_bullfrog_list[i]}:")
        bullfrog_params["nsteps"] = nsteps_bullfrog_list[i]
        file_ext = f"nsteps{nsteps_bullfrog_list[i]}"  # "bullfrog" is already in the filename
        generate_sim_params(bullfrog_params, ICs_path, wd, simdir, file_ext, force)

    print("> Generating the input power spectrum...")
    if not isfile(input_power_file):
        Pk = PowerSpectrum(L, L, L, N, N, N, cosmo_small_to_full_dict(cosmo))
        Pk.write(input_power_file)

    print("> Generating the k grid...")
    Pinit = 100
    trim_threshold = 100  # Min number of modes required per bin for the summaries
    logkmin = np.log10(4 * np.pi / (np.sqrt(3) * L))
    kmax = np.pi * N / L
    Pbins_left_bnds = np.logspace(logkmin, np.log10(kmax), Pinit + 1, dtype=np.float32)
    Pbins_left_bnds = Pbins_left_bnds[:-1]
    input_ss_file = simdir + "input_ss_k_grid.h5"
    if not isfile(input_ss_file):
        Gk = FourierGrid(
            L,
            L,
            L,
            N,
            N,
            N,
            k_modes=Pbins_left_bnds,
            kmax=kmax,
            trim_bins=True,
            trim_threshold=trim_threshold,
        )
        Gk.write(input_ss_file)

    print("\n###########################")
    print("## Generate white noises ##")
    print("###########################")

    if not isfile(input_white_noise_file):
        rng = np.random.default_rng(seedphase)

        generate_white_noise_Field(
            seedphase=seedphase,
            fname_whitenoise=input_white_noise_file,
            seedname_whitenoise=input_seed_phase_file,
            N=N,
            L=L,
            corner=corner,
            force_phase=False,
        )

    print("\n#############################")
    print("## Running the simulations ##")
    print("#############################")

    if verbose < 2:
        from io import BytesIO
        from low_level import stdout_redirector, stderr_redirector

    print("> Starting LPT simulation...")
    if (
        not isfile(ICs_path)
        or not isfile(simdir + "lpt_density.h5")
        or not isfile(simdir + "lpt_particles.gadget3")
    ):
        fname_simparfile = f"{wd}example_lpt.sbmy"
        fname_simlogs = f"{logdir}lpt.txt"
        if verbose < 2:
            f = BytesIO()
            g = BytesIO()
            with stdout_redirector(f):
                with stderr_redirector(g):
                    pySbmy(fname_simparfile, fname_simlogs)
                g.close()
            f.close()
        else:
            pySbmy(fname_simparfile, fname_simlogs)

    print("> Starting PM simulations...")
    for i in range(nsim_pm):
        file_ext = f"nsteps{nsteps_pm_list[i]}"  # "pm" is already in the filename
        if not isfile(simdir + f"{file_ext}_final_density_pm.h5") or force:
            print(f">> Starting PM simulation {i}: {file_ext}")
            fname_simparfile = f"{wd}{file_ext}_example_pm.sbmy"
            fname_simlogs = f"{logdir}{file_ext}_pm.txt"
            if verbose < 2:
                f = BytesIO()
                g = BytesIO()
                with stdout_redirector(f):
                    with stderr_redirector(g):
                        pySbmy(fname_simparfile, fname_simlogs)
                    g.close()
                f.close()
            else:
                pySbmy(fname_simparfile, fname_simlogs)
        else:
            print(f">> PM simulation {i}: {file_ext} already completed")

    print("> Starting COLA simulations...")
    for i in range(nsim_cola):
        file_ext = f"nsteps{nsteps_cola_list[i]}"  # "cola" is already in the filename
        if not isfile(simdir + f"{file_ext}_final_density_cola.h5") or force:
            print(f">> Starting COLA simulation {i}: {file_ext}")
            fname_simparfile = f"{wd}{file_ext}_example_cola.sbmy"
            fname_simlogs = f"{logdir}{file_ext}_cola.txt"
            if verbose < 2:
                f = BytesIO()
                g = BytesIO()
                with stdout_redirector(f):
                    with stderr_redirector(g):
                        pySbmy(fname_simparfile, fname_simlogs)
                    g.close()
                f.close()
            else:
                pySbmy(fname_simparfile, fname_simlogs)
        else:
            print(f">> COLA simulation {i}: {file_ext} already completed")

    print("> Starting BullFrog simulations...")
    for i in range(nsim_bullfrog):
        file_ext = f"nsteps{nsteps_bullfrog_list[i]}"  # "bullfrog" is already in the filename
        if not isfile(simdir + f"{file_ext}_final_density_bullfrog.h5") or force:
            print(f">> Starting BullFrog simulation {i}: {file_ext}")
            fname_simparfile = f"{wd}{file_ext}_example_bullfrog.sbmy"
            fname_simlogs = f"{logdir}{file_ext}_bullfrog.txt"
            if verbose < 2:
                f = BytesIO()
                g = BytesIO()
                with stdout_redirector(f):
                    with stderr_redirector(g):
                        pySbmy(fname_simparfile, fname_simlogs)
                    g.close()
                f.close()
            else:
                pySbmy(fname_simparfile, fname_simlogs)
        else:
            print(f">> BullFrog simulation {i}: {file_ext} already completed")

    print("> Done!\n")
