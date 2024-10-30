from os.path import isfile
from pathlib import Path
import numpy as np
import io

from pysbmy.power import PowerSpectrum
from pysbmy.fft import FourierGrid
from pysbmy import pySbmy

from tools import stdout_redirector as redirect_stdout

import argparse
    
parser = argparse.ArgumentParser(
                    prog='RunPMCOLABullFrog',
                    description='Runs simulations of PM, COLA and BullFrog to compare them',)

parser.add_argument('i', metavar='i', type=int, default=0, help='run id')
parser.add_argument('-pm','--nsteps_pm', type=int, nargs="*", default=[20, 50, 100, 200], help='number of time steps to simulate for PM')
parser.add_argument('-co','--nsteps_cola', type=int, nargs="*", default=[2, 5, 10, 20, 50], help='number of time steps to simulate for COLA')
parser.add_argument('-bf','--nsteps_BullFrog', type=int, nargs="*", default=[2, 5, 10, 20, 50], help='number of time steps to simulate for BullFrog')
parser.add_argument('-L', type=float, default=250.0, help='Length of the box (Mpc)')
parser.add_argument('-N', type=int, default=128, help='LPT mesh size per dimension')
parser.add_argument('-Np', type=int, default=128, help='Number of particles per dimension')
parser.add_argument('-Npm', type=int, default=128, help='Density/Potential mesh size per dimension')
parser.add_argument('--RedshiftLPT', type=float, default=24.0, help='Redshift of the LPT pre-simulation')
parser.add_argument('--RedshiftFCs', type=float, default=0.0, help='Redshift of the final conditions (output)')
parser.add_argument('-V','--verbose_smby',action='store_true',help='Display the Simbelmynë stdout in the console')
parser.add_argument('-S', '--seed', type=int, default=None, help='User provided seed for the white noise')


args = parser.parse_args()

print("\n================================")
print(f"\tRunPMCOLABullFrog\t Run ID: {args.i}")
print("\n================================\n")

print("\n> Getting parameters...\n")

WORKDIR = "runs_params/"
SIMDIR = (
    "runs_sims/"  # need not be on same disk
)

from tools import generate_sim_params, generate_white_noise_Field
from params import (
    cosmo,
    z2a,
    cosmo_small_to_full_dict,
)

run_id = f"run{args.i}"
force = False

# Parameters for the pm simulations
nsteps_pm_list = args.nsteps_pm # , 30]

# Parameters for the cola simulations
# nsteps_cola_list = [2, 5, 10]  # , 15]
nsteps_cola_list = args.nsteps_cola

# Parameters for the cola simulations
nsteps_bullfrog_list = args.nsteps_BullFrog  # , 15]


wd = WORKDIR + run_id + "/"
simdir = SIMDIR + run_id + "/"
logdir = simdir + "logs/"
Path(wd).mkdir(parents=True, exist_ok=True)
Path(logdir).mkdir(parents=True, exist_ok=True)

nsim_pm = len(nsteps_pm_list)
nsim_cola = len(nsteps_cola_list)
nsim_bullfrog = len(nsteps_bullfrog_list)

ICs_path = simdir + "initial_density.h5"
simpath = simdir

# Path to the input matter power spectrum (generated later)
input_power_file = simdir + "input_power.h5"

# path to the initial white noise (generated later)
input_white_noise_file = simdir + "input_white_noise.h5"
input_seed_phase_file = simdir + "seed"

L0 = L1 = L2 = args.L
N0=args.N
Np0=args.Np
Npm0=args.Npm
RedshiftLPT=args.RedshiftLPT
ai=z2a(RedshiftLPT)
RedshiftFCs=args.RedshiftFCs
af=z2a(RedshiftFCs)


common_params = {
    "Np0": Np0,
    "N0": N0,
    "L0": L0,
    "corner0": -L0/2,
    "corner1": -L1/2,
    "corner2": -L2/2,
    "h": cosmo["h"],
    "Omega_m": cosmo["Omega_m"],
    "Omega_b": cosmo["Omega_b"],
    "n_s": cosmo["n_s"],
    "sigma8": cosmo["sigma8"],
}

lpt_params = common_params.copy()
lpt_params["method"] = "lpt"
lpt_params["InputPowerSpectrum"] = input_power_file
lpt_params["InputWhiteNoise"] = input_white_noise_file

pm_params = common_params.copy()
pm_params["method"] = "pm"
pm_params["TimeStepDistribution"] = 0
pm_params["ai"] = ai
pm_params["af"] = af
pm_params["RedshiftLPT"] = RedshiftLPT
pm_params["RedshiftFCs"] = RedshiftFCs
pm_params["Npm0"] = Npm0

cola_params = common_params.copy()
cola_params["method"] = "cola"
cola_params["TimeStepDistribution"] = 0
cola_params["ai"] = ai
cola_params["af"] = af
cola_params["RedshiftLPT"] = RedshiftLPT
cola_params["RedshiftFCs"] = RedshiftFCs
cola_params["Npm0"] = Npm0

bullfrog_params = common_params.copy()
bullfrog_params["method"] = "bullfrog"
bullfrog_params["TimeStepDistribution"] = 0
bullfrog_params["ai"] = ai
bullfrog_params["af"] = af
bullfrog_params["RedshiftLPT"] = RedshiftLPT
bullfrog_params["RedshiftFCs"] = RedshiftFCs
bullfrog_params["Npm0"] = Npm0


if not isfile(input_white_noise_file):
    ## Get the seed for the white noise
    if args.seed is None:
        seed=np.random.randint(10**6)
    else:
        seed=args.seed

    print(">> Generating the initial white noise...\n")

    generate_white_noise_Field(seedphase=seed,
                            fname_whitenoise=input_white_noise_file,
                            seedname_whitenoise=input_seed_phase_file,
                            force_phase=False,
                            N=N0,
                            L=L0,
                            corner=-L0/2)


print("\n> Generating simulations cards...\n")


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
    

print("\n> Generating input power spectrum and ss file...\n")
    
# If cosmo["WhichSpectrum"] == "class", then the module classy is
# required.
if not isfile(input_power_file):
    Pk = PowerSpectrum(L0, L1, L2, N0, N0, N0, cosmo_small_to_full_dict(cosmo))
    Pk.write(input_power_file)
    
kmin = 0.0
kmax = 1.0
Nk = 64
k_modes = np.linspace(kmin, kmax, Nk + 1)[:-1]
input_ss_file = simdir + "input_ss_k_grid.h5"
if not isfile(input_ss_file):
    Gk = FourierGrid(L0, L1, L2, N0, N0, N0, k_modes=k_modes, kmax=kmax)
    Gk.write(input_ss_file)
    

print("\n> Starting LPT simulation...\n")

if not isfile(ICs_path) or not isfile(simdir + "lpt_density.h5") or not isfile(simdir + "lpt_particles.gadget3"):
    if not args.verbose_smby:
        with io.BytesIO() as f:
            with redirect_stdout(f):
                pySbmy(f"{wd}example_lpt.sbmy", f"{logdir}lpt.txt")
    else:
        pySbmy(f"{wd}example_lpt.sbmy", f"{logdir}lpt.txt")

print("\n> Starting PM simulations...\n")
for i in range(nsim_pm):
    file_ext = f"nsteps{nsteps_pm_list[i]}"  # "pm" is already in the filename
    if not isfile(simdir + f"{file_ext}_final_density_pm.h5"):
        print(f"\n>> Starting PM simulation {i}: {file_ext}\n")
        if not args.verbose_smby:
            with io.BytesIO() as f:
                with redirect_stdout(f):
                    pySbmy(f"{wd}{file_ext}_example_pm.sbmy", f"{logdir}{file_ext}_pm.txt")
        else:
            pySbmy(f"{wd}{file_ext}_example_pm.sbmy", f"{logdir}{file_ext}_pm.txt")
    else:
        print(f"\n>> PM simulation {i}: {file_ext} already completed\n")

print("\n> Starting COLA simulations...\n")
for i in range(nsim_cola):
    file_ext = f"nsteps{nsteps_cola_list[i]}"  # "cola" is already in the filename
    if not isfile(simdir + f"{file_ext}_final_density_cola.h5"):
        print(f"\n>> Starting COLA simulation {i}: {file_ext}\n")
        if not args.verbose_smby:
            with io.BytesIO() as f:
                with redirect_stdout(f):
                    pySbmy(f"{wd}{file_ext}_example_cola.sbmy", f"{logdir}{file_ext}_cola.txt")
        else:
            pySbmy(f"{wd}{file_ext}_example_cola.sbmy", f"{logdir}{file_ext}_cola.txt")
    else:
        print(f"\n>> COLA simulation {i}: {file_ext} already completed\n")

print("\n> Starting BullFrog simulations...\n")
for i in range(nsim_bullfrog):
    file_ext = f"nsteps{nsteps_bullfrog_list[i]}"  # "bullfrog" is already in the filename
    if not isfile(simdir + f"{file_ext}_final_density_bullfrog.h5"):
        print(f"\n>> Starting BullFrog simulation {i}: {file_ext}\n")
        if not args.verbose_smby:
            with io.BytesIO() as f:
                with redirect_stdout(f):
                    pySbmy(f"{wd}{file_ext}_example_bullfrog.sbmy", f"{logdir}{file_ext}_bullfrog.txt")
        else:
            pySbmy(f"{wd}{file_ext}_example_bullfrog.sbmy", f"{logdir}{file_ext}_bullfrog.txt")
    else:
        print(f"\n>> BullFrog simulation {i}: {file_ext} already completed\n")


print("\n> Done!\n")
