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

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cmocean

fs = 18
fs_titles = fs + 2
fs_legend = fs - 4
cmap = cmocean.cm.thermal
cols = ["C{}".format(i) for i in range(10)]




def plot_mode_runs_full_comparison(ax,
                                   Pk,
                                   Pk_ref,
                                   k,
                                   ts_keys:list[str],
                                   modes:dict[str,list[str]],
                                   nsteps_all:dict[str,dict[str,list[int]]],
                                   line_styles:list[str]|None=None,
                                   markers:list[str]|None=None,
                                   title:str|None=None,
                                   xlims:tuple[float,float]|None=None,
                                   ylims:tuple[float,float]|None=None,
                                   yticks:list[float]|None=None,
                                   bnd1:float|None=None,
                                   bnd2:float|None=None,
                                   y_log:bool=False,
                                   nstep_min:int=-1,
                                   lintresh:float=1e-3,
                                   cross:bool=False,):
    
    ax.set_xscale("log")
    if xlims is not None:
        ax.set_xlim(xlims)

    if y_log:
        ax.set_yscale("symlog", linthresh=lintresh)
    else:
        if ylims is not None:
            ax.set_ylim(ylims)
        if yticks is not None:
            ax.set_yticks(yticks)
    
    if bnd1 is not None:
        ax.axhspan(1 - bnd1, 1 + bnd1, color="grey", alpha=0.2)
    if bnd2 is not None:
        ax.axhspan(1 - bnd2, 1 + bnd2, color="grey", alpha=0.1)

    ax.grid(which="major",alpha=0.5)
    ax.grid(which="minor",alpha=0.2)

    ax.set_xlabel("$k$ [$h/\\mathrm{Mpc}$]", fontsize=fs)
    if not cross:
        ax.set_ylabel("$P(k)/P_\\mathrm{ref}(k)$" + (" - 1 " if y_log else ""), fontsize=fs)
    else:
        ax.set_ylabel("Cross correlation" + (" - 1 " if y_log else ""), fontsize=fs)

    ax.tick_params(which="both", direction="in")
    ax.tick_params(axis="both", which="major", labelsize=fs)
    ax.tick_params(axis="both", which="minor", labelsize=fs)

    if line_styles is None:
        line_styles={"PM":"--", "COLA":":", "BullFrog":"-", "INI":"-", "LPT":"-."}

    if markers is None:
        markers=dict(zip(ts_keys,["o","s","^","*","+","x","D","v","<",">"]))

    max_ts = max([sum([(n>=nstep_min) for n in nsteps_all[mode][key]]) for key in ts_keys for mode in modes[key]])
    


    for key in ts_keys:
        for mode in modes[key]:
            for i, nsteps in enumerate(nsteps_all[mode][key]):
                if nsteps>=nstep_min:
                    ax.plot(
                        k,
                        Pk[mode][key][i] / Pk_ref - y_log,
                        label=f"{key} {mode} $n_t={nsteps}$",
                        linestyle=line_styles[mode],
                        marker=markers[key],
                        markersize=4,
                        color=cols[i],
                    )
                

    _handles, _labels = plt.gca().get_legend_handles_labels()     

    handles, labels = [],[]
    i=0
    for key in ts_keys:
        for mode in modes[key]:
            j=i
            for nsteps in nsteps_all[mode][key]:
                if nsteps>=nstep_min:
                    handles.append(_handles[i])
                    labels.append(_labels[i])
                    i+=1
            handles+=[mpatches.Patch(color="none", label="") for _ in range(max(0, max_ts-(i-j)))]
            labels+=["" for _ in range(max(0, max_ts-(i-j)))]

    if title is not None:
        ax.set_title(title, fontsize=fs_titles)

    ax.legend(
        handles,
        labels,
        loc="upper center",
        ncol=len([mode for key in ts_keys for mode in modes[key]]),
        bbox_to_anchor=(0.5, -0.15),
        fontsize=fs_legend,
    )