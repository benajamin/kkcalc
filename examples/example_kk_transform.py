"""
This example demonstrates how to perform a simple Kramers-Kronig transform on a dataset.
"""

from kkcalc2.models import (
    asf_im,
    asp_db_im,
    asp_db_im_extended,
)
from kkcalc2 import stoichiometry
import numpy as np
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    # This file is part of the Kramers-Kronig Calculator software package.

    # Create a merge of physical data and database data
    POLYSTYRENE = "CH"
    PS_NAME = "Polystyrene"
    ps_stoich = stoichiometry(POLYSTYRENE)
    asp_db_PS = asp_db_im(ps_stoich, name="PS Database")

    # Import Data
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    data_file = os.path.normpath(os.path.join(data_dir, "PS_004_-dc.txt"))
    data_PS = np.genfromtxt(data_file, skip_header=4)
    assert data_PS.shape[1] == 2, "Data file must have two columns"

    # Create the atomic scattering factors from NEXAFS data
    asf_PS: asf_im = asf_im.from_NEXAFS(
        energies=data_PS[:, 0],
        NEXAFS=data_PS[:, 1],
        stoichiometry=ps_stoich,
        name=PS_NAME,
    )

    # Combine the data with the database into polynomials
    PS_imag = asp_db_im_extended(
        data_asf=asf_PS, database=asp_db_PS, merge_domain=(280, 320)
    )

    # Setup plotting to demonstrate
    fig: plt.Figure
    ax0: plt.Axes
    ax1: plt.Axes
    fig, ax0 = plt.subplots(1, 1)
    ax1 = ax0.twinx()
    ax0.set_ylabel("Real")
    ax1.set_ylabel("Imag")
    ax0.set_xlabel("Energy (eV)")
    # ax1.set_xlabel('Energy')
    ax0.set_xscale("log")
    # ax1.set_xscale('log')

    # Perform the transfom with i improvement iterations
    PS_real = PS_imag.kk_transform(stoichiometry=ps_stoich, tolerance=0.01)

    # Add a line plot of the extended and transformed results
    l1 = ax0.plot(
        PS_real.energies,
        PS_real.factors,
        marker=".",
        c="teal",
        label=f"Real, Num Points: {PS_real.energies.shape[0]}",
    )
    l2 = ax1.plot(
        PS_imag.energies,
        PS_imag.asf,
        marker=".",
        c="orange",
        label=f"Imag, Num Points: {PS_imag.energies.shape[0]}",
    )
    # Scale the data to the same range
    idx_overlap = (PS_imag.energies >= asf_PS.energies.min()) & (
        PS_imag.energies <= asf_PS.energies.max()
    )
    data_y_scaled = (PS_imag.asf.max() - PS_imag.asf.min()) / (
        asf_PS.factors.max() - asf_PS.factors.min()
    ) * (asf_PS.factors - asf_PS.factors.min()) + PS_imag.asf.min()
    l0 = ax1.plot(
        asf_PS.energies,
        data_y_scaled,
        marker=".",
        alpha=0.4,
        c="black",
        label=f"Imag Raw, Num Points: {asf_PS.energies.shape[0]}",
    )

    # Limit the plot to the range of data
    domain = (asf_PS.energies.min(), asf_PS.energies.max())
    ax0.set_xlim(*domain)

    # Finalize UI elements
    lines = l1 + l2 + l0
    labels = [line.get_label() for line in lines]
    ax0.legend(lines, labels)
    fig.tight_layout()
    plt.show()
