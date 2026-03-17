"""
An example demonstrating the use of `improve_accuracy` in kk_transform.

This example demonstrates how to use the `improve_accuracy` option in the `kk_transform` method
of the `asp_db_extended` model to improve the accuracy of the Kramers-Kronig transform at a K-Edge.
The dataset is the Polystryene Carbon K-Edge.
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
    asp_db_PS = asp_db_im(ps_stoich)

    # Import Data
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    data_file = os.path.normpath(os.path.join(data_dir, "PS_004_-dc.txt"))
    data_PS = np.genfromtxt(data_file, skip_header=4)
    assert data_PS.shape[1] == 2, "Data file must have two columns"

    # Create the atomic scattering factors from NEXAFS data
    asf_PS: asf_im = asf_im.from_NEXAFS(
        energies=data_PS[:, 0], NEXAFS=data_PS[:, 1], stoichiometry=ps_stoich
    )

    # Combine the data with the database into polynomials
    PS_imag = asp_db_im_extended(
        data_asf=asf_PS, database=asp_db_PS, merge_domain=(280, 320)
    )

    # Setup plotting to demonstrate
    fig: plt.Figure
    ax: list[plt.Axes]
    fig, axs = plt.subplots(2, 2, figsize=(16, 9), sharex=True)
    for ax in axs:
        ax[0].set_ylabel("Real")
        ax[1].set_ylabel("Imag")
        ax[0].set_xlabel("Energy")
        ax[1].set_xlabel("Energy")
        ax[0].set_xscale("log")
        ax[1].set_xscale("log")

    # Instead of using the extended database, use the data to observe what the kk transform does.
    # PS_imag = asf_PS.to_ASP()

    # Use the KK transform on two sets of objects, the extended data, and the raw data.
    for j, asp_obj in enumerate([PS_imag, asf_PS.to_ASP()]):
        past_points = []
        ax = axs[j]

        for i in range(0, 5):
            # Perform the transfom with i improvement iterations
            PS_extended_real = asp_obj.kk_transform(
                relativistic_correction=ps_stoich.relativistic_correction,
                tolerance=0.01,
                max_iter=i,
            )

            # Convert the imaginary coefs to imaginary factors
            PS_imag_factors = PS_imag.eval_asf(PS_extended_real.energies)

            # Find new points
            new_energies = np.setdiff1d(
                PS_extended_real.energies, past_points, assume_unique=True
            )
            new_points = np.array(
                [
                    i
                    for i, e in enumerate(PS_extended_real.energies)
                    if e in new_energies
                ]
            )
            past_points = np.concatenate((past_points, new_energies))

            # Plot a line of the unimproved data
            if i == 0:
                ax[0].plot(
                    PS_extended_real.energies,
                    PS_extended_real.factors,
                    alpha=0.4,
                    label="Real Unimproved",
                )
                ax[1].plot(
                    PS_extended_real.energies,
                    PS_imag_factors,
                    alpha=0.4,
                    label="Imag Unimproved",
                )

            # Scatter plot the new points
            l1 = ax[0].scatter(
                PS_extended_real.energies[new_points],
                PS_extended_real.factors[new_points],
                s=(i + 1) ** 2,
                label=f"Real Iteration {i}, {len(new_points)} points",
            )
            l2 = ax[1].scatter(
                PS_extended_real.energies[new_points],
                PS_imag_factors[new_points],
                s=(i + 1) ** 2,
                label=f"Imag Iteration {i}, {len(new_points)} points",
            )

        # Add a line plot of the final result
        ax[0].plot(
            PS_extended_real.energies,
            PS_extended_real.factors,
            alpha=0.4,
            c=l1.get_facecolor(),
            label="Real Final",
        )
        ax[1].plot(
            PS_extended_real.energies,
            PS_imag_factors,
            alpha=0.4,
            c=l2.get_facecolor(),
            label="Imag Final",
        )

        # Set x-axis limits to the K-Edge
        lims = [270, 330]
        ax[0].set_xlim(*lims)
        ax[1].set_xlim(*lims)

        # Finalize UI elements
        ax[0].legend()
        ax[1].legend()

    fig.suptitle(
        "Polystyrene Carbon K-Edge, Extended Database vs Raw Data KK Transforms"
    )
    axs[0][0].set_title("Extended Transform")
    axs[0][1].set_title("Extended Imag")
    axs[1][0].set_title("Raw Transform")
    axs[1][1].set_title("Raw Imag")

    fig.tight_layout()
    plt.show()
