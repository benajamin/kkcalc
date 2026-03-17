"""
An example script showing the KK transform inverse process.
"""

if __name__ == "__main__":
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from kkcalc2.models import (
        asf_im,
        asp_db_re,
        asp_db_re_extended,
    )
    from kkcalc2 import stoichiometry

    POLYSTYRENE = "CH"
    PS_NAME = "Polystyrene"
    ps_stoich = stoichiometry(POLYSTYRENE)

    # Import raw data
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    data_file = os.path.normpath(os.path.join(data_dir, "PS_004_-dc.txt"))
    data_PS = np.genfromtxt(data_file, skip_header=4)

    # Create the atomic scattering factors from NEXAFS data
    asf_PS: asf_im = asf_im.from_NEXAFS(
        energies=data_PS[:, 0],
        NEXAFS=data_PS[:, 1],
        stoichiometry=ps_stoich,
        scale_to_database=True,
    )

    # Also collect the real database values
    asp_PS_realdb = asp_db_re(stoichiometry=ps_stoich)

    # Calculate the real
    asp_PS_real = asf_PS.kk_transform(improve_accuracy=True, max_iter=5, tolerance=1e-1)

    # Extend the real
    asp_PS_real_extended = asp_db_re_extended(data_asf=asf_PS, database=ps_stoich)

    # Plot the objects
    fig, axs = plt.subplots(2, 1, figsize=(16, 9))
    axs: list[plt.Axes]

    axs[0].plot(
        asf_PS.energies,
        asf_PS.factors,
        marker=".",
        label=f"Imag, Num Points: {asf_PS.energies.shape[0]}",
    )
    axs[1].plot(
        asp_PS_realdb.energies,
        asp_PS_realdb.asf,
        marker=".",
        label=f"Real DB, Num Points: {asp_PS_realdb.energies.shape[0]}",
    )
    axs[1].plot(
        asp_PS_real.energies,
        asp_PS_real.factors,
        marker=".",
        label=f"Real KK, Num Points: {asp_PS_real.energies.shape[0]}",
    )
    axs[1].plot(
        asp_PS_real_extended.energies,
        asp_PS_real_extended.asf,
        marker=".",
        label=f"Real KK Ext., Num Points: {asp_PS_real_extended.energies.shape[0]}",
    )

    # Finalise
    axs[0].set_ylabel("Imaginary")
    axs[0].set_xlabel("Energy (eV)")
    axs[0].set_xscale("log")
    axs[0].legend()
    axs[1].set_ylabel("Real")
    axs[1].set_xlabel("Energy (eV)")
    axs[1].set_xscale("log")
    axs[1].legend()
    plt.show()
