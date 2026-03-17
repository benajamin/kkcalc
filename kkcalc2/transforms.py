"""
This module contains the Kramers-Kronig transform methods.
"""

import math
import numpy as np
import numpy.typing as npt
import warnings

DEF_ITER: int = 50
"""The default number of iterations to use in improving the accuracy of the Kramers-Kronig transform."""
DEF_TOL: float = 1e-2
"""The default tolerance to use in improving the accuracy of the Kramers-Kronig transform."""

from kkcalc2 import conversions  # noqa E402


def KK_General_PP(
    target_energies: npt.NDArray,
    energies: npt.NDArray,
    imag_coefs: npt.NDArray,
    orders: npt.NDArray | None = None,
    relativistic_correction: float = 0,
) -> npt.NDArray:
    r"""
    Apply the Kramers-Kronig transform on general coefficients.

    Converts from `f1` defined with general coefficients `orders` to `f2`,
    with the 'Piecewise Polynomial' algorithm by Watts et. al. (2014).

    .. math::
        f_2 (E) = \frac{2}{\pi} P \int_{0}^{\infty}\frac{x f_1(x)}{x^2 - E^2} dx + \mathcal{Z}^\star

    Parameters
    ----------
    target_energies : array_like
        An array of energies with length `L` at which to evaluate the real spectrum.
    energies : array_like
        An array of energies with length `M+1` describing the spans on which the `imag_coefs` are defined.
    imag_coefs : array_like
        A 2D array with shape `(M, N)`, consisting `M` sets of `N` polynomial coefficients,
        belonging to the power terms indicated by 'order'.
    orders : array_like, optional
        A list of length `N`, listing orders corresponding to the imag_coef indices.
        By default assumes [1, 0, -1, ...] for the columns of imaginary_spectrum.
    relativistic_correction : float
        The relativistic correction to the Kramers-Kronig transform.
        .. math::
            \mathcal{Z}^\star = \sum_i (Z_i - (Z_i/82.5)^{2.37}) \cdot n_i
        Calculated via the `stoich` module, either by the  `relativistic_correction_eq`
        or the `stoichiometry` object property `stoichiometry.relativistic_correction`.

    Returns
    -------
    np.ndarray
        The real component of the atomic scattering factors, evaluated at `target_energies`.
    """
    # Convert inputs to numpy arrays is not already
    target_energies = np.asarray(target_energies)
    energies = np.asarray(energies)
    imag_coefs = np.asarray(imag_coefs)
    if orders is None:
        orders = np.linspace(1, 2 - imag_coefs.shape[1], imag_coefs.shape[1], dtype=int)
    else:
        orders = np.asarray(orders)

    # Check input dimensions match
    if imag_coefs.shape[0] != energies.shape[0] - 1:
        raise ValueError(
            f"First axis of imaginary_spectrum must have one less row ({imag_coefs.shape[0]}) than energies ({energies.shape[0]})"
        )
    if imag_coefs.shape[1] != orders.shape[0]:
        raise ValueError(
            f"Second axis of imaginary_spectrum must have the same number of columns ({imag_coefs.shape[1]}) as orders ({orders.shape[0]})"
        )

    # Need to build arrays with dimensions X-E-n [Energies, Evaluation Energies, Orders]
    C = np.tile(imag_coefs[:, np.newaxis, :], (1, len(target_energies), 1))  # 2D to 3D
    N = np.tile(
        orders[np.newaxis, np.newaxis, :], (len(energies) - 1, len(target_energies), 1)
    )  # 1D to 3D
    X = np.tile(
        energies[:, np.newaxis, np.newaxis], (1, len(target_energies), len(orders))
    )  # 1D to 3D
    E = np.tile(
        target_energies[np.newaxis, :, np.newaxis], (len(energies) - 1, 1, len(orders))
    )  # 1D TO 3D
    # Calculate when evaluating energies are equal to the data energies, in a boolean array (1 or 0)
    poles = np.equal(
        X,
        np.tile(
            target_energies[np.newaxis, :, np.newaxis], (len(energies), 1, len(orders))
        ),
    )

    # Basic integral, the resulting shape matches the evaluation energies
    integral: npt.NDArray = np.sum(
        -C
        * (-E) ** N
        * np.log(np.abs((X[1:, :, :] + E) / (X[:-1, :, :] + E)))  # X_{i+1} / X_i
        - C
        * (E) ** N
        * (
            1 - poles[1:, :, :]
        )  # If a pole in numerator or denominator, then zero contribution
        # Need to multiply by something to prevent contribution at i=0.
        # * np.array([0] + [1] * (poles.shape[0]-2))[:, np.newaxis, np.newaxis] #* (1-poles[:-1,:,:])
        * np.log(
            np.abs(
                # `(X_{i+1}-E) / (X_i-E)` if not a pole, else `(X_{i+1}-E) / (X_{i-1}-E)`
                (
                    X[1:, :, :] - E + poles[1:, :, :]
                )  # Non-zero if a pole, just to avoid log(0)
                / (
                    # Alternative
                    X[:-1, :, :] * (1 - poles[:-1, :, :])  # If not a pole, then use X
                    + poles[:-1, :, :]
                    * np.r_[
                        X[np.newaxis, 0], X[0:-2]
                    ]  # If a pole, then use X-1 value prior (unless at X[0]).
                    - E
                    # + np.r_[(poles[0] * poles[1])[np.newaxis,:], np.zeros(poles[3:,].shape)] # If poles[0] or poles[1] is a pole, then add small amount.
                    # Original
                    # X[:-1,:,:] * (1-poles[:-1,:,:]) # If not a pole, then use X
                    # + poles[:-1,:,:] * X[[0, *range(energies.shape[0]-2)],:,:] #If a pole, then use X-1 value prior (unless at X[0]).
                    # - E
                    # TODO Prevent divide by zero at i = 0
                    # + np.array([1] + [0] * (poles.shape[0]-2))[:, np.newaxis, np.newaxis] # Prevent divide by zero at pole at i = 0
                )
            )
        )
        * (
            1
            - np.r_[
                (poles[0] * poles[1])[np.newaxis, :],
                np.zeros(poles[2:,].shape),
            ]
        ),  # If poles[0] or poles[1] is a pole, zero value.
        axis=(0, 2),
    )

    ### Calculate the Kramers-Kronig integral additional terms
    if np.any(orders <= -2):  # N<=-2, ln(x) terms
        i = (slice(None, None, None), slice(None, None, None), orders <= -2)
        integral += np.sum(
            C[*i]
            * ((-E[*i]) ** N[*i] + E[*i] ** N[*i])
            * np.log(np.absolute((X[1:, :, orders <= -2]) / (X[:-1, :, orders <= -2]))),
            axis=(0, 2),
        )

    if np.any(orders >= 0):  # N>=0,  x^k terms
        for ni in np.where(orders >= 0)[0]:
            i = [slice(None, None, None), slice(None, None, None), ni]
            n = orders[ni]
            for k in range(n, 0, -2):
                integral += np.sum(
                    C[*i]
                    / float(-k)
                    * 2
                    * E[*i] ** (n - k)
                    * (X[1:, :, ni] ** k - X[:-1, :, ni] ** k),
                    axis=0,
                )

    if np.any(orders <= -3):  # N<=-3, x^k terms
        for ni in np.where(orders <= -3)[0]:
            i = [slice(None, None, None), slice(None, None, None), ni]
            n = orders[ni]
            for k in range(n + 2, 0, 2):
                n = n.astype(float)
                integral += np.sum(
                    C[*i]
                    / float(k)
                    * ((-1) ** (n - k) + 1)
                    * E[*i] ** (n - k)
                    * (X[1:, :, ni] ** k - X[:-1, :, ni] ** k),
                    axis=0,
                )
    return integral / math.pi + relativistic_correction


def KK_General_PP_inv(
    target_energies: npt.NDArray,
    energies: npt.NDArray,
    real_coefs: npt.NDArray,
    relativistic_correction: float,
    orders: npt.NDArray | None = None,
) -> npt.NDArray:
    r"""
    Apply the Kramers-Kronig transform on general coefficients.

    Converts from `f2` defined with general coefficients `orders` to `f1`,
    with the 'Piecewise Polynomial' algorithm by Watts et. al. (2014).
    Uses the `KK_General_PP` algorithm to calculate the inverse Kramers-Kronig transform.

    .. math::
        f_1 (E) = Z^\star - \frac{2}{\pi} P \int_{0}^{\infty}\frac{x f_2(x)}{x^2 - E^2} dx

    Parameters
    ----------
    target_energies : array_like
        An array of energies with length `L` at which to evaluate the imaginary spectrum.
    energies : array_like
        An array of energies with length `M+1` describing the spans on which the `real_coefs` are defined.
    real_coefs : array_like
        A 2D array with shape `(M, N)`, consisting `M` sets of `N` polynomial coefficients,
        belonging to the power terms indicated by 'order'.
    relativistic_correction : float
        The relativistic correction to the Kramers-Kronig transform.
        .. math::
            \mathcal{Z}^\star = \sum_i (Z_i - (Z_i/82.5)^{2.37}) \cdot n_i
        Calculated via the `stoich` module, either by the  `relativistic_correction_eq`
        or the `stoichiometry` object property `stoichiometry.relativistic_correction`.
    orders : array_like, optional
        A list of length `N`, listing orders corresponding to the real_coef indices.
        By default assumes [1, 0, -1, ...] for the columns of real_spectrum.

    Returns
    -------
    np.ndarray
        Evaluated real coefficients of the atomic scattering factors at `target_energies`.
    """
    # Use the
    return (
        -target_energies
        * KK_General_PP(
            target_energies=target_energies,
            energies=energies,
            # Moves coefficients one place to the right (equivalent to moving orders one to the left).
            # This implies a change in coefficient with energy order.
            # The new energy order is instead seen as [0, -1, -2, -3, 1] when multiplying by coefficients.
            # TODO: Why?
            imag_coefs=np.roll(real_coefs, 1, axis=1),
            orders=orders,
            relativistic_correction=-relativistic_correction,  # inverse the relativistic correction.
        )
    )


def KK_PP(
    target_energies: npt.NDArray,
    energies: npt.NDArray,
    imag_coefs: npt.NDArray,
    relativistic_correction: float,
) -> npt.NDArray:
    r"""
    Apply the Kramers-Kronig transform on imaginary polynomials (f2) to calculate real factors (f1).

    Converts from `f2` to `f1` with 'Piecewise Polynomial' algorithm by Watts et. al. (2014).

    .. math::
        f_2 (E) = \frac{2}{\pi} P \int_{0}^{\infty}\frac{x f_1(x)}{x^2 - E^2} dx + \mathcal{Z}^\star

    Parameters
    ----------
    target_energies : array_like
        An array of energies with length `L` at which to evaluate the real spectrum.
    energies : array_like
        An array of energies with length `M+1` describing the spans on which the `imag_coefs` are defined.
    imag_coefs : array_like
        A 2D array with shape `(M, 5)`, consisting `M` sets of 5 polynomial coefficients
        for the imaginary part of the scattering factors defined between the `M+1` energies.
        The 5 coefficients correspond to energy powers [1, 0, -1, -2, -3].
    relativistic_correction : float
        The relativistic correction to the Kramers-Kronig transform.
        .. math::
            \mathcal{Z}^\star = \sum_i (Z_i - (Z_i/82.5)^{2.37}) \cdot n_i
        Calculated via the `stoich` module, either by the  `relativistic_correction_eq`
        or the `stoichiometry` object property `stoichiometry.relativistic_correction`.

    Returns
    -------
    np.ndarray
        The real part of the scattering factors evaluated at `target_energies`.
    """
    target_energies = np.asarray(target_energies)
    energies = np.asarray(energies)
    imag_coefs = np.asarray(imag_coefs)

    # if np.all(target_energies == energies):
    #     # If every target energy is already in the energy list, ... TODO
    #     raise NotImplementedError(
    #         "This function is not yet implemented for the case where every target energy is in the energy list."
    #     )
    # else:
    # M is the number of polynomial energy spans, N is the number of target energies.

    X1 = energies[0:-1]  # M
    X2 = energies[1:]  # M
    E = np.tile(
        target_energies, (len(energies) - 1, 1)
    ).T  # Results in a 2D of shape (N, M)
    coefs_T = imag_coefs.T
    #
    Symb_1 = (
        (coefs_T[0, :] * E + coefs_T[1, :]) * (X2 - X1)
        + 0.5 * coefs_T[0, :] * (X2**2 - X1**2)
        - (coefs_T[3, :] / E + coefs_T[4, :] * E**-2) * np.log(np.abs(X2 / X1))
        + coefs_T[4, :] / E * (X2**-1 - X1**-1)
    )
    #
    Symb_2 = (
        (-coefs_T[0, :] * E + coefs_T[1, :]) * (X2 - X1)
        + 0.5 * coefs_T[0, :] * (X2**2 - X1**2)
        + (coefs_T[3, :] / E - coefs_T[4, :] * E**-2) * np.log(np.abs(X2 / X1))
        - coefs_T[4, :] / E * (X2 ** (-1) - X1 ** (-1))
    ) + (
        coefs_T[0, :] * E**2
        - coefs_T[1, :] * E
        + coefs_T[2, :]
        - coefs_T[3, :] * E**-1
        + coefs_T[4, :] * E**-2
    ) * np.log(np.abs((X2 + E) / (X1 + E)))
    #
    Symb_3 = (
        (1 - 1 * ((X2 == E) | (X1 == E)))
        * (
            coefs_T[0, :] * E**2
            + coefs_T[1, :] * E
            + coefs_T[2, :]
            + coefs_T[3, :] * E**-1
            + coefs_T[4, :] * E**-2
        )
        * np.log(np.abs((X2 - E + 1 * (X2 == E)) / (X1 - E + 1 * (X1 == E))))
    )
    # Sum areas for approximate integral
    Symb_B = np.sum(Symb_2 - Symb_1 - Symb_3, axis=1)

    # Patch Poles
    poles = energies[1:-1] == E[:, 0:-1]
    E_sing = np.append(np.insert(np.any(poles, axis=0), [0, 0], False), [False, False])
    Eval_sing = np.any(poles, axis=1)

    X1 = energies[E_sing[2:]]
    XE = energies[E_sing[1:-1]]
    X2 = energies[E_sing[:-2]]
    # C1 = Full_coeffs[:, E_sing[2:-1]] # Not used... why?
    C2 = coefs_T[:, E_sing[1:-2]]
    Symb_singularities = np.zeros(len(target_energies))
    val = (
        C2[0, :] * XE**2
        + C2[1, :] * XE
        + C2[2, :]
        + C2[3, :] * XE**-1
        + C2[4, :] * XE**-2
    ) * np.log(np.abs((X2 - XE) / (X1 - XE)))
    # print(val)
    # print(Eval_sing)
    Symb_singularities[Eval_sing] = val
    # Finish things off
    KK_Re = (Symb_B - Symb_singularities) / (
        math.pi * target_energies
    ) + relativistic_correction
    return KK_Re


def KK_PP_inv(
    target_energies: npt.NDArray,
    energies: npt.NDArray,
    real_coefs: npt.NDArray,
    relativistic_correction: float,
) -> npt.NDArray:
    r"""
    Apply the Kramers-Kronig transform on real polynomials (f1) to calculate imaginary factors (f2).

    Converts from `f1` to `f2` with 'Piecewise Polynomial' algorithm by Watts et. al. (2014).
    Uses the `KK_PP` algorithm to calculate the inverse Kramers-Kronig transform.

    .. math::
        f_2 (E) = \frac{2}{\pi} P \int_{0}^{\infty}\frac{x f_1(x)}{x^2 - E^2} dx + \mathcal{Z}^\star

    Parameters
    ----------
    target_energies : array_like
        An array of energies with length `L` at which to evaluate the real spectrum.
    energies : array_like
        An array of energies with length `M+1` describing the spans on which the `imag_coefs` are defined.
    real_coefs : array_like
        A 2D array with shape `(M, 5)`, consisting `M` sets of 5 polynomial coefficients
        for the real part of the scattering factors defined between the `M+1` energies.
        The 5 coefficients correspond to energy powers [1, 0, -1, -2, -3].
    relativistic_correction : float
        The relativistic correction to the Kramers-Kronig transform.
        .. math::
            \mathcal{Z}^\star = \sum_i (Z_i - (Z_i/82.5)^{2.37}) \cdot n_i
        Calculated via the `stoich` module, either by the  `relativistic_correction_eq`
        or the `stoichiometry` object property `stoichiometry.relativistic_correction`.

    Returns
    -------
    np.ndarray
        The real part of the scattering factors evaluated at `target_energies`.
    """
    ## Inverse KK is only a minor modification of the forward algorithm
    return -target_energies * KK_PP(
        target_energies=target_energies,
        energies=energies,
        imag_coefs=np.roll(real_coefs, 1, axis=1),
        relativistic_correction=-relativistic_correction,
    )


def improve_accuracy(
    energies: npt.ArrayLike,
    real_asf: npt.ArrayLike,
    imag_coefs: npt.ArrayLike,
    relativistic_correction: float,
    tolerance: float = DEF_TOL,
    max_iter: int = DEF_ITER,
    orders: npt.NDArray[np.integer] | None = None,
) -> tuple[npt.NDArray, npt.NDArray]:
    r"""
    Calculate extra data points so that the Kramers-Kronig transform is more accurate.

    Parameters
    ----------
    energies : npt.NDArray
        The photon energies with length `N`.
    real_asf : npt.NDArray
        The real part of the scattering factors with length `N`.
    imag_coefs : npt.NDArray
        Polynomial coefficients of shape `(N+1, M)` representing the imaginary
        part of the scattering factors, where `M` is the number of coefficients.
    relativistic_correction : float
        The relativistic correction to the Kramers-Kronig transform.
        .. math::
            \mathcal{Z}^\star = \sum_i (Z_i - (Z_i/82.5)^{2.37}) \cdot n_i
        Calculated via the `stoich` module, either by the  `relativistic_correction_eq`
        or the `stoichiometry` object property `stoichiometry.relativistic_correction`.
    tolerance : float
        The level of error allowed in the linear extrapolation.
    max_iter : int, optional
        The maximum number of iterations an energy interval can be split.
        By default 50.
    orders : npt.NDArray, optional
        The polynomial orders corresponding to the columns of `imag_coefs`.
        Must be the same length (`M`). By default None, which assumes `[1, 0, -1, -2, ...]`.

    Returns
    -------
    energies : npt.NDArray
        The photon energies of length N.
    real : npt.NDArray
        The real part of the scattering factors, of length N.

    Notes
    -----
    Instability occurs near the low energy values, where more and more points are added to the data set.
    Setting a tolerance too small can result in errors in this region.
    """
    assert np.all(energies[:-1] <= energies[1:]), "Energies must be in ascending order."

    # TODO: Add a try except statement for when too many iterations occur,
    # and the midpoint calculation fails to find a midpoint value. For improve_acc_inv as well.

    energies = np.asarray(energies)
    real_asf = np.asarray(real_asf)
    imag_coefs = np.asarray(imag_coefs)
    if orders is not None:
        orders = np.asarray(orders)

    # List has N items, indexed from 0 to N-1
    # idx_extra is an array of indexes 1 to N-1 (representing the insertion index of every single midpoint between indexes 0 and N-1).
    # Midpoints are calculated as x_mid = (x[i] + x[i-1]) / 2.
    idx_extra = np.arange(1, energies.shape[0], dtype=int)  # 1, 2 ... N-1

    # Imag polynomial to factors
    imag = conversions.ASP_to_ASF(energies, imag_coefs, orders)

    # Iterate until the error is below the tolerance, or the maximum number of iterations is reached
    for i in range(max_iter):
        # Check if any energy delta is zero
        if np.any(energies[idx_extra] == energies[idx_extra - 1]):
            idx_extra = np.delete(
                idx_extra, np.where(energies[idx_extra] == energies[idx_extra - 1])
            )

        # Get energy midpoints
        en_mid = (energies[idx_extra] + energies[idx_extra - 1]) / 2
        # Ensure midpoints are not equal to existing energies
        if np.any(np.isin(en_mid, energies)):
            idx_extra = np.delete(idx_extra, np.where(np.isin(en_mid, energies)))
        # If no midpoints left, break the loop
        if idx_extra.shape[0] == 0:
            warnings.warn(
                f"Midpoints beyond precision at iteration '{i}'.",
                UserWarning,
            )
            break

        # Calculate new midpoint imag values
        im_mid = conversions.ASP_to_ASF(
            energies=en_mid,  # midpoint energy between i-1 and i uses coefs at i-1
            coefs=imag_coefs[idx_extra - 1, :],
            orders=orders,
        )

        # Calculate new midpoint real values
        re_mid = KK_PP(
            target_energies=en_mid,
            energies=energies,
            imag_coefs=imag_coefs,
            relativistic_correction=relativistic_correction,
        )

        # Evaluate new (polynomial) values to the average of the old (linear) values. If coefs are linear, this will be zero.
        # Difference from linear is the error, bigger is better for finding new corrections.
        im_err = np.abs(im_mid - (imag[idx_extra] + imag[idx_extra - 1]) / 2)
        re_err = np.abs(re_mid - (real_asf[idx_extra] + real_asf[idx_extra - 1]) / 2)

        # Boolean for improvement - newly evaluated points have a change greater than the tolerance
        improved = (re_err > tolerance) | (im_err > tolerance)

        # Manual override for the first midpoint index near 10 eV, which doesn't converge.
        if improved[0] and idx_extra[0] == 1 and i > 20:
            improved[0] = False

        # Check if at satisfactory level
        if np.sum(improved) == 0:
            # # Return values if no improvements are made
            return energies, real_asf  # , imag

        else:  # some improvements are made
            idx_improved = idx_extra[
                improved
            ]  # insertion indexes for midpoints where improvements are made

            # Insert new points and values where errors are big. Energies length becomes M = N + sum(improved)
            energies = np.insert(energies, idx_improved, en_mid[improved])
            imag = np.insert(imag, idx_improved, im_mid[improved])
            real_asf = np.insert(real_asf, idx_improved, re_mid[improved])

            ### Create new indexes and insert duplicate coefficients for the next iteration
            # Create new array of midpoint locations to evaluate.
            new_value_locs = np.insert(
                arr=np.zeros(
                    imag_coefs.shape[0], dtype=bool
                ),  # Copy existing midpoint list
                obj=idx_extra[
                    improved
                ],  # Add the locations where improvements were made
                values=True,  # Insert True
            )
            new_midpoint_locs = np.where(new_value_locs)[
                0
            ]  # Locate the indexes of new midpoint values after insertions.

            # Duplicate coefficients at the improved indexes
            for j in range(idx_improved.shape[0] - 1, -1, -1):
                # Iterate from the last improved coefficient to the start
                idx = (
                    idx_improved[j] - 1
                )  # Get the midpoint index, move back one to get the coefficient index that defines that region
                # Duplicate the coefficients at the index
                imag_coefs = np.r_[
                    imag_coefs[:idx], [imag_coefs[idx]], imag_coefs[idx:]
                ]  # Insert the new coefficients

            # Update the set of indexes where new midpoints need to be calculated.
            # Add the midpoint after (in addition to before) the new points.
            # Transpose required otherwise indexes out of order.
            idx_extra = np.vstack(
                (new_midpoint_locs, new_midpoint_locs + 1)
            ).T.flatten()

    # Return values if the maximum number of iterations is reached
    return energies, real_asf


def improve_accuracy_inv(
    energies: npt.ArrayLike,
    imag_asf: npt.ArrayLike,
    real_coefs: npt.ArrayLike,
    relativistic_correction: float,
    tolerance: float = DEF_TOL,
    max_iter: int = DEF_ITER,
    orders: npt.ArrayLike | None = None,
) -> tuple[npt.NDArray, npt.NDArray]:
    r"""
    Calculate extra data points so that the Kramers-Kronig transform is more accurate.

    Parameters
    ----------
    energies : npt.NDArray
        The photon energies with length `N`.
    imag_asf : npt.NDArray
        The imaginary part of the scattering factors with length `N`.
    real_coefs : npt.NDArray
        Polynomial coefficients of shape `(N+1, M)` representing the real
        part of the scattering factors, where `M` is the number of coefficients.
    relativistic_correction : float
        The relativistic correction to the Kramers-Kronig transform.
        .. math::
            \mathcal{Z}^\star = \sum_i (Z_i - (Z_i/82.5)^{2.37}) \cdot n_i
        Calculated via the `stoich` module, either by the  `relativistic_correction_eq`
        or the `stoichiometry` object property `stoichiometry.relativistic_correction`.
    tolerance : float
        The level of error allowed in the linear extrapolation.
    max_iter : int, optional
        The maximum number of iterations an energy interval can be split.
        By default 50.
    orders : npt.NDArray, optional
        The polynomial orders corresponding to the columns of `imag_coefs`.
        Must be the same length (`M`). By default None, which assumes `[1, 0, -1, -2, ...]`.

    Returns
    -------
    energies : npt.NDArray
        The photon energies of length N.
    imag : npt.NDArray
        The imaginary part of the scattering factors, of length N.

    Notes
    -----
    Instability occurs near the low energy values, where more and more points are added to the data set.
    Setting a tolerance too small can result in errors in this region.
    """
    assert np.all(energies[:-1] <= energies[1:]), "Energies must be in ascending order."

    energies = np.asarray(energies)
    imag_asf = np.asarray(imag_asf)
    real_coefs = np.asarray(real_coefs)
    if orders is not None:
        orders = np.asarray(orders)

    # List has N items, indexed from 0 to N-1
    # idx_extra is an array of indexes 1 to N-1 (representing the insertion index of every single midpoint between indexes 0 and N-1).
    # Midpoints are calculated as x_mid = (x[i] + x[i-1]) / 2.
    idx_extra = np.linspace(1, energies.shape[0] - 1, energies.shape[0] - 1, dtype=int)

    # # Imag polynomial to factors
    # real = conversions.ASP_to_ASF(energies, real_coefs, orders)

    # Iterate until the error is below the tolerance, or the maximum number of iterations is reached
    for i in range(max_iter):
        # Check if any energy delta is zero
        if np.any(energies[idx_extra] == energies[idx_extra - 1]):
            idx_extra = np.delete(
                idx_extra, np.where(energies[idx_extra] == energies[idx_extra - 1])
            )

        # Get energy midpoints
        en_mid = (energies[idx_extra] + energies[idx_extra - 1]) / 2

        # # Calculate new midpoint imag values
        # re_mid = conversions.ASP_to_ASF(
        #     energies=en_mid, # midpoint energy between i-1 and i uses coefs at i-1
        #     coefs=real_coefs[idx_extra-1, :],
        #     orders=orders
        # )

        # Calculate new midpoint real values
        im_mid = KK_PP_inv(
            target_energies=en_mid,
            energies=energies,
            real_coefs=real_coefs,
            relativistic_correction=relativistic_correction,
        )

        # Evaluate new (polynomial) values to the average of the old (linear) values. If coefs are linear, this will be zero.
        # Difference from linear is the error, bigger is better for finding new corrections.
        # re_err = np.abs(re_mid - (real[idx_extra] + real[idx_extra-1]) / 2)
        im_err = np.abs(im_mid - (imag_asf[idx_extra] + imag_asf[idx_extra - 1]) / 2)

        # Boolean for improvement - newly evaluated points have a change greater than the tolerance
        improved = im_err > tolerance  # | (re_err > tolerance)

        # Manual override for the first midpoint index near 10 eV, which doesn't converge.
        if improved[0] and idx_extra[0] == 1 and i > 20:
            improved[0] = False

        # Check if at satisfactory level
        if np.sum(improved) == 0:
            # # Return values if no improvements are made
            return energies, imag_asf  # , real_coefs

        else:  # some improvements are made
            idx_improved = idx_extra[
                improved
            ]  # insertion indexes for midpoints where improvements are made

            # Insert new points and values where errors are big. Energies length becomes M = N + sum(improved)
            energies = np.insert(energies, idx_improved, en_mid[improved])
            # real = np.insert(real, idx_improved, re_mid[improved])
            imag_asf = np.insert(imag_asf, idx_improved, im_mid[improved])

            ### Create new indexes and insert duplicate coefficients for the next iteration
            # Create new array of midpoint locations to evaluate.
            new_value_locs = np.insert(
                arr=np.zeros(
                    real_coefs.shape[0], dtype=bool
                ),  # Copy existing midpoint list
                obj=idx_extra[
                    improved
                ],  # Add the locations where improvements were made
                values=True,  # Insert True
            )
            new_midpoint_locs = np.where(new_value_locs)[
                0
            ]  # Locate the indexes of new midpoint values after insertions.

            # Duplicate coefficients at the improved indexes
            for j in range(idx_improved.shape[0] - 1, -1, -1):
                # Iterate from the last improved coefficient to the start
                idx = (
                    idx_improved[j] - 1
                )  # Get the midpoint index, move back one to get the coefficient index that defines that region
                # Duplicate the coefficients at the index
                real_coefs = np.r_[
                    real_coefs[:idx], [real_coefs[idx]], real_coefs[idx:]
                ]  # Insert the new coefficients

            # Update the set of indexes where new midpoints need to be calculated.
            # Add the midpoint after (in addition to before) the new points.
            # Transpose required otherwise indexes out of order.
            idx_extra = np.vstack(
                (new_midpoint_locs, new_midpoint_locs + 1)
            ).T.flatten()

    # Return values if the maximum number of iterations is reached
    return energies, imag_asf  # , real_coefs
