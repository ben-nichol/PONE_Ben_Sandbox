import numpy as np
from scipy.special import loggamma


def asymmetric_gauss_pdf(x, mu, sigma, r):
    norm = 2.0 / (np.sqrt(2 * np.pi * sigma**2) * (r + 1))
    exp = np.where(
        x < mu,
        np.exp(-0.5 * ((x - mu) / sigma) ** 2),
        np.exp(-0.5 * ((x - mu) / (sigma * r)) ** 2),
    )
    return norm * exp


def gamma_pdf(x, k, theta):
    if x <= 0.0:
        y = 0.0
    else:
        y = (
            1.0
            / (np.exp(loggamma(k)) * theta**k)
            * x ** (k - 1.0)
            * np.exp(-x / theta)
        )
    return y


# --------------------------------------------------


def uniform_weight(d):
    weight = 0.001
    return weight


def gamma_weight(d):
    weight = -4.84e-8 * d * d * d + 6.80e-6 * d * d + 2.06e-3 * d + 3.17e-2
    return weight


def anorm0_weight(d):
    weight = -6.57e-8 * d * d * d + 2.83e-5 * d * d - 6.80e-3 * d + 7.86e-1
    return weight


def anorm1_weight(d):
    weight = 3.65e-8 * d * d * d - 1.84e-5 * d * d + 1.59e-3 * d + 1.12e-1
    return weight


def anorm2_weight(d):
    weight = 8.03e-8 * d * d * d - 1.74e-5 * d * d + 3.17e-3 * d + 6.95e-2
    return weight


def anorm0_mu(d):
    mu_parm = 9.70e-7 * d * d * d - 3.05e-4 * d * d + 1.15e-1 * d - 2.46e0
    return mu_parm


def anorm1_mu(d):
    mu_parm = 8.79e-6 * d * d * d - 1.42e-3 * d * d + 2.29e-1 * d + 2.31e-1
    return mu_parm


def anorm2_mu(d):
    mu_parm = -5.87e-6 * d * d * d + 1.72e-3 * d * d - 1.99e-2 * d + 9.83e-0
    return mu_parm


def anorm0_r(d):
    r_parm = 8.35e-7 * d * d * d - 2.38e-4 * d * d + 2.72e-2 * d + 2.10e-0
    return r_parm


def anorm1_r(d):
    r_parm = 4.28e-6 * d * d * d - 1.04e-3 * d * d + 9.57e-2 * d + 7.81e-0
    return r_parm


def anorm2_r(d):
    r_parm = 0.001 * d + 10.0
    return r_parm


def anorm0_sigma(d):
    sigma_parm = 1.08e-7 * d * d * d - 1.25e-5 * d * d + 7.15e-3 * d + 8.45e-1
    return sigma_parm


def anorm1_sigma(d):
    sigma_parm = 4.15e-7 * d * d * d - 1.01e-4 * d * d + 1.30e-2 * d + 7.14e-1
    return sigma_parm


def anorm2_sigma(d):
    sigma_parm = 1.30e-6 * d * d * d - 2.23e-4 * d * d + 4.19e-2 * d + 2.78e-0
    return sigma_parm


def gamma_scale(d):
    if d < 200.0:
        gamma_parm = 2.854e-3 * d * d - 2.55e-1 * d + 105.0
    else:
        gamma_parm = 300.0
    return gamma_parm


# --------------------------------------------------


def mixture_model_pdf(t_res, dist, tmin=-200.0, tmax=2000.0):
    pdf_value = 0.0
    if (t_res < tmin) or (t_res > tmax) or (dist <= 0) or (dist > 200.0):
        pdf_value += 0.0

    else:
        # uniform
        pdf_value += uniform_weight(dist) * 1.0 / (tmax - tmin)

        # gamma
        pdf_value += gamma_weight(dist) * gamma_pdf(t_res, 2.0, gamma_scale(dist))

        # 3 assymmetric gaussians
        pdf_value += anorm0_weight(dist) * asymmetric_gauss_pdf(
            t_res, anorm0_mu(dist), anorm0_sigma(dist), anorm0_r(dist)
        )

        pdf_value += anorm1_weight(dist) * asymmetric_gauss_pdf(
            t_res, anorm1_mu(dist), anorm1_sigma(dist), anorm1_r(dist)
        )

        pdf_value += anorm2_weight(dist) * asymmetric_gauss_pdf(
            t_res, anorm2_mu(dist), anorm2_sigma(dist), anorm2_r(dist)
        )

    return pdf_value
