import numpy as np

from scipy.special import betaln

MAX_CONVERSIONS = 800000

def search_for_barrier(
        z_low: int,
        z_high: int,
        alpha: float,
        power_level: float,
        null_p: float,
        alt_p: float,

):
    log_null_p = np.log(null_p)
    log_null_1_p = np.log(1 - null_p)

    log_alt_p  = np.log(alt_p)
    log_alt_1_p = np.log(1 - alt_p)

    z = z_low + 2 * np.floor((z_high - z_low) / 4)

    while z_low < z_high:
        null_cdf = 0
        alt_cdf = 0

        old_low = z_low
        old_high = z_high

        for n in range(z, MAX_CONVERSIONS + 1, 2):
            k = 0.5 * (n + z)
            prefix = z / n / k
            lbeta_k = betaln(k, n + 1 - k)

            null_cdf += prefix * np.exp(-lbeta_k + (k - z) * log_null_1_p + k * log_null_p)
            alt_cdf += prefix * np.exp(-lbeta_k + (k-z)*log_alt_p + k*log_alt_1_p)

            if (np.isnan(null_cdf) | np.isnan(alt_cdf)):
                break

            if alt_cdf > power_level:
                if null_cdf < alpha:
                    z_high = z
                else:
                    z_low = z + 2

                break
            elif null_cdf > alpha:
                z_low = z + 2

            if (np.isnan(null_cdf)) | (np.isnan(alt_cdf)) | (n >= MAX_CONVERSIONS):
                print("NaN...")
                break

            print(f"High: {z_high}, Low: {z_low}, Z: {Z}, null: {null_cdf}, alt: {alt_cdf}")
            z = z_low + 2 * np.floor((z_high - z_low) / 4)

    return z
