#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

CATALOG = "SPARC_Lelli2016c.mrt.txt"
MASS_MODELS = "MassModels_Lelli2016c.mrt.txt"

G = 4.30091e-6  # kpc (km/s)^2 Msun^-1
LOWER_BOUNDS = np.array([4.0, -2.0])
UPPER_BOUNDS = np.array([12.0, 3.0])
STARTS = [(a, b) for a in [5.5, 7.5, 9.5, 11.0]
                 for b in [-1.0, 0.0, 1.0, 2.0]]

BOUNDARY_RMS = 33.68304987427586
BOUNDARY_K = 3


def parse_catalog(path: Path) -> pd.DataFrame:
    rows = []
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for index, line in enumerate(handle):
            if index < 98:
                continue
            parts = line.split()
            if len(parts) < 18:
                continue
            try:
                rows.append(
                    {
                        "Galaxy": parts[0],
                        "L36": float(parts[7]),
                        "MHI": float(parts[13]),
                        "RHI": float(parts[14]),
                    }
                )
            except (ValueError, IndexError):
                pass
    return pd.DataFrame(rows).drop_duplicates("Galaxy")


def parse_mass_models(path: Path) -> pd.DataFrame:
    rows = []
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for index, line in enumerate(handle):
            if index < 25:
                continue
            parts = line.split()
            if len(parts) < 10:
                continue
            try:
                rows.append(
                    {
                        "Galaxy": parts[0],
                        "R": float(parts[2]),
                        "Vobs": float(parts[3]),
                        "Vgas": float(parts[5]),
                        "Vdisk": float(parts[6]),
                        "Vbul": float(parts[7]),
                    }
                )
            except (ValueError, IndexError):
                pass
    return pd.DataFrame(rows)


def nfw_velocity_squared(
    radius_kpc: np.ndarray,
    log10_rho_s: float,
    log10_r_s: float,
) -> np.ndarray:
    rho_s = 10.0 ** log10_rho_s
    r_s = 10.0 ** log10_r_s
    x = np.maximum(radius_kpc / r_s, 1.0e-12)
    enclosed_mass = (
        4.0
        * np.pi
        * rho_s
        * r_s**3
        * (np.log1p(x) - x / (1.0 + x))
    )
    return G * enclosed_mass / np.maximum(radius_kpc, 1.0e-12)


def main() -> None:
    root = Path(__file__).resolve().parent
    outputs = root / "outputs"
    outputs.mkdir(exist_ok=True)

    catalog = parse_catalog(root / CATALOG)
    mass_models = parse_mass_models(root / MASS_MODELS)

    eligible = catalog.loc[catalog["RHI"] > 0].copy()
    data = mass_models.merge(eligible, on="Galaxy", how="inner")

    data["Vbar2"] = (
        data["Vdisk"] ** 2
        + data["Vbul"] ** 2
        + data["Vgas"] * np.abs(data["Vgas"])
    )
    excluded = data.loc[data["Vbar2"] < 0, ["Galaxy", "R", "Vbar2"]].copy()
    data = data.loc[data["Vbar2"] >= 0].copy()

    fit_rows = []
    prediction_rows = []

    for galaxy, group in data.groupby("Galaxy", sort=True):
        radius = group["R"].to_numpy(float)
        observed = group["Vobs"].to_numpy(float)
        vbar2 = group["Vbar2"].to_numpy(float)

        def predicted(theta: np.ndarray) -> np.ndarray:
            halo_v2 = nfw_velocity_squared(radius, theta[0], theta[1])
            return np.sqrt(np.maximum(vbar2 + halo_v2, 0.0))

        best_rss = np.inf
        best_theta = None

        for start in STARTS:
            result = least_squares(
                lambda theta: predicted(theta) - observed,
                x0=np.asarray(start, dtype=float),
                bounds=(LOWER_BOUNDS, UPPER_BOUNDS),
                method="trf",
                xtol=1.0e-10,
                ftol=1.0e-10,
                gtol=1.0e-10,
                max_nfev=2000,
            )
            rss = float(np.sum(result.fun**2))
            if rss < best_rss:
                best_rss = rss
                best_theta = result.x.copy()

        prediction = predicted(best_theta)
        halo_v2 = nfw_velocity_squared(radius, best_theta[0], best_theta[1])

        fit_rows.append(
            {
                "Galaxy": galaxy,
                "n_points": len(group),
                "log10_rho_s_Msun_kpc3": float(best_theta[0]),
                "log10_r_s_kpc": float(best_theta[1]),
                "rho_s_Msun_kpc3": float(10.0 ** best_theta[0]),
                "r_s_kpc": float(10.0 ** best_theta[1]),
                "rss": best_rss,
                "rms_km_s": float(np.sqrt(best_rss / len(group))),
            }
        )

        out = group.copy()
        out["Vhalo_NFW"] = np.sqrt(np.maximum(halo_v2, 0.0))
        out["Vpred_NFW"] = prediction
        out["Residual_NFW"] = out["Vobs"] - out["Vpred_NFW"]
        prediction_rows.append(out)

    fits = pd.DataFrame(fit_rows)
    predictions = pd.concat(prediction_rows, ignore_index=True)

    n = len(predictions)
    galaxy_count = predictions["Galaxy"].nunique()
    rss_nfw = float(np.sum(predictions["Residual_NFW"] ** 2))
    rms_nfw = float(np.sqrt(rss_nfw / n))
    tss = float(np.sum((predictions["Vobs"] - predictions["Vobs"].mean()) ** 2))
    r2_nfw = float(1.0 - rss_nfw / tss)

    k_nfw = 2 * galaxy_count
    aic_nfw = float(n * np.log(rss_nfw / n) + 2 * k_nfw)
    bic_nfw = float(n * np.log(rss_nfw / n) + k_nfw * np.log(n))

    rss_boundary = float(n * BOUNDARY_RMS**2)
    aic_boundary = float(
        n * np.log(rss_boundary / n) + 2 * BOUNDARY_K
    )
    bic_boundary = float(
        n * np.log(rss_boundary / n) + BOUNDARY_K * np.log(n)
    )

    summary = {
        "status": "LOCKED_REPRODUCIBLE_BENCHMARK",
        "sample": {
            "galaxies": int(galaxy_count),
            "points": int(n),
            "excluded_points": excluded.to_dict(orient="records"),
        },
        "nfw_definition": {
            "parameters": [
                "log10(rho_s / Msun kpc^-3)",
                "log10(r_s / kpc)",
            ],
            "bounds": {
                "log10_rho_s": [4.0, 12.0],
                "log10_r_s": [-2.0, 3.0],
            },
            "objective": "unweighted velocity residual sum of squares per galaxy",
            "optimizer": "scipy.optimize.least_squares",
            "method": "trf",
            "starts": STARTS,
            "xtol": 1.0e-10,
            "ftol": 1.0e-10,
            "gtol": 1.0e-10,
            "max_nfev": 2000,
        },
        "nfw_result": {
            "rms_km_s": rms_nfw,
            "r2": r2_nfw,
            "rss": rss_nfw,
            "k": int(k_nfw),
            "aic": aic_nfw,
            "bic": bic_nfw,
        },
        "boundary_result": {
            "rms_km_s": BOUNDARY_RMS,
            "k": BOUNDARY_K,
            "aic": aic_boundary,
            "bic": bic_boundary,
        },
        "comparison": {
            "delta_aic_nfw_minus_boundary": aic_nfw - aic_boundary,
            "delta_bic_nfw_minus_boundary": bic_nfw - bic_boundary,
        },
    }

    fits.to_csv(outputs / "nfw_per_galaxy_best_fits.csv", index=False)
    predictions.to_csv(outputs / "nfw_per_point_predictions.csv", index=False)
    excluded.to_csv(outputs / "excluded_points.csv", index=False)
    (outputs / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
