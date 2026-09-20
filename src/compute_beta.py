"""Calculate MMS1 electron, ion, and total plasma beta on the FPI ion grid."""

from datetime import datetime, timezone
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pyspedas
from scipy.constants import elementary_charge, mu_0

from load_fgm import load_fgm_data
from load_fpi import load_product, select_one


OUTPUT = Path(__file__).resolve().parents[1] / "figures" / "mms1_beta_overview.png"
INTERVALS = {
    "A": ("2015-10-16/13:05:00", "2015-10-16/13:05:35"),
    "B": ("2015-10-16/13:05:50", "2015-10-16/13:06:20"),
}


def thermal_pressure(density_cm3, temperature_ev):
    return density_cm3 * 1.0e6 * temperature_ev * elementary_charge


def magnetic_pressure(magnitude_nt):
    field_t = magnitude_nt * 1.0e-9
    return field_t**2 / (2.0 * mu_0)


def scalar_temperature(parallel, perpendicular):
    return (parallel + 2.0 * perpendicular) / 3.0


def get_series(name):
    data = pyspedas.get_data(name)
    return np.asarray(data.times), np.asarray(data.y, dtype=float)


def quality_mask(times, density, temperature):
    return (
        np.isfinite(times)
        & np.isfinite(density)
        & np.isfinite(temperature)
        & (density > 0.0)
        & (temperature > 0.0)
    )


def interpolate_no_extrapolation(source_times, source_values, target_times):
    result = np.full(target_times.shape, np.nan, dtype=float)
    valid = np.isfinite(source_times) & np.isfinite(source_values)
    x = source_times[valid]
    y = source_values[valid]
    if len(x) < 2:
        return result
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    inside = np.isfinite(target_times) & (target_times >= x[0]) & (target_times <= x[-1])
    result[inside] = np.interp(target_times[inside], x, y)
    return result


def utc_datetimes(times):
    return [datetime.fromtimestamp(float(t), tz=timezone.utc) for t in times]


def utc_string(timestamp):
    return datetime.fromtimestamp(float(timestamp), tz=timezone.utc).isoformat()


def validate_physics():
    density = 5.0
    temperature = 100.0
    field = 10.0
    pressure = thermal_pressure(density, temperature)
    p_b = magnetic_pressure(field)
    beta = pressure / p_b

    checks = {
        "doubling density doubles thermal pressure": np.isclose(
            thermal_pressure(2.0 * density, temperature), 2.0 * pressure
        ),
        "doubling temperature doubles thermal pressure": np.isclose(
            thermal_pressure(density, 2.0 * temperature), 2.0 * pressure
        ),
        "doubling B quadruples magnetic pressure": np.isclose(
            magnetic_pressure(2.0 * field), 4.0 * p_b
        ),
        "doubling B reduces beta by four": np.isclose(
            pressure / magnetic_pressure(2.0 * field), beta / 4.0
        ),
    }
    print("Numerical validation checks:")
    for description, passed in checks.items():
        print(f"  {'PASS' if passed else 'FAIL'}: {description}")
    if not all(checks.values()):
        raise RuntimeError("One or more plasma-beta validation checks failed.")


def interval_medians(label, bounds, common):
    start, end = (pyspedas.time_double(value) for value in bounds)
    times = common["times"]
    interval = (times >= start) & (times <= end)
    print(f"Interval {label} ({bounds[0]} to {bounds[1]} UTC):")
    fields = (
        ("|B|", "b", "nT"),
        ("ne", "ne", "cm^-3"),
        ("ni", "ni", "cm^-3"),
        ("Te", "te", "eV"),
        ("Ti", "ti", "eV"),
        ("p_thermal", "p_thermal_npa", "nPa"),
        ("p_B", "p_b_npa", "nPa"),
        ("beta_total", "beta_total", ""),
    )
    for display, key, unit in fields:
        values = common[key][interval]
        finite = values[np.isfinite(values)]
        median = np.median(finite) if len(finite) else np.nan
        print(f"  median {display}: {median:.6g}{' ' + unit if unit else ''}")


def main():
    validate_physics()

    _, fgm_times, fgm_field, _, _ = load_fgm_data(local_only=True, verbose=False)
    electron_variables = load_product("des-moms", verbose=False)
    ion_variables = load_product("dis-moms", verbose=False)

    ne_times, ne = get_series(select_one(electron_variables, "des", "numberdensity"))
    ne_para_times, te_para = get_series(select_one(electron_variables, "des", "temppara"))
    ne_perp_times, te_perp = get_series(select_one(electron_variables, "des", "tempperp"))
    ni_times, ni = get_series(select_one(ion_variables, "dis", "numberdensity"))
    ni_para_times, ti_para = get_series(select_one(ion_variables, "dis", "temppara"))
    ni_perp_times, ti_perp = get_series(select_one(ion_variables, "dis", "tempperp"))

    if not (np.array_equal(ne_times, ne_para_times) and np.array_equal(ne_times, ne_perp_times)):
        raise RuntimeError("Electron density and temperature moments use different native grids.")
    if not (np.array_equal(ni_times, ni_para_times) and np.array_equal(ni_times, ni_perp_times)):
        raise RuntimeError("Ion density and temperature moments use different native grids.")

    te = scalar_temperature(te_para, te_perp)
    ti = scalar_temperature(ti_para, ti_perp)
    electron_valid = quality_mask(ne_times, ne, te)
    ion_valid = quality_mask(ni_times, ni, ti)

    pe = np.full(ne.shape, np.nan)
    pi = np.full(ni.shape, np.nan)
    pe[electron_valid] = thermal_pressure(ne[electron_valid], te[electron_valid])
    pi[ion_valid] = thermal_pressure(ni[ion_valid], ti[ion_valid])

    common_times = ni_times.copy()
    pe_common = interpolate_no_extrapolation(ne_times, pe, common_times)
    ne_common = interpolate_no_extrapolation(ne_times, np.where(electron_valid, ne, np.nan), common_times)
    te_common = interpolate_no_extrapolation(ne_times, np.where(electron_valid, te, np.nan), common_times)

    b_magnitude = fgm_field[:, 3].astype(float)
    b_valid = np.isfinite(fgm_times) & np.isfinite(b_magnitude) & (b_magnitude > 0.0)
    b_common = interpolate_no_extrapolation(
        fgm_times[b_valid], b_magnitude[b_valid], common_times
    )
    p_b = np.full(common_times.shape, np.nan)
    valid_b_common = np.isfinite(b_common) & (b_common > 0.0)
    p_b[valid_b_common] = magnetic_pressure(b_common[valid_b_common])

    beta_e = pe_common / p_b
    beta_i = pi / p_b
    p_thermal = pi + pe_common
    beta_total = p_thermal / p_b

    common = {
        "times": common_times,
        "b": b_common,
        "ne": ne_common,
        "ni": np.where(ion_valid, ni, np.nan),
        "te": te_common,
        "ti": np.where(ion_valid, ti, np.nan),
        "p_thermal_npa": p_thermal * 1.0e9,
        "p_b_npa": p_b * 1.0e9,
        "beta_total": beta_total,
    }

    fig, axes = plt.subplots(5, 1, figsize=(11, 13), sharex=True, constrained_layout=True)
    fgm_dt = utc_datetimes(fgm_times)
    for column, label in enumerate(("Bx", "By", "Bz", "|B|")):
        axes[0].plot(fgm_dt, fgm_field[:, column], label=label, linewidth=0.75)
    axes[0].set_ylabel("B [nT]")
    axes[0].set_title("MMS1 FGM and FPI FAST L2 overview")
    axes[0].legend(loc="best", ncol=4)

    axes[1].plot(utc_datetimes(ne_times), ne, label="ne", linewidth=0.9)
    axes[1].plot(utc_datetimes(ni_times), ni, label="ni", linewidth=0.9)
    axes[1].set_ylabel("Density [cm$^{-3}$]")
    axes[1].legend(loc="best")

    axes[2].plot(utc_datetimes(ne_times), te, label="Te", linewidth=0.9)
    axes[2].plot(utc_datetimes(ni_times), ti, label="Ti", linewidth=0.9)
    axes[2].set_ylabel("Temperature [eV]")
    axes[2].legend(loc="best")

    common_dt = utc_datetimes(common_times)
    axes[3].plot(common_dt, p_thermal * 1.0e9, label="pe + pi", linewidth=0.9)
    axes[3].plot(common_dt, p_b * 1.0e9, label="pB", linewidth=0.9)
    axes[3].set_ylabel("Pressure [nPa]")
    axes[3].legend(loc="best")

    axes[4].plot(common_dt, beta_e, label="beta_e", linewidth=0.9)
    axes[4].plot(common_dt, beta_i, label="beta_i", linewidth=0.9)
    axes[4].plot(common_dt, beta_total, label="beta_total", linewidth=1.1)
    positive_beta = np.concatenate((beta_e, beta_i, beta_total))
    positive_beta = positive_beta[np.isfinite(positive_beta)]
    if len(positive_beta) and np.all(positive_beta > 0.0):
        axes[4].set_yscale("log")
    axes[4].axhline(1.0, color="gray", linestyle="--", linewidth=0.8, label="beta = 1")
    axes[4].set_ylabel("Plasma beta")
    axes[4].set_xlabel("UTC time")
    axes[4].legend(loc="best", ncol=2)

    for axis in axes:
        axis.grid(True, alpha=0.25)
    locator = mdates.AutoDateLocator()
    axes[4].xaxis.set_major_locator(locator)
    axes[4].xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator, tz=timezone.utc))

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=180)
    plt.close(fig)

    print("\nInterval medians:")
    for label, bounds in INTERVALS.items():
        interval_medians(label, bounds, common)

    transition = pyspedas.time_double("2015-10-16/13:05:45")
    index = int(np.searchsorted(common_times, transition))
    print("\nNative ion-grid beta_total samples surrounding 13:05:45 UTC:")
    for sample in (index - 2, index - 1, index):
        print(f"  {utc_string(common_times[sample])}: {beta_total[sample]:.6g}")
    print(f"Saved figure: {OUTPUT}")


if __name__ == "__main__":
    main()
