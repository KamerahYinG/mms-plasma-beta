"""Load MMS1 FPI FAST L2 plasma moments and create a diagnostic overview."""

from datetime import datetime, timezone
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pyspedas


TRANGE = ["2015-10-16/13:03:30", "2015-10-16/13:08:00"]
VARFORMAT = "*numberdensity*|*temppara*|*tempperp*"
OUTPUT = Path(__file__).resolve().parents[1] / "figures" / "mms1_fpi_overview.png"


def utc_string(unix_seconds):
    return datetime.fromtimestamp(float(unix_seconds), tz=timezone.utc).isoformat()


def metadata_units(name):
    metadata = pyspedas.get_data(name, metadata=True)
    data_att = metadata.get("data_att", {})
    if data_att.get("units"):
        return data_att["units"]
    return metadata.get("CDF", {}).get("VATT", {}).get("UNITS", "not reported")


def select_one(variables, species, quantity):
    matches = [
        name
        for name in variables
        if f"_{species}_" in name.lower() and f"_{quantity}_" in name.lower()
        and "_err_" not in name.lower()
        and "_bg_" not in name.lower()
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one returned {species} {quantity} variable; found {matches}."
        )
    return matches[0]


def describe(label, name):
    data = pyspedas.get_data(name)
    values = np.asarray(data.y)
    times = np.asarray(data.times)
    print(f"{label}:")
    print(f"  exact variable name: {name}")
    print(f"  units: {metadata_units(name)}")
    print(f"  array shape: {values.shape}")
    print(f"  number of samples: {len(times)}")
    print(f"  first timestamp: {utc_string(times[0])}")
    print(f"  last timestamp: {utc_string(times[-1])}")
    return times, values


def load_product(datatype, verbose=True):
    variables = pyspedas.projects.mms.fpi(
        trange=TRANGE,
        probe="1",
        data_rate="fast",
        level="l2",
        datatype=datatype,
        varformat=VARFORMAT,
        get_support_data=False,
        time_clip=True,
        spdf=False,
        no_update=True,
    )
    if not variables:
        raise RuntimeError(f"The MMS FPI loader returned no variables for {datatype}.")
    if verbose:
        print(f"Actual variables returned for {datatype.upper()}:")
        for name in variables:
            print(f"  {name}")
    return variables


def scalar_temperature(parallel_times, parallel, perpendicular_times, perpendicular, species):
    if not np.array_equal(parallel_times, perpendicular_times):
        raise RuntimeError(f"{species} temperature moments do not share a native time grid.")
    return parallel_times, (parallel + 2.0 * perpendicular) / 3.0


def finite_series(times, values):
    mask = np.isfinite(times) & np.isfinite(values)
    return times[mask], values[mask]


def plot_series(ne, ni, te, ti):
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True, constrained_layout=True)

    for times, values, label in (*ne, "ne"), (*ni, "ni"):
        axes[0].plot(
            [datetime.fromtimestamp(float(t), tz=timezone.utc) for t in times],
            values,
            label=label,
            linewidth=0.9,
        )
    axes[0].set_ylabel("Density [cm$^{-3}$]")
    axes[0].legend(loc="best")

    axes[1].plot(
        [datetime.fromtimestamp(float(t), tz=timezone.utc) for t in te[0]],
        te[1],
        color="tab:blue",
        linewidth=0.9,
    )
    axes[1].set_ylabel("Electron Te [eV]")

    axes[2].plot(
        [datetime.fromtimestamp(float(t), tz=timezone.utc) for t in ti[0]],
        ti[1],
        color="tab:red",
        linewidth=0.9,
    )
    axes[2].set_ylabel("Ion Ti [eV]")
    axes[2].set_xlabel("UTC time")

    for axis in axes:
        axis.grid(True, alpha=0.25)
    axes[0].set_title("MMS1 FPI FAST L2 plasma moments")
    locator = mdates.AutoDateLocator()
    axes[2].xaxis.set_major_locator(locator)
    axes[2].xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator, tz=timezone.utc))

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=180)
    plt.close(fig)


def main():
    electron_variables = load_product("des-moms")
    ion_variables = load_product("dis-moms")

    selected = {
        "electron number density": select_one(electron_variables, "des", "numberdensity"),
        "electron parallel temperature": select_one(electron_variables, "des", "temppara"),
        "electron perpendicular temperature": select_one(electron_variables, "des", "tempperp"),
        "ion number density": select_one(ion_variables, "dis", "numberdensity"),
        "ion parallel temperature": select_one(ion_variables, "dis", "temppara"),
        "ion perpendicular temperature": select_one(ion_variables, "dis", "tempperp"),
    }

    print("\nSelected moment variables:")
    loaded = {label: describe(label, name) for label, name in selected.items()}

    te = scalar_temperature(
        *loaded["electron parallel temperature"],
        *loaded["electron perpendicular temperature"],
        "electron",
    )
    ti = scalar_temperature(
        *loaded["ion parallel temperature"],
        *loaded["ion perpendicular temperature"],
        "ion",
    )
    ne = finite_series(*loaded["electron number density"])
    ni = finite_series(*loaded["ion number density"])
    te = finite_series(*te)
    ti = finite_series(*ti)

    plot_series(ne, ni, te, ti)

    print("\nFinite-value ranges:")
    for label, (_, values) in (("ne", ne), ("ni", ni), ("Te", te), ("Ti", ti)):
        print(f"  {label}: {np.min(values):.6g} to {np.max(values):.6g}")
    print(f"Saved figure: {OUTPUT}")


if __name__ == "__main__":
    main()
