"""Plot an overview of the cached MMS1 FGM survey magnetic field."""

from datetime import datetime, timezone
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from load_fgm import load_fgm_data


OUTPUT = Path(__file__).resolve().parents[1] / "figures" / "mms1_fgm_overview.png"


def utc_string(unix_seconds):
    return datetime.fromtimestamp(float(unix_seconds), tz=timezone.utc).isoformat()


def main():
    variable, times, field, units, _ = load_fgm_data(local_only=True, verbose=False)
    if variable != "mms1_fgm_b_gse_srvy_l2":
        raise RuntimeError(f"Unexpected magnetic-field variable selected: {variable}")

    valid = np.isfinite(times) & np.all(np.isfinite(field[:, :4]), axis=1)
    times = times[valid]
    field = field[valid, :4]
    datetimes = [datetime.fromtimestamp(float(t), tz=timezone.utc) for t in times]

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, constrained_layout=True)
    axes[0].plot(datetimes, field[:, 0], label="Bx", linewidth=0.8)
    axes[0].plot(datetimes, field[:, 1], label="By", linewidth=0.8)
    axes[0].plot(datetimes, field[:, 2], label="Bz", linewidth=0.8)
    axes[0].set_ylabel("Magnetic field [nT]")
    axes[0].set_title("MMS1 FGM magnetic field (GSE)")
    axes[0].legend(loc="best")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(datetimes, field[:, 3], color="black", linewidth=0.8)
    axes[1].set_ylabel("Magnetic field [nT]")
    axes[1].set_xlabel("UTC time")
    axes[1].grid(True, alpha=0.25)

    locator = mdates.AutoDateLocator()
    axes[1].xaxis.set_major_locator(locator)
    axes[1].xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator, tz=timezone.utc))

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=180)
    plt.close(fig)

    labels = ("Bx", "By", "Bz", "|B|")
    print(f"Loaded variable: {variable}")
    print(f"Units: {units}")
    for column, label in enumerate(labels):
        print(f"{label} minimum: {np.min(field[:, column]):.6g} nT")
        print(f"{label} maximum: {np.max(field[:, column]):.6g} nT")

    largest_index = int(np.argmax(field[:, 3]))
    print(f"Time of largest |B|: {utc_string(times[largest_index])}")
    print(f"Saved figure: {OUTPUT}")


if __name__ == "__main__":
    main()
