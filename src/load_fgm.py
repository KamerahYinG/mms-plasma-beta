"""Download and inspect MMS1 FGM survey L2 magnetic-field data."""

from datetime import datetime, timezone

import numpy as np
import pyspedas


TRANGE = ["2015-10-16/13:03:30", "2015-10-16/13:08:00"]


def utc_string(unix_seconds):
    """Format a Unix timestamp as UTC without changing its precision."""
    return datetime.fromtimestamp(float(unix_seconds), tz=timezone.utc).isoformat()


def metadata_units(metadata):
    """Read units from the metadata structures used by PySPEDAS/PyTplot."""
    data_att = metadata.get("data_att", {})
    if data_att.get("units"):
        return data_att["units"]

    cdf = metadata.get("CDF", {})
    vatt = cdf.get("VATT", {})
    return vatt.get("UNITS", "not reported")


def load_fgm_data(local_only=False, verbose=True):
    """Load the MMS1 FGM interval and return its selected GSE field data."""
    if local_only:
        variables = pyspedas.projects.mms.fgm(
            trange=TRANGE,
            probe="1",
            data_rate="srvy",
            level="l2",
            time_clip=True,
            spdf=False,
            no_update=True,
        )
        source = "project-local cache"
    else:
        variables = pyspedas.projects.mms.fgm(
            trange=TRANGE,
            probe="1",
            data_rate="srvy",
            level="l2",
            time_clip=True,
            spdf=True,
        )

        source = "NASA SPDF"
        if not variables:
            print("NASA SPDF returned no variables; trying the public MMS SDC endpoint.")
            variables = pyspedas.projects.mms.fgm(
                trange=TRANGE,
                probe="1",
                data_rate="srvy",
                level="l2",
                time_clip=True,
                spdf=False,
            )
            source = "public MMS SDC"

    if verbose:
        print(f"Data source used: {source}")
        print("Variables returned by the MMS FGM loader:")
        for name in variables or []:
            print(f"  {name}")

    if not variables:
        raise RuntimeError("The MMS FGM loader returned no variables.")

    candidates = [name for name in variables if "_fgm_b_gse_" in name.lower()]
    selected = None
    selected_data = None
    for name in candidates:
        data = pyspedas.get_data(name)
        if data is not None and np.asarray(data.y).ndim == 2 and data.y.shape[1] >= 4:
            selected = name
            selected_data = data
            break

    if selected is None:
        raise RuntimeError(
            "No returned GSE magnetic-field variable with Bx, By, Bz, and |B| was found."
        )

    values = np.asarray(selected_data.y)
    times = np.asarray(selected_data.times)
    metadata = pyspedas.get_data(selected, metadata=True)
    units = metadata_units(metadata)
    return selected, times, values, units, variables


def main():
    selected, times, values, units, _ = load_fgm_data()
    valid = np.isfinite(times) & np.all(np.isfinite(values[:, :4]), axis=1)
    valid_indices = np.flatnonzero(valid)[:5]

    print("\nSelected GSE magnetic-field variable:")
    print(f"  variable name: {selected}")
    print(f"  array shape: {values.shape}")
    print(f"  first timestamp: {utc_string(times[0])}")
    print(f"  last timestamp: {utc_string(times[-1])}")
    print(f"  units: {units}")
    print(f"  number of samples: {len(times)}")
    print(f"  contains Bx: {values.shape[1] >= 1}")
    print(f"  contains By: {values.shape[1] >= 2}")
    print(f"  contains Bz: {values.shape[1] >= 3}")
    print(f"  contains |B|: {values.shape[1] >= 4}")

    print("\nFirst 5 valid samples (UTC, Bx, By, Bz, |B|):")
    for index in valid_indices:
        bx, by, bz, magnitude = values[index, :4]
        print(
            f"  {utc_string(times[index])}  "
            f"{bx:.6g}  {by:.6g}  {bz:.6g}  {magnitude:.6g}"
        )

    if len(valid_indices) < 5:
        raise RuntimeError(f"Only {len(valid_indices)} valid samples were available.")


if __name__ == "__main__":
    main()
