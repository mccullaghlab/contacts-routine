"""Helpers for accepting every trajectory format supported by MDAnalysis."""

from pathlib import Path

import click


def _is_supported_trajectory(path):
    """Return whether MDAnalysis can select a coordinate reader for *path*."""
    from MDAnalysis.coordinates.core import get_reader_for

    try:
        get_reader_for(str(path))
    except (TypeError, ValueError):
        return False
    return True


def trajectory_files(trajectory=None, trajectory_list=None):
    """Return validated trajectory paths from a path, directory, or list file.

    Directory contents are selected through MDAnalysis' reader registry rather
    than from a hard-coded set of suffixes.  Consequently this automatically
    follows the formats supported by the installed MDAnalysis version.
    """
    if trajectory_list:
        source = Path(trajectory_list)
        if source.is_dir():
            candidates = sorted(path for path in source.iterdir() if path.is_file())
            paths = [str(path) for path in candidates if _is_supported_trajectory(path)]
            if not paths:
                raise click.UsageError(
                    f"No MDAnalysis-supported trajectory files found in folder: {source}"
                )
        else:
            with source.open() as stream:
                paths = [line.strip() for line in stream if line.strip()]
    elif trajectory:
        paths = [trajectory]
    else:
        raise click.UsageError(
            "Either --trajectory or --trajectory-list must be provided"
        )

    if not paths:
        raise click.UsageError("No trajectory paths were provided")

    missing = [path for path in paths if not Path(path).is_file()]
    if missing:
        raise click.UsageError(
            "Trajectory file(s) do not exist: " + ", ".join(missing)
        )

    unsupported = [path for path in paths if not _is_supported_trajectory(path)]
    if unsupported:
        raise click.UsageError(
            "Unsupported trajectory format(s) for the installed MDAnalysis "
            "version: " + ", ".join(unsupported)
        )

    return paths
