#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter
import numpy as np
import pandas as pd


PARENT_FOLDER = "DOMP1ComparisionRectangularEnclosureHotWallColdMedium"
CASE_TITLE = "Hot-wall cold-medium"
ABSORPTION_COEFFICIENT = 20.0
SCATTERING_COEFFICIENT = 0.0
DOMAIN_LENGTH = 1.0
SIGMA = 5.670374419e-8
MEDIUM_TEMPERATURE = 0.0
WALL_TEMPERATURE = 500.0

P1_PROFILE_TEST = "P1RadiativeEnergyAndFluxProfileTest"
DOM_PROFILE_TEST = "DOMRadiativeEnergyAndFluxProfileTest"
P1_MESH_TEST = "P1MeshConvergenceTest"
DOM_MESH_TEST = "DOMMeshConvergenceTest"
DOM_ANGULAR_TEST = "DOMAngularConvergenceTest"

SUPPORTED_ORDINATES = (2, 4, 6, 8)
SUPPORTED_MESHES = (20, 40, 80, 160, 320)


@dataclass
class ProfileData:
    label: str
    directory: Path
    x: np.ndarray
    g: np.ndarray
    qx: np.ndarray
    qmag: np.ndarray
    source: np.ndarray


@dataclass
class RectilinearField:
    x: np.ndarray
    y: np.ndarray
    values: np.ndarray


def find_project_root() -> Path:
    script_directory = Path(__file__).resolve().parent
    for candidate in (script_directory, *script_directory.parents):
        if (candidate / "PathConfig.hpp").exists() and (candidate / "ValidationTests").exists():
            return candidate
    return script_directory.parent.parent


def apply_report_style() -> None:
    plt.rcParams.update({
        "figure.figsize": (9.5, 6.2),
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "font.size": 16,
        "axes.titlesize": 17,
        "axes.labelsize": 18,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 11.5,
        "lines.linewidth": 2.4,
        "lines.markersize": 7.0,
        "axes.grid": True,
        "grid.alpha": 0.28,
        "grid.linestyle": "--",
        "axes.linewidth": 1.15,
        "legend.frameon": True,
        "legend.framealpha": 0.95,
        "legend.edgecolor": "0.75",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })



PROFILE_STYLES = {
    "P1": {"color": "#D81B60", "linestyle": "-", "linewidth": 2.8},
    "DOM S2": {"color": "#FF7F00", "linestyle": "--", "linewidth": 2.6},
    "DOM S4": {"color": "#00A087", "linestyle": "-.", "linewidth": 2.6},
    "DOM S6": {"color": "#0072B2", "linestyle": ":", "linewidth": 2.9},
    "DOM S8": {"color": "#7B2CBF", "linestyle": "-", "linewidth": 2.6},
    "Paper benchmark": {
        "color": "#111111", "linestyle": "None", "marker": "o",
        "markersize": 5.5, "markerfacecolor": "none", "markeredgewidth": 1.2, "markevery": 12,
    },
}


def profile_style(label: str, index: int) -> dict:
    normalized = label.replace("$", "").replace("_", "").strip()
    normalized = normalized.replace("DOM S", "DOM S")
    if "Paper" in label or "benchmark" in label.lower() or "Sarr and Groth" in label:
        return PROFILE_STYLES["Paper benchmark"].copy()
    for key in ("P1", "DOM S2", "DOM S4", "DOM S6", "DOM S8"):
        if key.replace(" ", "") in normalized.replace(" ", ""):
            return PROFILE_STYLES[key].copy()
    fallback_colors = ("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#00BFC4")
    return {"color": fallback_colors[index % len(fallback_colors)], "linestyle": "-", "linewidth": 2.4}


def place_legend_outside(ax: plt.Axes, columns: int = 1) -> None:
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        ncol=columns,
        frameon=True,
        handlelength=3.0,
        labelspacing=0.65,
    )

def finish_axes(ax: plt.Axes, title: str, xlabel: str, ylabel: str, log_x: bool = False, log_y: bool = False) -> None:
    ax.set_title(title, pad=14, fontweight="bold")
    ax.set_xlabel(xlabel, labelpad=9)
    ax.set_ylabel(ylabel, labelpad=9)
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax.tick_params(axis="both", which="major", length=7, width=1.1)
    ax.tick_params(axis="both", which="minor", length=4, width=0.85)
    if not log_x:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    if not log_y:
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    for spine in ax.spines.values():
        spine.set_linewidth(1.15)


def save_figure(fig: plt.Figure, directory: Path, stem: str, created_files: list[Path]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        output_path = directory / f"{stem}.{extension}"
        fig.savefig(output_path, bbox_inches="tight")
        created_files.append(output_path)
    plt.close(fig)


def warn_missing(path: Path) -> None:
    warnings.warn(f"Skipping missing result: {path}")


def clean_columns(dataframe: pd.DataFrame) -> pd.DataFrame:
    dataframe = dataframe.copy()
    dataframe.columns = [str(column).strip() for column in dataframe.columns]
    return dataframe


def read_standard_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    dataframe = clean_columns(pd.read_csv(path))
    for column in dataframe.columns:
        converted = pd.to_numeric(dataframe[column], errors="coerce")
        if converted.notna().sum() == dataframe[column].notna().sum():
            dataframe[column] = converted
    return dataframe


def read_metadata_csv(path: Path, first_header_column: str = "x") -> tuple[dict[str, str], pd.DataFrame]:
    if not path.exists():
        raise FileNotFoundError(path)

    lines = path.read_text().splitlines()
    header_index = None
    metadata: dict[str, str] = {}

    for index, line in enumerate(lines):
        stripped = line.strip()
        if not stripped:
            continue
        tokens = [token.strip() for token in stripped.split(",")]
        if tokens and tokens[0] == first_header_column:
            header_index = index
            break
        if len(tokens) >= 2:
            metadata[tokens[0]] = tokens[1]

    if header_index is None:
        raise ValueError(f"Could not find the '{first_header_column}' data header in {path}")

    dataframe = clean_columns(pd.read_csv(path, skiprows=header_index))
    dataframe = dataframe.dropna(how="all")
    for column in dataframe.columns:
        dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce")
    dataframe = dataframe.dropna(how="all")
    return metadata, dataframe


def first_existing_column(dataframe: pd.DataFrame, candidates: Sequence[str]) -> str:
    for candidate in candidates:
        if candidate in dataframe.columns:
            return candidate
    raise KeyError(f"None of {candidates} were found. Available columns: {list(dataframe.columns)}")


def ordinate_number(path: Path) -> int:
    match = re.search(r"ordinate(\d+)", str(path))
    if match is None:
        raise ValueError(f"Could not determine ordinate order from {path}")
    return int(match.group(1))


def mesh_number(path: Path) -> int:
    match = re.match(r"(\d+)x(\d+)x(\d+)", path.name)
    if match is None:
        raise ValueError(f"Could not determine mesh size from {path}")
    return int(match.group(1))

def load_profile_data(directory: Path, label: str) -> ProfileData:
    _, g_data = read_metadata_csv(directory / "centerlineG.csv")
    _, flux_data = read_metadata_csv(directory / "centerlineFlux.csv")
    _, source_data = read_metadata_csv(directory / "centerlineRadiativeSourceTerm.csv")

    x_column = first_existing_column(g_data, ("x_over_L", "x"))
    g_column = first_existing_column(g_data, ("G_P1", "G_DOM", "G"))
    flux_x_column = first_existing_column(flux_data, ("radiative_flux_x_P1", "radiative_flux_x_DOM"))
    flux_magnitude_column = first_existing_column(flux_data, ("radiative_flux_magnitude_P1", "radiative_flux_magnitude_DOM"))
    source_column = first_existing_column(source_data, ("radiative_source_term_P1", "radiative_source_term_DOM"))

    x = g_data[x_column].to_numpy(dtype=float)
    if x_column == "x":
        x = x / DOMAIN_LENGTH

    return ProfileData(
        label=label,
        directory=directory,
        x=x,
        g=g_data[g_column].to_numpy(dtype=float),
        qx=flux_data[flux_x_column].to_numpy(dtype=float),
        qmag=flux_data[flux_magnitude_column].to_numpy(dtype=float),
        source=source_data[source_column].to_numpy(dtype=float),
    )



def load_paper_benchmark(case_root: Path) -> ProfileData | None:

    script_directory = Path(__file__).resolve().parent
    candidates = (
        script_directory / "paper_benchmark.csv",
        script_directory / "sarr_groth_benchmark.csv",
        case_root / "paper_benchmark.csv",
        case_root / "sarr_groth_benchmark.csv",
    )
    path = next((candidate for candidate in candidates if candidate.exists()), None)
    if path is None:
        warnings.warn(
            "Paper benchmark CSV was not found. The three comparison plots will "
            "contain P1 and DOM S2/S4/S6/S8 only. Add paper_benchmark.csv with "
            "columns x_over_L,G,qx,source to include the published curve."
        )
        return None

    data = read_standard_csv(path)
    x_column = first_existing_column(data, ("x_over_L", "x"))
    g_column = first_existing_column(data, ("G", "radiative_energy", "irradiation"))
    qx_column = first_existing_column(data, ("qx", "radiative_flux_x", "radiative_heat_flux"))
    source_column = first_existing_column(data, ("source", "Sr", "radiative_source_term"))

    x = data[x_column].to_numpy(dtype=float)
    if x_column == "x":
        x = x / DOMAIN_LENGTH
    qx = data[qx_column].to_numpy(dtype=float)
    return ProfileData(
        label="Sarr and Groth DOM",
        directory=path.parent,
        x=x,
        g=data[g_column].to_numpy(dtype=float),
        qx=qx,
        qmag=np.abs(qx),
        source=data[source_column].to_numpy(dtype=float),
    )


def discover_profile_data(case_root: Path) -> list[ProfileData]:
    profiles: list[ProfileData] = []

    p1_directory = case_root / P1_PROFILE_TEST / "320x320x1"
    if (p1_directory / "centerlineG.csv").exists():
        profiles.append(load_profile_data(p1_directory, "P1"))
    else:
        warn_missing(p1_directory / "centerlineG.csv")

    dom_root = case_root / DOM_PROFILE_TEST
    if dom_root.exists():
        ordinate_directories = sorted(dom_root.glob("ordinate*/320x320x1"), key=ordinate_number)
        for directory in ordinate_directories:
            if all((directory / filename).exists() for filename in ("centerlineG.csv", "centerlineFlux.csv", "centerlineRadiativeSourceTerm.csv")):
                try:
                    profiles.append(load_profile_data(directory, f"DOM S{ordinate_number(directory)}"))
                except (ValueError, KeyError, pd.errors.EmptyDataError) as error:
                    warnings.warn(f"Skipping incomplete profile data in {directory}: {error}")

    benchmark = load_paper_benchmark(case_root)
    if benchmark is not None:
        profiles.append(benchmark)

    return profiles


def plot_profile_quantity(profiles: Sequence[ProfileData], attribute: str, title: str, ylabel: str, output_directory: Path, stem: str, created_files: list[Path], log_y: bool = False) -> None:
    if not profiles:
        return

    fig, ax = plt.subplots(figsize=(10.8, 6.4))

    for index, profile in enumerate(profiles):
        values = np.asarray(getattr(profile, attribute), dtype=float)
        if log_y:
            values = np.maximum(np.abs(values), 1.0e-30)
        style = profile_style(profile.label, index)
        ax.plot(profile.x, values, label=profile.label, **style)

    finish_axes(ax, title, "Normalised coordinate, $x/L$ [-]", ylabel, log_y=log_y)
    ax.set_xlim(0.0, 1.0)
    if len(profiles) > 1:
        place_legend_outside(ax)
        fig.subplots_adjust(right=0.76, top=0.88, bottom=0.15, left=0.13)
    else:
        ax.legend(loc="best")
        fig.subplots_adjust(top=0.88, bottom=0.15, left=0.13, right=0.96)
    save_figure(fig, output_directory, stem, created_files)


def plot_individual_profiles(profile: ProfileData, created_files: list[Path]) -> None:
    output_directory = profile.directory / "Plots"
    model_name = profile.label

    quantities = (
        ("g", f"{CASE_TITLE}: {model_name} centreline radiative energy", "Radiative energy, $G$ [W m$^{-2}$]", "centerline_radiative_energy"),
        ("qx", f"{CASE_TITLE}: {model_name} signed centreline radiative flux", "Radiative flux, $q_x$ [W m$^{-2}$]", "centerline_signed_radiative_flux"),
        ("qmag", f"{CASE_TITLE}: {model_name} centreline flux magnitude", "Radiative flux magnitude, $|q_x|$ [W m$^{-2}$]", "centerline_radiative_flux_magnitude"),
        ("source", f"{CASE_TITLE}: {model_name} centreline radiative source", "Radiative source, $S_r$ [W m$^{-3}$]", "centerline_radiative_source"),
    )

    for attribute, title, ylabel, stem in quantities:
        plot_profile_quantity((profile,), attribute, title, ylabel, output_directory, stem, created_files)


def plot_symmetry_file(path: Path, label: str, created_files: list[Path]) -> None:
    if not path.exists():
        return

    dataframe = read_standard_csv(path)
    if "maximum_relative_error" not in dataframe.columns:
        return

    quantities = dataframe["quantity"].astype(str).tolist()
    values = np.maximum(pd.to_numeric(dataframe["maximum_relative_error"], errors="coerce").to_numpy(dtype=float), 1.0e-18)

    fig, ax = plt.subplots()
    ax.bar(quantities, values)
    finish_axes(ax, f"{CASE_TITLE}: {label} centreline symmetry checks", "Quantity", "Maximum relative symmetry error [-]", log_y=True)
    ax.tick_params(axis="x", rotation=18)
    save_figure(fig, path.parent / "Plots", "centerline_symmetry_checks", created_files)


def plot_all_profiles(case_root: Path, created_files: list[Path]) -> list[ProfileData]:
    profiles = discover_profile_data(case_root)
    comparison_directory = case_root / "Plots"

    plot_profile_quantity(
        profiles,
        "g",
        f"{CASE_TITLE}: centreline irradiation",
        "Irradiation, $G$ [W m$^{-2}$]",
        comparison_directory,
        "centreline_irradiation_all_models",
        created_files,
    )
    plot_profile_quantity(
        profiles,
        "qx",
        f"{CASE_TITLE}: centreline radiative heat flux",
        "Radiative heat flux, $q_x$ [W m$^{-2}$]",
        comparison_directory,
        "centreline_radiative_heat_flux_all_models",
        created_files,
    )
    plot_profile_quantity(
        profiles,
        "source",
        f"{CASE_TITLE}: centreline radiative source term",
        "Radiative source term, $S_r$ [W m$^{-3}$]",
        comparison_directory,
        "centreline_radiative_source_term_all_models",
        created_files,
    )
    return profiles


def plot_mesh_convergence_csv(path: Path, label: str, created_files: list[Path]) -> pd.DataFrame | None:
    if not path.exists():
        warn_missing(path)
        return None

    dataframe = read_standard_csv(path)
    output_directory = path.parent / "Plots"
    x = pd.to_numeric(dataframe["cells_x"], errors="coerce").to_numpy(dtype=float)
    total_cells = pd.to_numeric(dataframe["total_cells_xy"], errors="coerce").to_numpy(dtype=float)
    relative_change = pd.to_numeric(dataframe["successive_grid_relative_L2_change"], errors="coerce").to_numpy(dtype=float)
    observed_order = pd.to_numeric(dataframe["observed_order"], errors="coerce").to_numpy(dtype=float)
    runtime = pd.to_numeric(dataframe["runtime"], errors="coerce").to_numpy(dtype=float)

    mask = np.isfinite(relative_change) & (relative_change > 0.0)
    if np.any(mask):
        fig, ax = plt.subplots()
        ax.plot(x[mask], relative_change[mask], marker="o")
        finish_axes(ax, f"{CASE_TITLE}: {label} spatial convergence", "Cells in each direction, $N$ [-]", "Successive-grid relative $L_2$ change [-]", log_x=True, log_y=True)
        save_figure(fig, output_directory, "mesh_convergence_relative_L2_change", created_files)

    mask = np.isfinite(observed_order)
    if np.any(mask):
        fig, ax = plt.subplots()
        ax.plot(x[mask], observed_order[mask], marker="o")
        finish_axes(ax, f"{CASE_TITLE}: {label} observed spatial order", "Cells in each direction, $N$ [-]", "Observed order [-]", log_x=True)
        save_figure(fig, output_directory, "mesh_convergence_observed_order", created_files)

    mask = np.isfinite(runtime) & (runtime > 0.0)
    if np.any(mask):
        fig, ax = plt.subplots()
        ax.plot(total_cells[mask], runtime[mask], marker="o")
        finish_axes(ax, f"{CASE_TITLE}: {label} runtime scaling", "Number of cells, $N_xN_y$ [-]", "Runtime [s]", log_x=True, log_y=True)
        save_figure(fig, output_directory, "mesh_convergence_runtime", created_files)

    return dataframe


def plot_combined_mesh_convergence(case_root: Path, datasets: Sequence[tuple[str, pd.DataFrame]], created_files: list[Path]) -> None:
    if not datasets:
        return

    output_directory = case_root / DOM_MESH_TEST / "Plots"

    fig, ax = plt.subplots()
    plotted = False
    for label, dataframe in datasets:
        x = pd.to_numeric(dataframe["cells_x"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(dataframe["successive_grid_relative_L2_change"], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(y) & (y > 0.0)
        if np.any(mask):
            ax.plot(x[mask], y[mask], marker="o", label=label)
            plotted = True
    if plotted:
        finish_axes(ax, f"{CASE_TITLE}: spatial-convergence comparison", "Cells in each direction, $N$ [-]", "Successive-grid relative $L_2$ change [-]", log_x=True, log_y=True)
        place_legend_outside(ax)
        fig.subplots_adjust(right=0.76)
        save_figure(fig, output_directory, "comparison_mesh_convergence_relative_L2_change", created_files)
    else:
        plt.close(fig)

    fig, ax = plt.subplots()
    plotted = False
    for label, dataframe in datasets:
        x = pd.to_numeric(dataframe["total_cells_xy"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(dataframe["runtime"], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(y) & (y > 0.0)
        if np.any(mask):
            ax.plot(x[mask], y[mask], marker="o", label=label)
            plotted = True
    if plotted:
        finish_axes(ax, f"{CASE_TITLE}: runtime comparison", "Number of cells, $N_xN_y$ [-]", "Runtime [s]", log_x=True, log_y=True)
        place_legend_outside(ax)
        fig.subplots_adjust(right=0.76)
        save_figure(fig, output_directory, "comparison_mesh_convergence_runtime", created_files)
    else:
        plt.close(fig)


def plot_mesh_convergence(case_root: Path, created_files: list[Path]) -> None:
    datasets: list[tuple[str, pd.DataFrame]] = []

    p1_path = case_root / P1_MESH_TEST / "mesh_convergence.csv"
    p1_data = plot_mesh_convergence_csv(p1_path, "P1", created_files)
    if p1_data is not None:
        datasets.append(("P1", p1_data))

    dom_root = case_root / DOM_MESH_TEST
    if dom_root.exists():
        for ordinate_directory in sorted(dom_root.glob("ordinate*"), key=ordinate_number):
            order = ordinate_number(ordinate_directory)
            dataframe = plot_mesh_convergence_csv(ordinate_directory / "mesh_convergence.csv", f"DOM S{order}", created_files)
            if dataframe is not None:
                datasets.append((f"DOM S{order}", dataframe))

    plot_combined_mesh_convergence(case_root, datasets, created_files)


def plot_angular_convergence(case_root: Path, created_files: list[Path]) -> None:
    angular_root = case_root / DOM_ANGULAR_TEST
    path = angular_root / "angular_convergence.csv"
    if not path.exists():
        warn_missing(path)
        return

    dataframe = read_standard_csv(path)
    output_directory = angular_root / "Plots"
    ordinate = pd.to_numeric(dataframe["ordinate_order"], errors="coerce").to_numpy(dtype=float)
    relative_change = pd.to_numeric(dataframe["successive_ordinate_relative_L2_change"], errors="coerce").to_numpy(dtype=float)
    observed_order = pd.to_numeric(dataframe["observed_angular_order"], errors="coerce").to_numpy(dtype=float)
    runtime = pd.to_numeric(dataframe["runtime"], errors="coerce").to_numpy(dtype=float)

    mask = np.isfinite(relative_change) & (relative_change > 0.0)
    if np.any(mask):
        fig, ax = plt.subplots()
        ax.plot(ordinate[mask], relative_change[mask], marker="o")
        finish_axes(ax, f"{CASE_TITLE}: DOM angular convergence", "Ordinate order, $S_N$ [-]", "Successive-ordinate relative $L_2$ change [-]", log_y=True)
        ax.set_xticks(ordinate)
        save_figure(fig, output_directory, "angular_convergence_relative_L2_change", created_files)

    mask = np.isfinite(observed_order)
    if np.any(mask):
        fig, ax = plt.subplots()
        ax.plot(ordinate[mask], observed_order[mask], marker="o")
        finish_axes(ax, f"{CASE_TITLE}: observed angular order", "Ordinate order, $S_N$ [-]", "Observed angular order [-]")
        ax.set_xticks(ordinate)
        save_figure(fig, output_directory, "angular_convergence_observed_order", created_files)

    mask = np.isfinite(runtime) & (runtime > 0.0)
    if np.any(mask):
        fig, ax = plt.subplots()
        ax.plot(ordinate[mask], runtime[mask], marker="o")
        finish_axes(ax, f"{CASE_TITLE}: DOM angular runtime", "Ordinate order, $S_N$ [-]", "Runtime [s]", log_y=True)
        ax.set_xticks(ordinate)
        save_figure(fig, output_directory, "angular_convergence_runtime", created_files)

    fig, ax = plt.subplots()
    plotted = False
    for order in SUPPORTED_ORDINATES:
        centerline_path = angular_root / f"ordinate{order}" / "320x320x1" / "centerlineG.csv"
        if not centerline_path.exists():
            continue
        _, centerline = read_metadata_csv(centerline_path)
        x_column = first_existing_column(centerline, ("x", "x_over_L"))
        g_column = first_existing_column(centerline, ("G_DOM", "G"))
        x = centerline[x_column].to_numpy(dtype=float)
        if x_column == "x":
            x = x / DOMAIN_LENGTH
        ax.plot(x, centerline[g_column].to_numpy(dtype=float), label=f"DOM S{order}")
        plotted = True

    if plotted:
        finish_axes(ax, f"{CASE_TITLE}: DOM angular centreline comparison", "Normalised coordinate, $x/L$ [-]", "Radiative energy, $G$ [W m$^{-2}$]")
        ax.set_xlim(0.0, 1.0)
        place_legend_outside(ax)
        fig.subplots_adjust(right=0.76)
        save_figure(fig, output_directory, "angular_convergence_centerline_G", created_files)
    else:
        plt.close(fig)

def _read_n_floats(stream, count: int) -> np.ndarray:
    values: list[float] = []
    while len(values) < count:
        line = stream.readline()
        if not line:
            raise ValueError("Unexpected end of VTK file")
        values.extend(float(value) for value in line.split())
    return np.asarray(values[:count], dtype=float)


def read_p1_vtk(path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    dimensions = None
    x_nodes = None
    y_nodes = None
    z_nodes = None
    number_of_cells = None
    fields: dict[str, np.ndarray] = {}

    with path.open("r") as stream:
        while True:
            line = stream.readline()
            if not line:
                break
            stripped = line.strip()
            if stripped.startswith("DIMENSIONS "):
                dimensions = tuple(int(value) for value in stripped.split()[1:4])
            elif stripped.startswith("X_COORDINATES "):
                count = int(stripped.split()[1])
                x_nodes = _read_n_floats(stream, count)
            elif stripped.startswith("Y_COORDINATES "):
                count = int(stripped.split()[1])
                y_nodes = _read_n_floats(stream, count)
            elif stripped.startswith("Z_COORDINATES "):
                count = int(stripped.split()[1])
                z_nodes = _read_n_floats(stream, count)
            elif stripped.startswith("CELL_DATA "):
                number_of_cells = int(stripped.split()[1])
            elif stripped.startswith("SCALARS "):
                if number_of_cells is None:
                    raise ValueError(f"CELL_DATA must appear before SCALARS in {path}")
                field_name = stripped.split()[1]
                lookup_line = stream.readline().strip()
                if not lookup_line.startswith("LOOKUP_TABLE"):
                    raise ValueError(f"Missing LOOKUP_TABLE after {field_name} in {path}")
                fields[field_name] = _read_n_floats(stream, number_of_cells)

    if dimensions is None or x_nodes is None or y_nodes is None or z_nodes is None:
        raise ValueError(f"Incomplete rectilinear-grid metadata in {path}")

    nx, ny, nz = dimensions[0]-1, dimensions[1]-1, dimensions[2]-1
    x = 0.5*(x_nodes[:-1]+x_nodes[1:])
    y = 0.5*(y_nodes[:-1]+y_nodes[1:])
    two_dimensional_fields = {name: values.reshape(nz, ny, nx)[nz//2] for name, values in fields.items()}
    return x, y, two_dimensional_fields


def ordinate_directions_and_weights(order: int) -> tuple[np.ndarray, np.ndarray]:
    base_sets = {
        2: ((0.5000000,0.7071068,0.5000000,1.5707963),),
        4: ((0.2958759,0.2958759,0.9082483,0.5235987),(0.2958759,0.9082483,0.2958759,0.5235987),(0.9082483,0.2958759,0.2958759,0.5235987)),
        6: ((0.1838670,0.1838670,0.9656013,0.1609517),(0.1838670,0.6950514,0.6950514,0.3626469),(0.1838670,0.9656013,0.1838670,0.1609517),(0.6950514,0.1838670,0.6950514,0.3626469),(0.6950514,0.6950514,0.1838670,0.3626469),(0.9656013,0.1838670,0.1838670,0.1609517)),
        8: ((0.1422555,0.1422555,0.9795543,0.1712359),(0.1422555,0.5773503,0.8040087,0.0992284),(0.1422555,0.8040087,0.5773503,0.0992284),(0.1422555,0.9795543,0.1422555,0.1712359),(0.5773503,0.1422555,0.8040087,0.0992284),(0.5773503,0.5773503,0.5773503,0.4617179),(0.5773503,0.8040087,0.1422555,0.0992284),(0.8040087,0.1422555,0.5773503,0.0992284),(0.8040087,0.5773503,0.1422555,0.0992284),(0.9795543,0.1422555,0.1422555,0.1712359)),
    }
    if order not in base_sets:
        raise ValueError(f"Unsupported ordinate order S{order}")

    directions: list[tuple[float, float, float]] = []
    weights: list[float] = []
    for x, y, z, weight in base_sets[order]:
        directions.extend(((x,y,z),(-x,y,z),(x,-y,z),(x,y,-z),(-x,-y,z),(x,-y,-z),(-x,y,-z),(-x,-y,-z)))
        weights.extend((weight,)*8)
    return np.asarray(directions, dtype=float), np.asarray(weights, dtype=float)


def read_dom_vtk(path: Path, order: int) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    directions, weights = ordinate_directions_and_weights(order)
    dimensions = None
    x_nodes = None
    y_nodes = None
    z_nodes = None
    number_of_cells = None
    g = None
    qx = None
    qy = None
    qz = None
    vector_index = 0

    with path.open("r") as stream:
        while True:
            line = stream.readline()
            if not line:
                break
            stripped = line.strip()
            if stripped.startswith("DIMENSIONS "):
                dimensions = tuple(int(value) for value in stripped.split()[1:4])
            elif stripped.startswith("X_COORDINATES "):
                x_nodes = _read_n_floats(stream, int(stripped.split()[1]))
            elif stripped.startswith("Y_COORDINATES "):
                y_nodes = _read_n_floats(stream, int(stripped.split()[1]))
            elif stripped.startswith("Z_COORDINATES "):
                z_nodes = _read_n_floats(stream, int(stripped.split()[1]))
            elif stripped.startswith("CELL_DATA "):
                number_of_cells = int(stripped.split()[1])
                g = np.zeros(number_of_cells)
                qx = np.zeros(number_of_cells)
                qy = np.zeros(number_of_cells)
                qz = np.zeros(number_of_cells)
            elif stripped.startswith("VECTORS "):
                if number_of_cells is None or g is None or qx is None or qy is None or qz is None:
                    raise ValueError(f"CELL_DATA must appear before VECTORS in {path}")
                if vector_index >= len(weights):
                    raise ValueError(f"More ordinate vectors than expected for S{order} in {path}")
                components = _read_n_floats(stream, 3*number_of_cells).reshape(number_of_cells, 3)
                direction = directions[vector_index]
                direction_norm_squared = float(np.dot(direction, direction))
                intensity = components @ direction / direction_norm_squared
                weight = weights[vector_index]
                g += weight*intensity
                qx += weight*components[:,0]
                qy += weight*components[:,1]
                qz += weight*components[:,2]
                vector_index += 1

    if vector_index != len(weights):
        raise ValueError(f"Expected {len(weights)} ordinate vectors for S{order}, found {vector_index} in {path}")
    if dimensions is None or x_nodes is None or y_nodes is None or z_nodes is None or g is None or qx is None or qy is None or qz is None:
        raise ValueError(f"Incomplete DOM VTK file: {path}")

    nx, ny, nz = dimensions[0]-1, dimensions[1]-1, dimensions[2]-1
    x = 0.5*(x_nodes[:-1]+x_nodes[1:])
    y = 0.5*(y_nodes[:-1]+y_nodes[1:])
    fields = {
        "G": g.reshape(nz,ny,nx)[nz//2],
        "qx": qx.reshape(nz,ny,nx)[nz//2],
        "qy": qy.reshape(nz,ny,nx)[nz//2],
        "qz": qz.reshape(nz,ny,nx)[nz//2],
    }
    return x, y, fields


def latest_dom_vtk(directory: Path) -> Path | None:
    simulation_directory = directory / "SimulationFolder"
    files = list(simulation_directory.glob("simulation*.vtk"))
    if not files:
        return None

    def simulation_number(path: Path) -> int:
        match = re.search(r"simulation(\d+)\.vtk$", path.name)
        return int(match.group(1)) if match else -1

    return max(files, key=simulation_number)


def temperature_field(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    X, Y = np.meshgrid(x, y)
    return np.full_like(X, MEDIUM_TEMPERATURE, dtype=float)


def calculate_derived_fields(x: np.ndarray, y: np.ndarray, g: np.ndarray, temperature: np.ndarray, model: str, qx: np.ndarray | None = None, qy: np.ndarray | None = None) -> dict[str, np.ndarray]:
    if model == "P1":
        gamma = 1.0/(3.0*(ABSORPTION_COEFFICIENT+SCATTERING_COEFFICIENT))
        dy = float(np.mean(np.diff(y)))
        dx = float(np.mean(np.diff(x)))
        dGdy, dGdx = np.gradient(g, dy, dx, edge_order=2)
        qx = -gamma*dGdx
        qy = -gamma*dGdy
    elif qx is None or qy is None:
        raise ValueError("DOM flux components are required")

    qmag = np.sqrt(qx*qx+qy*qy)
    source = ABSORPTION_COEFFICIENT*(4.0*SIGMA*np.power(temperature,4)-g)
    return {"G":g, "qx":qx, "qy":qy, "qmag":qmag, "source":source, "temperature":temperature}


def _collect_profile_contour_fields(case_root: Path) -> list[tuple[str, np.ndarray, np.ndarray, dict[str, np.ndarray]]]:
    datasets: list[tuple[str, np.ndarray, np.ndarray, dict[str, np.ndarray]]] = []

    p1_directory = case_root / P1_PROFILE_TEST / "320x320x1"
    p1_vtk = p1_directory / "simulation.vtk"
    if p1_vtk.exists():
        x, y, vtk_fields = read_p1_vtk(p1_vtk)
        if "G" not in vtk_fields:
            raise ValueError(f"G was not found in {p1_vtk}")
        temperature = vtk_fields.get("temperature", temperature_field(x, y))
        datasets.append(("P1", x, y, calculate_derived_fields(x, y, vtk_fields["G"], temperature, "P1")))
    else:
        warn_missing(p1_vtk)

    dom_root = case_root / DOM_PROFILE_TEST
    for order in SUPPORTED_ORDINATES:
        directory = dom_root / f"ordinate{order}" / "320x320x1"
        vtk_path = latest_dom_vtk(directory)
        if vtk_path is None:
            warn_missing(directory / "SimulationFolder" / "simulation*.vtk")
            continue
        x, y, vtk_fields = read_dom_vtk(vtk_path, order)
        temperature = temperature_field(x, y)
        fields = calculate_derived_fields(
            x, y, vtk_fields["G"], temperature, "DOM", vtk_fields["qx"], vtk_fields["qy"]
        )
        datasets.append((f"DOM $S_{order}$", x, y, fields))

    return datasets


def _model_filename(label: str) -> str:
    if label == "P1":
        return "P1"
    match = re.search(r"S_(\d+)", label)
    return f"DOM_S{match.group(1)}" if match else re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_")


def plot_individual_log_contours(
    datasets: Sequence[tuple[str, np.ndarray, np.ndarray, dict[str, np.ndarray]]],
    field_name: str,
    quantity_title: str,
    colorbar_label: str,
    output_directory: Path,
    stem_prefix: str,
    created_files: list[Path],
) -> None:
    if not datasets:
        return

    positive_values = []
    for _, _, _, fields in datasets:
        values = np.asarray(fields[field_name], dtype=float)
        finite_positive = values[np.isfinite(values) & (values > 0.0)]
        if finite_positive.size:
            positive_values.append(finite_positive)

    if not positive_values:
        raise ValueError(
            f"No positive finite values are available for logarithmic {field_name} contours"
        )

    all_positive = np.concatenate(positive_values)
    global_min = float(np.min(all_positive))
    global_max = float(np.max(all_positive))
    if not np.isfinite(global_min) or not np.isfinite(global_max) or global_max <= 0.0:
        raise ValueError(f"Invalid values encountered in {field_name} contour data")

    if field_name == "G":
        levels = np.asarray([
            1.0e-1, 3.16228e-1, 1.0, 3.16228, 10.0, 31.6228,
            100.0, 316.228, 1000.0, 3162.28, 10000.0,
        ])
    elif field_name == "qmag":
        levels = np.asarray([
            1.0e-1, 2.83525e-1, 8.32553e-1, 2.40225, 6.93545,
            20.0, 57.708, 166.511, 480.45, 1386.29, 4000.0,
        ])
    else:
        if math.isclose(global_min, global_max):
            global_min = global_max / 10.0
        levels = np.geomspace(global_min, global_max, 11)

    norm = LogNorm(vmin=float(levels[0]), vmax=float(levels[-1]))

    for label, x, y, fields in datasets:
        values = np.asarray(fields[field_name], dtype=float)
        masked_values = np.ma.masked_where(~np.isfinite(values) | (values <= 0.0), values)

        X, Y = np.meshgrid(x, y)

        fig, ax = plt.subplots(figsize=(7.8, 6.5))
        contour = ax.contourf(
            X,
            Y,
            masked_values,
            levels=levels,
            norm=norm,
            extend="both",
        )

        ax.contour(
            X,
            Y,
            masked_values,
            levels=levels,
            norm=norm,
            linewidths=0.45,
            colors="0.25",
            alpha=0.75,
        )
        ax.set_title(f"{label}", fontweight="bold", pad=12)
        ax.set_xlabel("$x$ [m]")
        ax.set_ylabel("$y$ [m]")
        ax.set_xlim(0.0, DOMAIN_LENGTH)
        ax.set_ylim(0.0, DOMAIN_LENGTH)
        ax.set_aspect("equal")
        ax.grid(False)

        colorbar = fig.colorbar(contour, ax=ax, pad=0.03, ticks=levels)
        colorbar.set_label(colorbar_label, labelpad=12)
        colorbar.ax.tick_params(labelsize=11)
        colorbar.ax.set_yticklabels([f"{value:g}" for value in levels])
        fig.subplots_adjust(top=0.88, bottom=0.14, left=0.12, right=0.84)

        model_name = _model_filename(label)
        save_figure(
            fig,
            output_directory,
            f"{stem_prefix}_{model_name}",
            created_files,
        )


def plot_profile_contours(case_root: Path, created_files: list[Path]) -> None:
    try:
        datasets = _collect_profile_contour_fields(case_root)
    except Exception as error:
        warnings.warn(f"Could not reconstruct contour fields: {error}")
        return

    expected_labels = {"P1", "DOM $S_2$", "DOM $S_4$", "DOM $S_6$", "DOM $S_8$"}
    found_labels = {label for label, _, _, _ in datasets}
    missing = expected_labels - found_labels
    if missing:
        warnings.warn("Contour comparison is incomplete; missing: " + ", ".join(sorted(missing)))

    output_directory = case_root / "Plots"
    plot_individual_log_contours(
        datasets,
        "G",
        f"",
        "Irradiation, $G$ [W m$^{-2}$]",
        output_directory,
        "irradiation_contour",
        created_files,
    )
    plot_individual_log_contours(
        datasets,
        "qmag",
        f"",
        "Heat-flux magnitude, $|\\mathbf{q}_r|$ [W m$^{-2}$]",
        output_directory,
        "heat_flux_magnitude_contour",
        created_files,
    )


def plot_p1_mesh_profile_overlay(case_root: Path, created_files: list[Path]) -> None:
    mesh_root = case_root / P1_MESH_TEST
    fig, ax = plt.subplots()
    plotted = False
    for mesh_size in SUPPORTED_MESHES:
        vtk_path = mesh_root / f"{mesh_size}x{mesh_size}x1" / "simulation.vtk"
        if not vtk_path.exists():
            continue
        try:
            x, y, fields = read_p1_vtk(vtk_path)
            if "G" not in fields:
                continue
            centreline_index = int(np.argmin(np.abs(y-0.5*DOMAIN_LENGTH)))
            ax.plot(x/DOMAIN_LENGTH, fields["G"][centreline_index,:], label=f"{mesh_size}x{mesh_size}")
            plotted = True
        except Exception as error:
            warnings.warn(f"Could not read {vtk_path}: {error}")
    if plotted:
        finish_axes(ax, f"{CASE_TITLE}: P1 mesh-centreline convergence", "Normalised coordinate, $x/L$ [-]", "Radiative energy, $G$ [W m$^{-2}$]")
        ax.set_xlim(0.0,1.0)
        place_legend_outside(ax)
        fig.subplots_adjust(right=0.76)
        save_figure(fig, mesh_root / "Plots", "mesh_convergence_centerline_G", created_files)
    else:
        plt.close(fig)


def plot_dom_mesh_profile_overlays(case_root: Path, created_files: list[Path]) -> None:
    mesh_root = case_root / DOM_MESH_TEST
    for order in SUPPORTED_ORDINATES:
        ordinate_root = mesh_root / f"ordinate{order}"
        fig, ax = plt.subplots()
        plotted = False
        for mesh_size in SUPPORTED_MESHES:
            directory = ordinate_root / f"{mesh_size}x{mesh_size}x1"
            vtk_path = latest_dom_vtk(directory)
            if vtk_path is None:
                continue
            try:
                x, y, fields = read_dom_vtk(vtk_path, order)
                centreline_index = int(np.argmin(np.abs(y-0.5*DOMAIN_LENGTH)))
                ax.plot(x/DOMAIN_LENGTH, fields["G"][centreline_index,:], label=f"{mesh_size}x{mesh_size}")
                plotted = True
            except Exception as error:
                warnings.warn(f"Could not read {vtk_path}: {error}")
        if plotted:
            finish_axes(ax, f"{CASE_TITLE}: DOM S{order} mesh-centreline convergence", "Normalised coordinate, $x/L$ [-]", "Radiative energy, $G$ [W m$^{-2}$]")
            ax.set_xlim(0.0,1.0)
            place_legend_outside(ax)
            fig.subplots_adjust(right=0.76)
            save_figure(fig, ordinate_root / "Plots", "mesh_convergence_centerline_G", created_files)
        else:
            plt.close(fig)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=f"Create all report plots for {CASE_TITLE}.")
    parser.add_argument("--output-root", type=Path, default=None, help="Path to the project Output directory.")
    parser.add_argument("--skip-vtk", action="store_true", help="Skip VTK contour reconstruction and VTK centreline overlays.")
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    project_root = find_project_root()
    output_root = arguments.output_root.resolve() if arguments.output_root is not None else project_root / "Output"
    case_root = output_root / PARENT_FOLDER

    if not case_root.exists():
        raise FileNotFoundError(
            f"Case output directory was not found: {case_root}\n"
            f"Run the C++ tests first or pass --output-root /path/to/Output."
        )

    apply_report_style()
    created_files: list[Path] = []

    plot_all_profiles(case_root, created_files)
    plot_mesh_convergence(case_root, created_files)
    plot_angular_convergence(case_root, created_files)

    if not arguments.skip_vtk:
        plot_profile_contours(case_root, created_files)
        plot_p1_mesh_profile_overlay(case_root, created_files)
        plot_dom_mesh_profile_overlays(case_root, created_files)

    print(f"Created {len(created_files)} plot files for {CASE_TITLE}.")
    for path in created_files:
        print(path)


if __name__ == "__main__":
    main()
