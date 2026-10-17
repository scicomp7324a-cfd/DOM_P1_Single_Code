#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


DEFAULT_OUTPUT_ROOT = Path("Output")
PARENT_CASE_FOLDER = "DOMPreliminaryValidation"
DEFAULT_FIG_DIR_NAME = "Plots"

CASE_CLOSED_CANDIDATES = (
    "ParallelPlateNoMediumVacuumClosedDOM",
    "ParallelPlateNoMediumVacuumClosed",
)
CASE_PARTICIPATING_CANDIDATES = (
    "ParallelPlateParticipatingMediumClosedDOM",
    "ParallelPlateParticipatingMediumClosed",
)
CASE_CONDUCTION_CANDIDATES = (
    "SquarePlateRadiationConductionClosedDOM",
    "SquarePlateRadiationConductionClosed",
    "SquarePlateRadiationConduction",
)

CLOSED_TEMPERATURE_CASES = ["300", "400", "500", "600", "700", "800", "900", "1000"]
PARTICIPATING_OPTICAL_LENGTHS = ["0.1", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "4.0", "5.0", "7.0", "9.0", "10.0"]
PARTICIPATING_ANALYTICAL_TAUS = ["0.1", "1.0", "2.0", "10.0"]

CONDUCTION_CASE_FOLDERS: list[str] = []
ORDINATES = (2, 4, 6, 8)

ORDINATE_STYLES = {
    2: {"color": "#D81B60", "linestyle": "-",  "marker": "o"},
    4: {"color": "#1E88E5", "linestyle": "--", "marker": "s"},
    6: {"color": "#00A087", "linestyle": "-.", "marker": "^"},
    8: {"color": "#7B2CBF", "linestyle": ":",  "marker": "D"},
}



def apply_report_style() -> None:
    plt.rcParams.update({
        "figure.figsize": (10.5, 7.2),
        "figure.dpi": 150,
        "savefig.dpi": 400,
        "font.size": 17,
        "axes.titlesize": 21,
        "axes.labelsize": 19,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 14,
        "lines.linewidth": 2.8,
        "lines.markersize": 7.5,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linestyle": "--",
        "axes.linewidth": 1.25,
        "legend.frameon": True,
        "legend.framealpha": 0.95,
        "legend.edgecolor": "0.75",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def finish_axes(ax: plt.Axes, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title, pad=14, fontweight="bold")
    ax.set_xlabel(xlabel, labelpad=10)
    ax.set_ylabel(ylabel, labelpad=10)
    ax.tick_params(axis="both", which="major", length=7, width=1.15)
    ax.tick_params(axis="both", which="minor", length=4, width=0.9)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    for spine in ax.spines.values():
        spine.set_linewidth(1.25)


def save_figure(fig: plt.Figure, fig_dir: Path, stem: str) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        out = fig_dir / f"{stem}.{ext}"
        fig.savefig(out, bbox_inches="tight")
    plt.close(fig)



VIBRANT_COLORS = (
    "#D81B60",  # magenta
    "#1E88E5",  # blue
    "#FFC107",  # amber
    "#004D40",  # dark teal
    "#7B1FA2",  # purple
    "#E64A19",  # deep orange
    "#00ACC1",  # cyan
    "#43A047",  # green
)

VIBRANT_MARKERS = ("o", "s", "^", "D", "v", "P", "X", "<")
VIBRANT_LINESTYLES = ("-", "--", "-.", ":", "-", "--", "-.", ":")


def place_external_legend(ax: plt.Axes, ncols: int = 1) -> None:
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        ncols=ncols,
        frameon=True,
        handlelength=3.0,
        labelspacing=0.65,
    )


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    return df


def read_standard_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return clean_columns(pd.read_csv(path))


def read_centerline_csv(path: Path) -> pd.DataFrame:

    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    lines = path.read_text().splitlines()
    header_idx = None
    header_tokens = (
        "x_over_L", "x_centerline_coordinate", "z_over_L",
        "y/l", "x/l"
    )

    for i, line in enumerate(lines):
        stripped = line.strip()
        if any(stripped.startswith(tok) for tok in header_tokens):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(f"Could not find centerline header in: {path}")

    df = pd.read_csv(path, skiprows=header_idx)
    df = clean_columns(df)
    df = df.dropna(how="all")
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(how="any")
    return df



def read_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            break
        parts = [part.strip() for part in line.split(",", 1)]
        if len(parts) == 2:
            metadata[parts[0]] = parts[1]
    return metadata


def ordinate_from_path_or_metadata(path: Path) -> int:
    for part in reversed(path.parts):
        match = re.fullmatch(r"S(2|4|6|8)", part, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))

    metadata = read_metadata(path)
    for key in ("ordinateS", "ordinate", "ordinate_order"):
        if key in metadata:
            value = metadata[key].strip().upper().removeprefix("S")
            return int(float(value))

    raise ValueError(f"Could not determine ordinate order for {path}")


def files_for_subcase_and_ordinates(
    case_dir: Path,
    subcase: str,
    filename: str,
) -> dict[int, Path]:
    found: dict[int, Path] = {}
    subcase_dir = case_dir / subcase

    search_root = subcase_dir if subcase_dir.is_dir() else case_dir
    for path in search_root.rglob(filename):
        try:
            ordinate = ordinate_from_path_or_metadata(path)
        except (ValueError, OSError):
            continue
        if ordinate in ORDINATES:
            found[ordinate] = path

    return dict(sorted(found.items()))


def normalise_ordinate_column(df: pd.DataFrame) -> pd.DataFrame:
    if "ordinate" not in df.columns:
        raise KeyError(
            "The diagnostics CSV has no 'ordinate' column. "
            "Regenerate the DOM output using the updated validation functions."
        )
    result = df.copy()
    result["ordinate"] = pd.to_numeric(result["ordinate"], errors="raise").astype(int)
    return result


def style_for_ordinate(ordinate: int) -> dict[str, str]:
    if ordinate not in ORDINATE_STYLES:
        raise KeyError(f"Unsupported ordinate order S{ordinate}")
    return ORDINATE_STYLES[ordinate]


def first_existing_column(df: pd.DataFrame, candidates: Sequence[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"None of these columns found: {candidates}. Available columns: {list(df.columns)}")


def sorted_case_folders(case_dir: Path) -> list[str]:
    folders = [p.name for p in case_dir.iterdir() if p.is_dir()]
    return sorted(folders, key=lambda s: float(s))


def resolve_validation_root(path: Path) -> Path:

    path = path.expanduser().resolve()

    candidates = [
        path,
        path / PARENT_CASE_FOLDER,
        path / "Output" / PARENT_CASE_FOLDER,
    ]

    for candidate in candidates:
        if candidate.is_dir() and any(
            (candidate / name).is_dir()
            for name in (
                *CASE_CLOSED_CANDIDATES,
                *CASE_PARTICIPATING_CANDIDATES,
                *CASE_CONDUCTION_CANDIDATES,
            )
        ):
            return candidate

    expected = path / "Output" / PARENT_CASE_FOLDER
    raise FileNotFoundError(
        "Could not locate the DOM preliminary-validation output tree. "
        f"Expected a folder such as: {expected}"
    )


def resolve_case_directory(validation_root: Path, candidates: Sequence[str]) -> Path:
    for name in candidates:
        candidate = validation_root / name
        if candidate.is_dir():
            return candidate

    wanted = {name.lower() for name in candidates}
    for child in validation_root.iterdir():
        if child.is_dir() and (
            child.name.lower() in wanted
            or any(child.name.lower().startswith(name.lower()) for name in candidates)
        ):
            return child

    raise FileNotFoundError(
        f"None of the expected case folders exists under {validation_root}: "
        + ", ".join(candidates)
    )


def find_output_file(case_dir: Path, filename: str, subcase: str | None = None) -> Path:

    direct_candidates: list[Path] = []
    if subcase is not None:
        direct_candidates.extend([
            case_dir / subcase / filename,
            case_dir / filename,
        ])
    else:
        direct_candidates.append(case_dir / filename)

    for candidate in direct_candidates:
        if candidate.exists():
            return candidate

    matches = list(case_dir.rglob(filename))
    if subcase is not None:
        exact = [p for p in matches if subcase in p.parts]
        if exact:
            return sorted(exact, key=lambda p: (len(p.parts), str(p)))[0]

    if matches:
        return sorted(matches, key=lambda p: (len(p.parts), str(p)))[0]

    suffix = f" for subcase {subcase}" if subcase is not None else ""
    raise FileNotFoundError(f"Missing {filename}{suffix} under {case_dir}")


def optional_output_file(case_dir: Path, filename: str, subcase: str | None = None) -> Path | None:
    try:
        return find_output_file(case_dir, filename, subcase)
    except FileNotFoundError:
        return None


def plot_closed_heat_flux_vs_east_temperature(output_root: Path, fig_dir: Path) -> None:
    case_dir = resolve_case_directory(output_root, CASE_CLOSED_CANDIDATES)
    df = normalise_ordinate_column(
        read_standard_csv(find_output_file(case_dir, "heat_diagnostics.csv"))
    )

    fig, ax = plt.subplots()

    analytical_drawn = False
    for ordinate in ORDINATES:
        group = df[df["ordinate"] == ordinate].copy()
        if group.empty:
            print(f"Warning: no closed-medium heat-flux data found for S{ordinate}")
            continue

        group = group.sort_values("east_wall_temperature")
        style = style_for_ordinate(ordinate)

        ax.plot(
            group["east_wall_temperature"].astype(float),
            np.abs(group["west_wall_heat_flux"].astype(float)),
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=2.8,
            marker=style["marker"],
            markersize=7.0,
            markerfacecolor="white",
            markeredgecolor=style["color"],
            markeredgewidth=1.4,
            label=f"DOM $S_{{{ordinate}}}$",
        )

        if not analytical_drawn:
            ax.plot(
                group["east_wall_temperature"].astype(float),
                np.abs(group["west_wall_analytical_heat_flux"].astype(float)),
                color="black",
                linestyle="None",
                marker="X",
                markersize=8.0,
                markerfacecolor="black",
                label="Analytical",
            )
            analytical_drawn = True

    finish_axes(
        ax,
        title="Closed Parallel Plates: Heat Flux vs East-Wall Temperature",
        xlabel="East-wall temperature, $T_E$ [K]",
        ylabel="Absolute west-wall heat flux, $|q_W|$ [W m$^{-2}$]",
    )
    place_external_legend(ax)
    fig.subplots_adjust(right=0.76)
    save_figure(fig, fig_dir, "DOM_plot1_closed_heat_flux_all_ordinates")


def plot_closed_centerline_G_profiles(output_root: Path, fig_dir: Path) -> None:

    case_dir = resolve_case_directory(output_root, CASE_CLOSED_CANDIDATES)

    for temp_case in CLOSED_TEMPERATURE_CASES:
        files = files_for_subcase_and_ordinates(
            case_dir, temp_case, "centerlineGInX.csv"
        )
        if not files:
            print(
                "Warning: no centreline irradiation files found for "
                f"east-wall temperature {temp_case} K"
            )
            continue

        fig, ax = plt.subplots()

        for ordinate, path in files.items():
            df = read_centerline_csv(path)
            x_col = first_existing_column(
                df, ["x_over_L", "x_centerline_coordinate"]
            )
            g_col = first_existing_column(df, ["G"])
            x = df[x_col].astype(float)

            style = style_for_ordinate(ordinate)
            ax.plot(
                x,
                df[g_col].astype(float),
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=2.8,
                marker=style["marker"],
                markersize=6.5,
                markerfacecolor="white",
                markeredgecolor=style["color"],
                markeredgewidth=1.3,
                markevery=max(1, len(df) // 10),
                label=f"DOM $S_{{{ordinate}}}$",
            )

        finish_axes(
            ax,
            title=(
                "Closed Parallel Plates: Centreline Irradiation "
                f"($T_E={temp_case}$ K)"
            ),
            xlabel="Normalised coordinate, $x/L$ [-]",
            ylabel="Irradiation, $G$ [W m$^{-2}$]",
        )
        ax.legend(loc="best")
        save_figure(
            fig,
            fig_dir,
            f"DOM_plot2_closed_centerline_G_TE_{temp_case}_all_ordinates",
        )


def plot_participating_dimensionless_heat_flux(output_root: Path, fig_dir: Path) -> None:
    case_dir = resolve_case_directory(output_root, CASE_PARTICIPATING_CANDIDATES)
    sim = normalise_ordinate_column(
        read_standard_csv(find_output_file(case_dir, "heat_diagnostics.csv"))
    )
    analytical = read_standard_csv(
        find_output_file(case_dir, "analytical_dimensionless_heat_flux.csv")
    )

    fig, ax = plt.subplots()

    ax.plot(
        analytical["absorption_coefficient"].astype(float),
        np.abs(analytical["dimension_heat_flux_east"].astype(float)),
        color="black",
        linestyle="None",
        marker="X",
        markersize=8.0,
        markerfacecolor="black",
        label="Analytical",
    )

    for ordinate in ORDINATES:
        group = sim[sim["ordinate"] == ordinate].copy()
        if group.empty:
            print(f"Warning: no participating-medium heat-flux data for S{ordinate}")
            continue
        group = group.sort_values("absorption_coefficient")
        style = style_for_ordinate(ordinate)

        ax.plot(
            group["absorption_coefficient"].astype(float),
            np.abs(group["dimensionless_heat_flux_west"].astype(float)),
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=2.8,
            marker=style["marker"],
            markersize=7.0,
            markerfacecolor="white",
            markeredgecolor=style["color"],
            markeredgewidth=1.3,
            label=f"DOM $S_{{{ordinate}}}$",
        )

    finish_axes(
        ax,
        title="Participating Medium: Dimensionless Heat Flux",
        xlabel="Optical thickness, $\\tau_L=aL$ [-]",
        ylabel="Dimensionless heat flux, $|\\psi|$ [-]",
    )
    place_external_legend(ax)
    fig.subplots_adjust(right=0.76)
    save_figure(
        fig,
        fig_dir,
        "DOM_plot3_participating_dimensionless_heat_flux_all_ordinates",
    )


def plot_participating_dimensionless_centerline_temperature(
    output_root: Path,
    fig_dir: Path,
) -> None:

    case_dir = resolve_case_directory(
        output_root,
        CASE_PARTICIPATING_CANDIDATES,
    )

    analytical = read_standard_csv(
        find_output_file(
            case_dir,
            "analytical_dimensionless_centerline_temperature.csv",
        )
    )

    x_ana_col = first_existing_column(
        analytical,
        ["x/l", "x_over_L"],
    )

    for tau in PARTICIPATING_ANALYTICAL_TAUS:
        files = files_for_subcase_and_ordinates(
            case_dir,
            tau,
            "dimensionlessCenterlineT.csv",
        )

        if not files:
            print(
                "Warning: no dimensionless centreline files found "
                f"for tau={tau}"
            )
            continue

        fig, ax = plt.subplots(figsize=(10.5, 6.5))

        if tau in analytical.columns:
            ax.plot(
                analytical[x_ana_col].astype(float),
                analytical[tau].astype(float),
                color="black",
                linestyle="None",
                marker="X",
                markersize=8.0,
                markerfacecolor="black",
                markeredgecolor="black",
                clip_on=False,
                zorder=5,
                label="Siegel and Howell reference",
            )

        for ordinate, path in sorted(files.items()):
            sim = read_centerline_csv(path)

            x_col = first_existing_column(
                sim,
                ["x_over_L", "x/l"],
            )
            y_col = first_existing_column(
                sim,
                ["dimensionless_T"],
            )

            style = style_for_ordinate(ordinate)

            ax.plot(
                sim[x_col].astype(float),
                sim[y_col].astype(float),
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=2.8,
                marker=style["marker"],
                markersize=6.0,
                markerfacecolor="white",
                markeredgecolor=style["color"],
                markeredgewidth=1.2,
                markevery=max(1, len(sim) // 10),
                clip_on=False,
                zorder=3,
                label=f"DOM $S_{{{ordinate}}}$",
            )

        finish_axes(
            ax,
            title=(
                "Participating Medium: Temperature Profile "
                f"($\\tau_L={tau}$)"
            ),
            xlabel="Normalised coordinate, $x/L$ [-]",
            ylabel="Dimensionless temperature, $T/T_W$ [-]",
        )

        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)

        ax.set_xticks(np.linspace(0.0, 1.0, 6))
        ax.set_yticks(np.linspace(0.0, 1.0, 6))

        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.015, 1.0),
            borderaxespad=0.0,
            frameon=True,
            fontsize=12,
            handlelength=2.5,
        )

        fig.subplots_adjust(
            left=0.12,
            right=0.73,
            bottom=0.15,
            top=0.88,
        )

        save_figure(
            fig,
            fig_dir,
            (
                "DOM_plot4_participating_temperature_tau_"
                f"{tau.replace('.', 'p')}_all_ordinates"
            ),
        )

def _collect_conduction_cases(case_dir: Path) -> dict[float, dict[int, Path]]:

    diagnostics_path = find_output_file(case_dir, "heat_diagnostics.csv")
    diagnostics = normalise_ordinate_column(read_standard_csv(diagnostics_path))

    if "stark_number" not in diagnostics.columns:
        raise KeyError(
            "heat_diagnostics.csv has no 'stark_number' column. "
            "Regenerate this case using the updated radiation-conduction function."
        )

    conductivity_to_stark: dict[float, float] = {}
    for _, row in diagnostics.iterrows():
        conductivity_to_stark[float(row["conductivity"])] = float(row["stark_number"])

    grouped: dict[float, dict[int, Path]] = {}
    for path in case_dir.rglob("dimensionlessCenterlineTInZ.csv"):
        metadata = read_metadata(path)
        ordinate = ordinate_from_path_or_metadata(path)

        if "stark_number" in metadata:
            stark = float(metadata["stark_number"])
        elif "conductivity" in metadata:
            conductivity = float(metadata["conductivity"])
            stark = min(
                conductivity_to_stark,
                key=lambda value: abs(value - conductivity),
            )
            stark = conductivity_to_stark[stark]
        else:
            try:
                conductivity = float(path.parent.parent.name)
            except ValueError as exc:
                raise ValueError(
                    f"Cannot infer conductivity/Stark number from {path}"
                ) from exc
            nearest = min(
                conductivity_to_stark,
                key=lambda value: abs(value - conductivity),
            )
            stark = conductivity_to_stark[nearest]

        grouped.setdefault(stark, {})[ordinate] = path

    return grouped


def plot_conduction_dimensionless_temperature_z(
    output_root: Path,
    fig_dir: Path,
) -> None:

    case_dir = resolve_case_directory(output_root, CASE_CONDUCTION_CANDIDATES)
    analytical = read_standard_csv(
        find_output_file(
            case_dir,
            "analytical_dimensionless_centerline_temperature.csv",
        )
    )
    x_col = first_existing_column(analytical, ["y/l", "z_over_L"])

    analytical_columns = {
        0.05: "N0.05",
        0.1: "N0.1",
        1.0: "N1",
    }

    grouped = _collect_conduction_cases(case_dir)
    if not grouped:
        raise FileNotFoundError(
            f"No dimensionlessCenterlineTInZ.csv files found under {case_dir}"
        )

    for stark in (0.05, 0.1, 1.0):
        files = grouped.get(stark, {})
        if not files:
            print(f"Warning: no radiation-conduction files found for N={stark:g}")
            continue

        fig, ax = plt.subplots()
        analytical_col = analytical_columns[stark]

        if analytical_col in analytical.columns:
            ax.plot(
                analytical[x_col].astype(float),
                analytical[analytical_col].astype(float),
                color="black",
                linestyle="None",
                marker="X",
                markersize=8.0,
                markerfacecolor="black",
                label="Rousse (2000)",
            )

        for ordinate, path in sorted(files.items()):
            sim = read_centerline_csv(path)
            z_col = first_existing_column(sim, ["z_over_L"])
            t_col = first_existing_column(sim, ["dimensionless_T"])
            style = style_for_ordinate(ordinate)

            ax.plot(
                sim[z_col].astype(float),
                sim[t_col].astype(float),
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=2.8,
                marker=style["marker"],
                markersize=6.0,
                markerfacecolor="white",
                markeredgecolor=style["color"],
                markeredgewidth=1.2,
                markevery=max(1, len(sim) // 10),
                label=f"DOM $S_{{{ordinate}}}$",
            )

        finish_axes(
            ax,
            title=f"Radiation--Conduction Square Plate ($N={stark:g}$)",
            xlabel="Normalised vertical coordinate, $z/L$ [-]",
            ylabel="Dimensionless temperature, $T/T_W$ [-]",
        )
        ax.legend(loc="best")
        save_figure(
            fig,
            fig_dir,
            f"DOM_plot5_conduction_temperature_N_{str(stark).replace('.', 'p')}_all_ordinates",
        )


def make_all_plots(output_root: Path, fig_dir: Path) -> None:
    apply_report_style()

    plot_closed_heat_flux_vs_east_temperature(output_root, fig_dir)
    plot_closed_centerline_G_profiles(output_root, fig_dir)
    plot_participating_dimensionless_heat_flux(output_root, fig_dir)
    plot_participating_dimensionless_centerline_temperature(output_root, fig_dir)
    plot_conduction_dimensionless_temperature_z(output_root, fig_dir)

    print(f"Saved all-ordinate DOM preliminary-validation figures to: {fig_dir.resolve()}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate report-quality DOM preliminary-validation plots from CSV output files.")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Project root, Output directory, or Output/DOMPreliminaryValidation directory.",
    )
    parser.add_argument(
        "--fig-dir",
        type=Path,
        default=None,
        help="Directory where figures should be saved. Default: Output/DOMPreliminaryValidation/Plots",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    validation_root = resolve_validation_root(args.output_root)
    fig_dir = (
        args.fig_dir.expanduser().resolve()
        if args.fig_dir
        else validation_root / DEFAULT_FIG_DIR_NAME
    )

    make_all_plots(validation_root, fig_dir)


if __name__ == "__main__":
    main()
