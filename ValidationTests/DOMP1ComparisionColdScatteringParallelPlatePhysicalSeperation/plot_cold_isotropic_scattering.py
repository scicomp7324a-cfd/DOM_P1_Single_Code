#!/usr/bin/env python3
import argparse
from pathlib import Path
import math

import pandas as pd
import matplotlib.pyplot as plt


CASE_NAME = "DOMP1ComparisionColdScatteringParallelPlatePhysicalSeperation"
DOM_ORDINATE = 8
FIGSIZE = (10.0, 6.0)
DPI = 300


PAPER_X = [
    0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
    0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
]
PAPER_G_DOM = [
    8160.0, 5990.0, 4430.0, 3440.0, 2760.0, 2160.0, 1710.0, 1480.0,
    1310.0, 1220.0, 1200.0, 1220.0, 1310.0, 1480.0, 1710.0, 2160.0,
    2760.0, 3440.0, 4430.0, 5990.0, 8160.0,
]

PAPER_Q_DOM = [
    3020.0, 2440.0, 1940.0, 1530.0, 1190.0, 900.0, 660.0, 460.0,
    290.0, 140.0, 0.0, -140.0, -290.0, -460.0, -660.0, -900.0,
    -1190.0, -1530.0, -1940.0, -2440.0, -3020.0,
]


PAPER_SR_DOM = [-2.0 * value for value in PAPER_G_DOM]



def find_case(root):
    root = Path(root).expanduser().resolve()
    candidates = [
        root,
        root / CASE_NAME,
        root / "Output" / CASE_NAME,
    ]

    for candidate in candidates:
        if (
            candidate.is_dir()
            and candidate.name == CASE_NAME
            and (candidate / "P1MeshConvergenceTest" / "mesh_convergence.csv").is_file()
        ):
            return candidate

    for candidate in root.rglob(CASE_NAME):
        if (
            candidate.is_dir()
            and "__MACOSX" not in candidate.parts
            and (candidate / "P1MeshConvergenceTest" / "mesh_convergence.csv").is_file()
        ):
            return candidate

    raise FileNotFoundError(
        "Could not find the generated cold-scattering output directory below: "
        + str(root)
    )


def read_profile(path, required_columns):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError("Missing required profile file: " + str(path))

    lines = path.read_text(encoding="utf-8").splitlines()
    header_index = None

    for index, line in enumerate(lines):
        columns = [item.strip() for item in line.split(",")]
        if all(column in columns for column in required_columns):
            header_index = index
            break

    if header_index is None:
        raise ValueError(
            "Could not find columns "
            + str(required_columns)
            + " in "
            + str(path)
        )

    frame = pd.read_csv(path, skiprows=header_index)
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        raise ValueError("Missing columns " + str(missing) + " in " + str(path))

    return frame.dropna(subset=required_columns).copy()


def require_columns(frame, columns, path):
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError("Missing columns " + str(missing) + " in " + str(path))


def save_figure(fig, output_base):
    fig.tight_layout()
    fig.savefig(str(output_base) + ".pdf", bbox_inches="tight")
    fig.savefig(str(output_base) + ".png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def escape_latex(value):
    text = str(value)
    replacements = [
        ("\\", r"\textbackslash{}"),
        ("&", r"\&"),
        ("%", r"\%"),
        ("$", r"\$"),
        ("#", r"\#"),
        ("_", r"\_"),
        ("{", r"\{"),
        ("}", r"\}"),
    ]
    for old, new in replacements:
        text = text.replace(old, new)
    return text


def format_value(value):
    if pd.isna(value):
        return "--"
    if isinstance(value, float):
        return "{:.6g}".format(value)
    return escape_latex(value)


def write_latex_table(frame, path, caption, label):
    columns = list(frame.columns)
    alignment = "l" + "r" * (len(columns) - 1)

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{" + caption + "}",
        r"\label{" + label + "}",
        r"\begin{tabular}{" + alignment + "}",
        r"\hline",
        " & ".join(escape_latex(column) for column in columns) + r" \\",
        r"\hline",
    ]

    for _, row in frame.iterrows():
        lines.append(
            " & ".join(format_value(row[column]) for column in columns) + r" \\"
        )

    lines.extend([r"\hline", r"\end{tabular}", r"\end{table}"])
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_mesh_table(case_dir, output_dir):
    columns = [
        "number_of_cells",
        "cell_size",
        "optical_cell_thickness",
        "successive_grid_relative_L2_change",
        "observed_order",
        "runtime",
    ]

    p1_path = case_dir / "P1MeshConvergenceTest" / "mesh_convergence.csv"
    dom_path = (
        case_dir
        / "DOMMeshConvergenceTest"
        / ("ordinate" + str(DOM_ORDINATE))
        / "mesh_convergence.csv"
    )

    p1 = pd.read_csv(p1_path)
    dom = pd.read_csv(dom_path)

    require_columns(p1, columns, p1_path)
    require_columns(dom, ["number_of_ordinates"] + columns, dom_path)

    p1_table = p1[columns].copy()
    p1_table.insert(0, "ordinate", pd.NA)
    p1_table.insert(0, "model", "P1")

    dom_table = dom[["number_of_ordinates"] + columns].copy()
    dom_table = dom_table.rename(columns={"number_of_ordinates": "ordinate"})
    dom_table.insert(0, "model", "DOM")

    table = pd.concat([p1_table, dom_table], ignore_index=True)
    table = table.rename(
        columns={
            "number_of_cells": "cells",
            "cell_size": "cell_size_m",
            "successive_grid_relative_L2_change": "successive_relative_L2_change_G",
            "runtime": "runtime_s",
        }
    )

    table.to_csv(output_dir / "cold_scattering_mesh_convergence_table.csv", index=False)
    write_latex_table(
        table,
        output_dir / "cold_scattering_mesh_convergence_table.tex",
        "Mesh convergence for the cold isotropically scattering parallel-plate case.",
        "tab:cold_scattering_mesh_convergence",
    )


def plot_profiles(case_dir, output_dir):
    p1_dir = case_dir / "P1RadiativeEnergyAndFluxProfileTest" / "300x1x1_1"
    dom_dir = (
        case_dir
        / "DOMRadiativeEnergyAndFluxProfileTest"
        / ("ordinate" + str(DOM_ORDINATE))
        / "300x1x1_1"
    )

    specifications = [
        (
            "centerlineG.csv",
            "G_P1",
            "G_DOM",
            r"Irradiation, $G$ (W m$^{-2}$)",
            "cold_scattering_centerline_irradiation",
        ),
        (
            "centerlineFlux.csv",
            "radiative_flux_x_P1",
            "radiative_flux_x_DOM",
            r"Radiative heat flux, $q_x$ (W m$^{-2}$)",
            "cold_scattering_centerline_heat_flux",
        ),
        (
            "centerlineRadiativeSourceTerm.csv",
            "radiative_source_term_P1",
            "radiative_source_term_DOM",
            r"Radiative source term, $S_r$ (W m$^{-3}$)",
            "cold_scattering_centerline_source",
        ),
    ]

    for filename, p1_column, dom_column, ylabel, output_name in specifications:
        p1 = read_profile(p1_dir / filename, ["x_over_L", p1_column])
        dom = read_profile(dom_dir / filename, ["x_over_L", dom_column])

        p1 = p1.sort_values("x_over_L")
        dom = dom.sort_values("x_over_L")

        fig, ax = plt.subplots(figsize=FIGSIZE)
        ax.plot(
            p1["x_over_L"],
            p1[p1_column],
            linestyle="--",
            linewidth=2.2,
            label="P1",
        )
        ax.plot(
            dom["x_over_L"],
            dom[dom_column],
            linestyle="-",
            linewidth=2.2,
            label=r"DOM ($S_8$)",
        )

        if filename == "centerlineG.csv":
            paper_values = PAPER_G_DOM
        elif filename == "centerlineFlux.csv":
            paper_values = PAPER_Q_DOM
        else:
            paper_values = PAPER_SR_DOM

        ax.scatter(
            PAPER_X,
            paper_values,
            marker="o",
            s=28,
            facecolors="none",
            edgecolors="black",
            linewidths=1.0,
            label="Sarr and Groth DOM",
            zorder=4,
        )

        ax.set_xlabel(r"Dimensionless position, $x/L$")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0.0, 1.0)
        ax.grid(True, which="major", alpha=0.25)
        ax.legend(frameon=False)

        save_figure(fig, output_dir / output_name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "root",
        nargs="?",
        default="Output",
        help="Output directory, case directory, or project root.",
    )
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()

    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.labelsize": 14,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 13,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )

    case_dir = find_case(args.root)
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else case_dir / "PlotsAndTables"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    make_mesh_table(case_dir, output_dir)
    plot_profiles(case_dir, output_dir)

    print("Cold-scattering plots and table written to:")
    print(output_dir)


if __name__ == "__main__":
    main()
