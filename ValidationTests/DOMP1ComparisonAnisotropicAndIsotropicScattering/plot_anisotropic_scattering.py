#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


CASE_NAME = "DOMP1ComparisonAnisotropicAndIsotropicScattering"
DOM_ORDINATE = 8
ANISOTROPY_FACTORS = [1.0, 0.5, 0.0, -0.5, -1.0]
PROFILE_FACTORS = [1.0, 0.0, -1.0]
PROFILE_BETAS = [0.1, 1.0, 3.0, 5.0, 10.0]
FIGSIZE = (10.5, 6.5)
DPI = 300


PAPER_DOM_TABLE = {
    (0.0, 0.1): {6: 0.9142, 12: 0.9141, 24: 0.9150, 36: 0.9153, 48: 0.9155},
    (0.0, 1.0): {6: 0.5471, 12: 0.5517, 24: 0.5530, 36: 0.5532, 48: 0.5533},
    (0.0, 10.0): {6: 0.1155, 12: 0.1164, 24: 0.1166, 36: 0.1167, 48: 0.1167},
    (1.0, 0.1): {6: 0.9357, 12: 0.9355, 24: 0.9364, 36: 0.9368, 48: 0.9370},
    (1.0, 1.0): {6: 0.6338, 12: 0.6400, 24: 0.6417, 36: 0.6419, 48: 0.6420},
    (1.0, 10.0): {6: 0.1624, 12: 0.1642, 24: 0.1647, 36: 0.1648, 48: 0.1648},
}

PAPER_EXACT = {
    (0.0, 0.1): 0.9157,
    (0.0, 1.0): 0.5532,
    (0.0, 10.0): 0.1167,
    (1.0, 0.1): 0.9371,
    (1.0, 1.0): 0.6423,
    (1.0, 10.0): 0.1649,
}

PAPER_PROFILE_X = [0.0, 0.25, 0.5, 0.75, 1.0]
PAPER_EMISSIVE = {
    0.0: {
        0.1: [0.43, 0.465, 0.50, 0.535, 0.57],
        1.0: [0.32, 0.405, 0.50, 0.595, 0.68],
        3.0: [0.23, 0.355, 0.50, 0.645, 0.77],
        5.0: [0.17, 0.325, 0.50, 0.675, 0.83],
        10.0: [0.10, 0.285, 0.50, 0.715, 0.90],
    },
    1.0: {
        0.1: [0.44, 0.470, 0.50, 0.530, 0.56],
        1.0: [0.34, 0.415, 0.50, 0.585, 0.66],
        3.0: [0.25, 0.365, 0.50, 0.635, 0.75],
        5.0: [0.19, 0.335, 0.50, 0.665, 0.81],
        10.0: [0.11, 0.295, 0.50, 0.705, 0.89],
    },
    -1.0: {
        0.1: [0.42, 0.460, 0.50, 0.540, 0.58],
        1.0: [0.29, 0.395, 0.50, 0.605, 0.71],
        3.0: [0.20, 0.345, 0.50, 0.655, 0.80],
        5.0: [0.14, 0.315, 0.50, 0.685, 0.86],
        10.0: [0.08, 0.275, 0.50, 0.725, 0.92],
    },
}


SINGH_REFERENCE = {
    (0.0, 0.1): 0.9157,
    (0.0, 1.0): 0.5532,
    (0.0, 10.0): 0.1167,
    (1.0, 0.1): 0.9371,
    (1.0, 1.0): 0.6423,
    (1.0, 10.0): 0.1649,
}


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
            and (candidate / "DOMAngularConvergence" / "angular_convergence.csv").is_file()
        ):
            return candidate

    for candidate in root.rglob(CASE_NAME):
        if (
            candidate.is_dir()
            and "__MACOSX" not in candidate.parts
            and (candidate / "DOMAngularConvergence" / "angular_convergence.csv").is_file()
        ):
            return candidate

    raise FileNotFoundError(
        "Could not find the generated anisotropic-scattering output directory below: "
        + str(root)
    )


def require_columns(frame, columns, path):
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError("Missing columns " + str(missing) + " in " + str(path))


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
    require_columns(frame, required_columns, path)
    return frame.dropna(subset=required_columns).copy()


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


def create_angular_table(case_dir, output_dir):
    path = case_dir / "DOMAngularConvergence" / "angular_convergence.csv"
    columns = [
        "anisotropy_factor",
        "extinction_coefficient",
        "ordinate_order",
        "dimensionless_east_wall_heat_flux",
        "runtime",
    ]

    present = pd.read_csv(path)
    require_columns(present, columns, path)
    present = present[columns].copy()

    present_rows = []
    for _, row in present.iterrows():
        key = (
            float(row["anisotropy_factor"]),
            float(row["extinction_coefficient"]),
        )
        paper_48 = PAPER_DOM_TABLE.get(key, {}).get(48, float("nan"))
        exact = PAPER_EXACT.get(key, float("nan"))
        present_rows.append(
            {
                "source": "Present DOM",
                "anisotropy_parameter": row["anisotropy_factor"],
                "beta": row["extinction_coefficient"],
                "angular_resolution": "S{}".format(int(row["ordinate_order"])),
                "dimensionless_heat_flux": row["dimensionless_east_wall_heat_flux"],
                "paper_DOM_48_rays": paper_48,
                "paper_exact_reference": exact,
                "percentage_difference_from_paper_DOM_48": (
                    100.0
                    * (row["dimensionless_east_wall_heat_flux"] - paper_48)
                    / paper_48
                    if not pd.isna(paper_48)
                    else float("nan")
                ),
                "runtime_s": row["runtime"],
            }
        )

    paper_rows = []
    for (factor, beta), ray_values in PAPER_DOM_TABLE.items():
        for rays, value in ray_values.items():
            paper_rows.append(
                {
                    "source": "Singh DOM",
                    "anisotropy_parameter": factor,
                    "beta": beta,
                    "angular_resolution": "{} rays".format(rays),
                    "dimensionless_heat_flux": value,
                    "paper_DOM_48_rays": ray_values[48],
                    "paper_exact_reference": PAPER_EXACT[(factor, beta)],
                    "percentage_difference_from_paper_DOM_48": (
                        100.0 * (value - ray_values[48]) / ray_values[48]
                    ),
                    "runtime_s": float("nan"),
                }
            )

    table = pd.DataFrame(present_rows + paper_rows)
    table = table.sort_values(
        ["anisotropy_parameter", "beta", "source", "angular_resolution"]
    )

    table.to_csv(output_dir / "dom_angular_convergence_table.csv", index=False)
    write_latex_table(
        table,
        output_dir / "dom_angular_convergence_table.tex",
        "Angular convergence of the present DOM and the published Singh and Mishra DOM results.",
        "tab:dom_scattering_angular_convergence",
    )

def plot_wall_flux(case_dir, output_dir):
    p1_path = case_dir / "P1Scattering" / "dimensionless_heat_flux.csv"
    dom_path = case_dir / "DOMScattering" / "dimensionless_heat_flux.csv"

    required = [
        "anisotropy_factor",
        "extinction_coefficient",
        "dimensionless_east_wall_heat_flux",
    ]

    p1 = pd.read_csv(p1_path)
    dom = pd.read_csv(dom_path)

    require_columns(p1, required, p1_path)
    require_columns(dom, required + ["ordinate_order"], dom_path)

    dom = dom[dom["ordinate_order"] == DOM_ORDINATE].copy()

    fig, ax = plt.subplots(figsize=FIGSIZE)
    colors = plt.get_cmap("tab10").colors

    for index, factor in enumerate(ANISOTROPY_FACTORS):
        color = colors[index]
        p1_case = p1[p1["anisotropy_factor"] == factor].sort_values(
            "extinction_coefficient"
        )
        dom_case = dom[dom["anisotropy_factor"] == factor].sort_values(
            "extinction_coefficient"
        )

        if p1_case.empty or dom_case.empty:
            raise ValueError("Missing wall-flux data for anisotropy factor " + str(factor))

        ax.plot(
            p1_case["extinction_coefficient"],
            p1_case["dimensionless_east_wall_heat_flux"],
            linestyle="--",
            linewidth=2.1,
            color=color,
        )
        ax.plot(
            dom_case["extinction_coefficient"],
            dom_case["dimensionless_east_wall_heat_flux"],
            linestyle="-",
            linewidth=2.1,
            color=color,
        )

        paper_points = [
            (beta, PAPER_DOM_TABLE[(factor, beta)][48])
            for beta in [0.1, 1.0, 10.0]
            if (factor, beta) in PAPER_DOM_TABLE
        ]
        if paper_points:
            ax.scatter(
                [point[0] for point in paper_points],
                [point[1] for point in paper_points],
                marker="o",
                s=42,
                facecolors="none",
                edgecolors=color,
                linewidths=1.4,
                zorder=5,
            )

    factor_handles = [
        Line2D(
            [0],
            [0],
            color=colors[index],
            linewidth=2.1,
            label=r"$g={:g}$".format(factor),
        )
        for index, factor in enumerate(ANISOTROPY_FACTORS)
    ]
    model_handles = [
        Line2D([0], [0], color="black", linestyle="--", linewidth=2.1, label="P1"),
        Line2D([0], [0], color="black", linestyle="-", linewidth=2.1, label=r"DOM ($S_8$)"),
        Line2D(
            [0],
            [0],
            color="black",
            marker="o",
            markerfacecolor="none",
            linestyle="None",
            label="Singh DOM (48 rays)",
        ),
    ]

    first_legend = ax.legend(
        handles=factor_handles,
        title="Anisotropy parameter",
        frameon=False,
        loc="upper right",
    )
    ax.add_artist(first_legend)
    ax.legend(handles=model_handles, title="Model", frameon=False, loc="lower left")

    ax.set_xlabel(r"Extinction coefficient, $\beta$ (m$^{-1}$)")
    ax.set_ylabel(r"Dimensionless wall heat flux, $\Psi$")
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, which="major", alpha=0.25)

    save_figure(fig, output_dir / "singh_dimensionless_wall_heat_flux")


def factor_directory(value):
    return "anisotropyFactor{:.6f}".format(value)


def extinction_directory(value):
    return "extinctionCoefficient{0:.6f}_sigmaS{0:.6f}".format(value)


def load_emissive_profile(case_dir, model, factor, beta):
    if model == "P1":
        path = (
            case_dir
            / "P1Scattering"
            / factor_directory(factor)
            / extinction_directory(beta)
            / "300x1x1"
            / "dimensionless_emissive_power.csv"
        )
        return read_profile(
            path,
            ["x_over_X", "dimensionless_emissive_power_P1"],
        )

    path = (
        case_dir
        / "DOMScattering"
        / factor_directory(factor)
        / extinction_directory(beta)
        / ("ordinate" + str(DOM_ORDINATE))
        / "dimensionless_emissive_power.csv"
    )
    return read_profile(
        path,
        ["x_over_X", "dimensionless_emissive_power_DOM"],
    )


def plot_emissive_power(case_dir, output_dir, factor):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    colors = plt.get_cmap("tab10").colors

    for index, beta in enumerate(PROFILE_BETAS):
        color = colors[index]
        p1 = load_emissive_profile(case_dir, "P1", factor, beta).sort_values("x_over_X")
        dom = load_emissive_profile(case_dir, "DOM", factor, beta).sort_values("x_over_X")

        ax.plot(
            p1["x_over_X"],
            p1["dimensionless_emissive_power_P1"],
            linestyle="--",
            linewidth=2.1,
            color=color,
        )
        ax.plot(
            dom["x_over_X"],
            dom["dimensionless_emissive_power_DOM"],
            linestyle="-",
            linewidth=2.1,
            color=color,
        )

        paper_values = PAPER_EMISSIVE[factor][beta]
        ax.scatter(
            PAPER_PROFILE_X,
            paper_values,
            marker="o",
            s=28,
            facecolors="none",
            edgecolors=color,
            linewidths=1.0,
            zorder=5,
        )

    beta_handles = [
        Line2D(
            [0],
            [0],
            color=colors[index],
            linewidth=2.1,
            label=r"$\beta={:g}$".format(beta),
        )
        for index, beta in enumerate(PROFILE_BETAS)
    ]
    model_handles = [
        Line2D([0], [0], color="black", linestyle="--", linewidth=2.1, label="P1"),
        Line2D([0], [0], color="black", linestyle="-", linewidth=2.1, label=r"DOM ($S_8$)"),
        Line2D(
            [0],
            [0],
            color="black",
            marker="o",
            markerfacecolor="none",
            linestyle="None",
            label="Singh et al.",
        ),
    ]

    first_legend = ax.legend(
        handles=beta_handles,
        title=r"Extinction coefficient (m$^{-1}$)",
        frameon=False,
        loc="upper left",
    )
    ax.add_artist(first_legend)
    ax.legend(handles=model_handles, title="Model", frameon=False, loc="lower right")

    ax.set_xlabel(r"Dimensionless position, $x/L$")
    ax.set_ylabel(r"Dimensionless emissive power, $\Phi=(T/T_E)^4$")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, which="major", alpha=0.25)

    if factor > 0:
        name = "forward"
    elif factor < 0:
        name = "backward"
    else:
        name = "isotropic"

    save_figure(
        fig,
        output_dir / ("singh_emissive_power_" + name),
    )


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
            "legend.fontsize": 12,
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

    create_angular_table(case_dir, output_dir)
    plot_wall_flux(case_dir, output_dir)

    for factor in PROFILE_FACTORS:
        plot_emissive_power(case_dir, output_dir, factor)

    print("Anisotropic-scattering plots and table written to:")
    print(output_dir)


if __name__ == "__main__":
    main()
