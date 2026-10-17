#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Any


DEFAULT_CASE_DIRECTORY = "DOMP1FiniteParallelPlateMemoryAccuracy"
DEFAULT_DOMP1_FILENAME = "unified_memory_accuracy_domp1.csv"
DEFAULT_DTM_FILENAME = "unified_memory_accuracy_dtm.csv"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot relative plate-A heat-rate error against peak host memory "
            "for P1, DOM and DTM, with each point labelled by its numerical "
            "configuration."
        )
    )

    parser.add_argument(
        "--input-directory",
        type=Path,
        default=None,
        help=(
            "Directory containing the DOM/P1 and DTM CSV files. By default, "
            "the script uses "
            "<project-root>/output/DOMP1FiniteParallelPlateMemoryAccuracy."
        ),
    )

    parser.add_argument(
        "--figure-directory",
        type=Path,
        default=None,
        help=(
            "Directory in which figures are written. By default, figures are "
            "written to a 'figures' directory beside the CSV files."
        ),
    )

    return parser.parse_args()


def project_root_from_script() -> Path:
    script_path = Path(__file__).resolve()

    try:
        return script_path.parents[2]
    except IndexError as exc:
        raise RuntimeError(
            "Could not determine the project root from the script location."
        ) from exc


def find_csv(directory: Path, preferred_name: str, stem_prefix: str) -> Path:
    preferred_path = directory / preferred_name

    if preferred_path.is_file():
        return preferred_path

    candidates = sorted(directory.glob(f"{stem_prefix}*.csv"))

    if len(candidates) == 1:
        return candidates[0]

    if not candidates:
        raise FileNotFoundError(
            f"No CSV matching '{stem_prefix}*.csv' was found in {directory}."
        )

    candidate_names = ", ".join(path.name for path in candidates)

    raise RuntimeError(
        f"Several CSV files match '{stem_prefix}*.csv' in {directory}: "
        f"{candidate_names}. Rename the required file to '{preferred_name}'."
    )


def safe_float(value: str, column_name: str, csv_path: Path) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid value '{value}' in column '{column_name}' of {csv_path}."
        ) from exc


def safe_int(value: str, column_name: str, csv_path: Path) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid value '{value}' in column '{column_name}' of {csv_path}."
        ) from exc


def read_results(csv_path: Path) -> list[dict[str, Any]]:
    required_columns = {
        "process_name",
        "configuration",
        "resolution",
        "number_of_cells",
        "ordinate_order",
        "number_of_ordinates",
        "peak_RSS_MiB",
        "relative_error_plate_A",
    }

    results: list[dict[str, Any]] = []

    with csv_path.open("r", newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file)

        if reader.fieldnames is None:
            raise ValueError(f"{csv_path} does not contain a CSV header.")

        missing_columns = required_columns.difference(reader.fieldnames)

        if missing_columns:
            missing_text = ", ".join(sorted(missing_columns))
            raise ValueError(
                f"{csv_path} is missing required columns: {missing_text}"
            )

        for row_number, row in enumerate(reader, start=2):
            process_name = row["process_name"].strip().upper()
            configuration = row["configuration"].strip()

            peak_memory = safe_float(
                row["peak_RSS_MiB"],
                "peak_RSS_MiB",
                csv_path,
            )

            relative_error = safe_float(
                row["relative_error_plate_A"],
                "relative_error_plate_A",
                csv_path,
            )

            resolution = safe_int(
                row["resolution"],
                "resolution",
                csv_path,
            )

            number_of_cells = safe_int(
                row["number_of_cells"],
                "number_of_cells",
                csv_path,
            )

            ordinate_order = safe_int(
                row["ordinate_order"],
                "ordinate_order",
                csv_path,
            )

            number_of_ordinates = safe_int(
                row["number_of_ordinates"],
                "number_of_ordinates",
                csv_path,
            )

            if peak_memory <= 0.0 or relative_error <= 0.0:
                print(
                    f"Skipping row {row_number} of {csv_path.name}: "
                    "memory and error must both be positive for logarithmic "
                    "axes."
                )
                continue

            results.append(
                {
                    "process_name": process_name,
                    "configuration": configuration,
                    "resolution": resolution,
                    "number_of_cells": number_of_cells,
                    "ordinate_order": ordinate_order,
                    "number_of_ordinates": number_of_ordinates,
                    "peak_RSS_MiB": peak_memory,
                    "relative_error_plate_A": relative_error,
                }
            )

    return results


def format_mesh(configuration: str) -> str:
    mesh_match = re.search(r"(\d+)x(\d+)x(\d+)", configuration)

    if mesh_match is None:
        return configuration

    return (
        f"{mesh_match.group(1)}"
        f"\u00d7{mesh_match.group(2)}"
        f"\u00d7{mesh_match.group(3)}"
    )


def make_point_label(result: dict[str, Any]) -> str:
    method = result["process_name"]
    configuration = result["configuration"]

    if method == "P1":
        return format_mesh(configuration)

    if method == "DOM":
        mesh = format_mesh(configuration)
        ordinate_order = result["ordinate_order"]
        return rf"{mesh}, $S_{{{ordinate_order}}}$"

    if method == "DTM":
        surface_match = re.search(r"(\d+)x(\d+)_gap", configuration)
        theta_match = re.search(r"_theta(\d+)", configuration)
        psi_match = re.search(r"_psi(\d+)", configuration)

        if surface_match is not None:
            surface_mesh = (
                f"{surface_match.group(1)}"
                f"\u00d7{surface_match.group(2)}"
            )
        else:
            surface_mesh = configuration

        if theta_match is not None and psi_match is not None:
            return (
                rf"{surface_mesh}, "
                rf"$N_\theta={theta_match.group(1)}$, "
                rf"$N_\psi={psi_match.group(1)}$"
            )

        return surface_mesh

    return configuration


def annotation_offset(
    method: str,
    method_index: int,
) -> tuple[int, int]:
    offset_patterns = {
        "P1": [
            (7, 7),
            (7, -13),
            (7, 7),
            (7, -13),
        ],
        "DOM": [
            (7, 7),
            (7, -13),
            (-7, 7),
            (-7, -13),
            (9, 11),
            (9, -17),
            (-9, 11),
            (-9, -17),
        ],
        "DTM": [
            (7, 7),
            (7, -13),
            (7, 7),
            (7, -13),
            (7, 7),
        ],
    }

    pattern = offset_patterns.get(method, [(7, 7)])
    return pattern[method_index % len(pattern)]


def configure_matplotlib() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.grid": True,
            "grid.alpha": 0.3,
        }
    )


def create_annotated_scatter_plot(
    results: list[dict[str, Any]],
    output_base: Path,
) -> None:
    import matplotlib.pyplot as plt

    marker_by_method = {
        "P1": "o",
        "DOM": "s",
        "DTM": "^",
    }

    display_name_by_method = {
        "P1": r"$P_1$",
        "DOM": "DOM",
        "DTM": "DTM",
    }

    figure, axis = plt.subplots(figsize=(9.2, 6.2))

    for method in ["P1", "DOM", "DTM"]:
        method_results = [
            result
            for result in results
            if result["process_name"] == method
        ]

        method_results.sort(
            key=lambda result: (
                result["peak_RSS_MiB"],
                result["relative_error_plate_A"],
            )
        )

        if not method_results:
            continue

        x_values = [
            result["peak_RSS_MiB"]
            for result in method_results
        ]

        y_values = [
            result["relative_error_plate_A"]
            for result in method_results
        ]

        axis.scatter(
            x_values,
            y_values,
            marker=marker_by_method[method],
            s=58,
            label=display_name_by_method[method],
            zorder=3,
        )

        for point_index, result in enumerate(method_results):
            x_offset, y_offset = annotation_offset(method, point_index)

            horizontal_alignment = (
                "left" if x_offset >= 0 else "right"
            )

            vertical_alignment = (
                "bottom" if y_offset >= 0 else "top"
            )

            axis.annotate(
                make_point_label(result),
                xy=(
                    result["peak_RSS_MiB"],
                    result["relative_error_plate_A"],
                ),
                xytext=(x_offset, y_offset),
                textcoords="offset points",
                fontsize=8,
                ha=horizontal_alignment,
                va=vertical_alignment,
                annotation_clip=False,
            )

    axis.set_xscale("log")
    axis.set_yscale("log")

    axis.set_xlabel(
        r"Peak resident-set memory, $M_{\mathrm{peak}}$ (MiB)"
    )

    axis.set_ylabel(
        r"Relative plate-A heat-rate error, $E_{Q_A}$"
    )

    axis.set_title(
        "Plate-A heat-rate accuracy versus peak host memory"
    )

    axis.legend(
        title="Radiation model",
        loc="best",
    )

    axis.margins(x=0.16, y=0.18)

    figure.tight_layout()

    output_base.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_base.with_suffix(".pdf"))
    figure.savefig(output_base.with_suffix(".png"))
    plt.close(figure)


def main() -> None:
    arguments = parse_arguments()
    project_root = project_root_from_script()

    input_directory = (
        arguments.input_directory.resolve()
        if arguments.input_directory is not None
        else (
            project_root
            / "output"
            / DEFAULT_CASE_DIRECTORY
        )
    )

    figure_directory = (
        arguments.figure_directory.resolve()
        if arguments.figure_directory is not None
        else input_directory / "figures"
    )

    domp1_csv = find_csv(
        input_directory,
        DEFAULT_DOMP1_FILENAME,
        "unified_memory_accuracy_domp1",
    )

    dtm_csv = find_csv(
        input_directory,
        DEFAULT_DTM_FILENAME,
        "unified_memory_accuracy_dtm",
    )

    results = read_results(domp1_csv)
    results.extend(read_results(dtm_csv))

    methods_present = {
        result["process_name"]
        for result in results
    }

    expected_methods = {"P1", "DOM", "DTM"}
    missing_methods = expected_methods.difference(methods_present)

    if missing_methods:
        missing_text = ", ".join(sorted(missing_methods))
        raise ValueError(
            f"No valid rows were found for: {missing_text}"
        )

    configure_matplotlib()

    output_base = (
        figure_directory
        / "error_vs_peak_memory_annotated"
    )

    create_annotated_scatter_plot(
        results,
        output_base,
    )

    print(f"Read DOM/P1 results from: {domp1_csv}")
    print(f"Read DTM results from:    {dtm_csv}")
    print(f"PDF figure written to:    {output_base.with_suffix('.pdf')}")
    print(f"PNG figure written to:    {output_base.with_suffix('.png')}")


if __name__ == "__main__":
    main()
