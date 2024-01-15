from itertools import chain
from pathlib import Path

from matplotlib import pyplot as plt
from pyryotype import plot_ideogram
from pyryotype.ideogram import GENOME, Detail, Orientation

OUT_DIR = Path(__file__).parent.parent / "example_outputs"


def test_simple_vertical_chr1():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(3, 25),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr1", left_margin=0, y_label="", vertical=Orientation.VERTICAL)

    fig.savefig(OUT_DIR / "testing_vert.png", bbox_inches="tight")


def test_23_vertical_chm13():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=24, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        _ax = plot_ideogram(
            ax, target=f"chr{i}", y_label="Chr. {i}", left_margin=0, vertical=Orientation.VERTICAL, genome=genome
        )

    fig.savefig(OUT_DIR / "testing_vert_23.png", bbox_inches="tight")


def test_23_vertical_hg38():
    genome = GENOME.HG38
    fig, axes = plt.subplots(ncols=22, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23)), strict=True):
        _ax = plot_ideogram(
            ax, target=f"chr{i}", y_label="Chr. {i}", left_margin=0, vertical=Orientation.VERTICAL, genome=genome
        )

    fig.savefig(OUT_DIR / "testing_vert_23_hg38.png", bbox_inches="tight")


def test_23_vertical_chm13_bare():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=24, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            y_label="Chr. {i}",
            left_margin=0,
            vertical=Orientation.VERTICAL,
            cytobands=Detail.BARE,
            genome=genome,
        )

    fig.savefig(OUT_DIR / "testing_vert_bare_23.png", bbox_inches="tight")


def test_23_horizontal_chm13_bare():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=1, nrows=24, figsize=(25, 15), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            y_label="Chr. {i}",
            left_margin=0,
            vertical=Orientation.HORIZONTAL,
            cytobands=Detail.BARE,
            genome=genome,
        )

    fig.savefig(OUT_DIR / "testing_horz_bare_23.png", bbox_inches="tight")


def test_23_vertical_chm13_regions():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=24, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    per_chr_regions = {"chr1": [(0, 1000000, "black"), (20_000_000, 25_000_000, "red")]}
    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        regions = per_chr_regions.get(f"chr{i}")
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            y_label=f"Chr. {i}",
            left_margin=0,
            vertical=Orientation.VERTICAL,
            genome=genome,
            regions=regions,
        )

    fig.savefig(OUT_DIR / "testing_vert_23_regions.png", bbox_inches="tight")


def test_23_horz_chm13_regions():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=1, nrows=24, figsize=(15, 25), facecolor="white", sharey=True)

    per_chr_regions = {"chr1": [(0, 1000000, "black"), (20_000_000, 25_000_000, "red")]}
    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        regions = per_chr_regions.get(f"chr{i}", None)
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            y_label=f"Chr. {i}",
            left_margin=0,
            height=0.99,
            vertical=Orientation.HORIZONTAL,
            genome=genome,
            regions=regions,
        )

    fig.savefig(OUT_DIR / "testing_horz_23_regions.png", bbox_inches="tight")
