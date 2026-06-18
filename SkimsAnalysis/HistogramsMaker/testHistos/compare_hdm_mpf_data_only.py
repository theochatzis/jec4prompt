#!/usr/bin/env python3
"""
Compare MPF/DB with fixed-Rn/Ru HDM in the barrel using data only.

Expected input ROOT objects:
  - MPF_2D
  - DB_2D
  - HDM_r0_2D
  - HDM_r1_2D
  - HDM_rn_2D
  - HDM_ru_2D

Default input:
  test_closure_data.root

Example:
  python3 compare_hdm_mpf_data_only.py \
    --data test_closure_data.root \
    --eta-min -1.305 \
    --eta-max 1.305 \
    --rn 1.0 \
    --ru 0.92 \
    --out hdm_mpf_data_only.root \
    --plot-dir hdm_mpf_data_only_plots
"""

import argparse
import os
import math
import ROOT


def setup_root_style():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gROOT.ForceStyle()


def ensure_dir(path):
    if path and not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)


def require_object(root_file, name):
    obj = root_file.Get(name)
    if not obj:
        raise RuntimeError(f"Could not find object '{name}' in {root_file.GetName()}")
    return obj


def format_x_axis(obj):
    obj.GetXaxis().SetNoExponent(True)
    obj.GetXaxis().SetMoreLogLabels(True)


def remove_stat_boxes(canvas):
    """Remove TPaveStats objects that ROOT may attach to primitives."""
    canvas.Update()
    prims = canvas.GetListOfPrimitives()
    to_remove = []
    for prim in prims:
        if prim.InheritsFrom("TPaveStats"):
            to_remove.append(prim)
    for prim in to_remove:
        prims.Remove(prim)
    canvas.Modified()
    canvas.Update()


def project_barrel_profile_y(profile2d, name, eta_min, eta_max):
    """
    Follow the L2L3Res.C style:
      - choose eta bins on X axis
      - ProfileY over that eta range
      - ProjectionX to TH1D
    """
    xaxis = profile2d.GetXaxis()

    bin_start = xaxis.FindBin(eta_min + 1e-5)
    bin_end = xaxis.FindBin(eta_max - 1e-5)

    prof = profile2d.ProfileY(f"{name}_prof", bin_start, bin_end)
    hist = prof.ProjectionX(name)

    hist.SetDirectory(0)
    hist.SetStats(False)
    format_x_axis(hist)

    return hist


def hdm_mpf_value(r0, rn, ru, Rn, Ru):
    denom = 1.0 - rn / Rn - ru / Ru
    if abs(denom) < 1e-12:
        return 0.0
    return (r0 - rn - ru) / denom


def hdm_db_value(r1, rn, ru, Rn, Ru):
    denom = 1.0 - rn / Rn - ru / Ru
    if abs(denom) < 1e-12:
        return 0.0
    return r1 / denom


def hdm_mpf_error(r0, rn, ru, er0, ern, eru, Rn, Ru):
    """
    Error propagation for:
      H = (r0 - rn - ru) / (1 - rn/Rn - ru/Ru)

    Correlations are ignored here.
    """
    D = 1.0 - rn / Rn - ru / Ru
    if abs(D) < 1e-12:
        return 0.0

    N = r0 - rn - ru

    d_dr0 = 1.0 / D
    d_drn = (-D + N / Rn) / (D * D)
    d_dru = (-D + N / Ru) / (D * D)

    err2 = (
        (d_dr0 * er0) ** 2
        + (d_drn * ern) ** 2
        + (d_dru * eru) ** 2
    )

    return math.sqrt(max(err2, 0.0))


def hdm_db_error(r1, rn, ru, er1, ern, eru, Rn, Ru):
    """
    Error propagation for:
      H = r1 / (1 - rn/Rn - ru/Ru)

    Correlations are ignored here.
    """
    D = 1.0 - rn / Rn - ru / Ru
    if abs(D) < 1e-12:
        return 0.0

    d_dr1 = 1.0 / D
    d_drn = r1 / (Rn * D * D)
    d_dru = r1 / (Ru * D * D)

    err2 = (
        (d_dr1 * er1) ** 2
        + (d_drn * ern) ** 2
        + (d_dru * eru) ** 2
    )

    return math.sqrt(max(err2, 0.0))


def build_hdm_hist(name, h_r0, h_r1, h_rn, h_ru, Rn, Ru, branch):
    """
    Build fixed-Rn/Ru HDM response from component histograms.

    branch = "mpf":
      H = (r0 - rn - ru) / (1 - rn/Rn - ru/Ru)

    branch = "db":
      H = r1 / (1 - rn/Rn - ru/Ru)
    """
    h = h_r0.Clone(name)
    h.SetDirectory(0)
    h.Reset()
    h.SetStats(False)
    format_x_axis(h)

    for b in range(1, h.GetNbinsX() + 1):
        r0 = h_r0.GetBinContent(b)
        r1 = h_r1.GetBinContent(b)
        rn = h_rn.GetBinContent(b)
        ru = h_ru.GetBinContent(b)

        er0 = h_r0.GetBinError(b)
        er1 = h_r1.GetBinError(b)
        ern = h_rn.GetBinError(b)
        eru = h_ru.GetBinError(b)

        if branch == "mpf":
            value = hdm_mpf_value(r0, rn, ru, Rn, Ru)
            error = hdm_mpf_error(r0, rn, ru, er0, ern, eru, Rn, Ru)
        elif branch == "db":
            value = hdm_db_value(r1, rn, ru, Rn, Ru)
            error = hdm_db_error(r1, rn, ru, er1, ern, eru, Rn, Ru)
        else:
            raise ValueError(f"Unknown HDM branch: {branch}")

        h.SetBinContent(b, value)
        h.SetBinError(b, error)

    return h


def clone_ratio(name, numerator, denominator):
    h = numerator.Clone(name)
    h.SetDirectory(0)
    h.SetStats(False)
    format_x_axis(h)
    h.Divide(denominator)
    return h


def set_hist_style(hist, color, marker_style, line_style=1):
    hist.SetStats(False)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(marker_style)
    hist.SetLineStyle(line_style)
    hist.SetLineWidth(2)
    hist.SetMarkerSize(0.9)
    format_x_axis(hist)


def hist_y_range(hists, default=(0.5, 1.5), pad=0.15):
    ymin = float("inf")
    ymax = -float("inf")

    for h in hists:
        if not h:
            continue
        for b in range(1, h.GetNbinsX() + 1):
            y = h.GetBinContent(b)
            e = h.GetBinError(b)
            if y == 0 and e == 0:
                continue
            ymin = min(ymin, y - e)
            ymax = max(ymax, y + e)

    if not math.isfinite(ymin) or not math.isfinite(ymax) or ymin == ymax:
        return default

    dy = ymax - ymin
    return ymin*0.7 - pad * dy, ymax*1.3 + pad * dy


def draw_overlay(
    output_path,
    hists,
    labels,
    title,
    ytitle,
    logx=True,
    y_range=None,
    draw_line_at_one=True,
):
    colors = [ROOT.kBlack, ROOT.kGray + 2, ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1]
    markers = [20, 24, 21, 22, 23, 25]

    non_null_hists = [h for h in hists if h]
    if not non_null_hists:
        return

    for i, h in enumerate(non_null_hists):
        set_hist_style(h, colors[i % len(colors)], markers[i % len(markers)])

    xmin = min(h.GetXaxis().GetXmin() for h in non_null_hists)
    xmax = max(h.GetXaxis().GetXmax() for h in non_null_hists)

    if y_range is None:
        y_range = hist_y_range(non_null_hists)

    canvas = ROOT.TCanvas(os.path.basename(output_path), "", 900, 700)
    canvas.SetRightMargin(0.04)
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.12)

    if logx:
        canvas.SetLogx(True)

    frame = canvas.DrawFrame(xmin, y_range[0], xmax, y_range[1], f"{title};p^{{tag}}_{{T}} (GeV);{ytitle}")
    frame.SetStats(False)
    format_x_axis(frame)
    frame.GetYaxis().SetTitleOffset(1.25)

    if draw_line_at_one:
        line = ROOT.TLine(xmin, 1.0, xmax, 1.0)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(ROOT.kDashed)
        line.Draw("same")
    else:
        line = None

    leg = ROOT.TLegend(0.58, 0.68, 0.92, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)

    for h, label in zip(non_null_hists, labels):
        h.Draw("E1 SAME")
        leg.AddEntry(h, label, "pe")

    leg.Draw()

    remove_stat_boxes(canvas)
    canvas.SaveAs(output_path)

    # Keep objects alive until after SaveAs
    if line:
        line.Draw()


def write_all(output_file, objects):
    output_file.cd()
    for obj in objects:
        if obj:
            obj.Write(obj.GetName(), ROOT.TObject.kOverwrite)


def main():
    parser = argparse.ArgumentParser(description="Data-only MPF vs HDM barrel comparison.")
    parser.add_argument("--data", default="test_closure_data.root", help="Input data ROOT file.")
    parser.add_argument("--out", default="hdm_mpf_data_only.root", help="Output ROOT file.")
    parser.add_argument("--plot-dir", default="hdm_mpf_data_only_plots", help="Directory for PNG plots.")
    parser.add_argument("--eta-min", type=float, default=-1.305, help="Minimum probe eta for barrel projection.")
    parser.add_argument("--eta-max", type=float, default=1.305, help="Maximum probe eta for barrel projection.")
    parser.add_argument("--rn", type=float, default=1.0, help="Fixed Rn value.")
    parser.add_argument("--ru", type=float, default=0.92, help="Fixed Ru value.")
    args = parser.parse_args()

    setup_root_style()
    ensure_dir(args.plot_dir)

    f_data = ROOT.TFile.Open(args.data, "READ")
    if not f_data or f_data.IsZombie():
        raise RuntimeError(f"Could not open input file: {args.data}")

    # Load 2D profiles
    p2_mpf = require_object(f_data, "MPF_2D")
    p2_db = require_object(f_data, "DB_2D")
    p2_r0 = require_object(f_data, "HDM_r0_2D")
    p2_r1 = require_object(f_data, "HDM_r1_2D")
    p2_rn = require_object(f_data, "HDM_rn_2D")
    p2_ru = require_object(f_data, "HDM_ru_2D")

    # Barrel projections
    h_mpf = project_barrel_profile_y(p2_mpf, "barrel_MPF", args.eta_min, args.eta_max)
    h_db = project_barrel_profile_y(p2_db, "barrel_DB", args.eta_min, args.eta_max)

    h_r0 = project_barrel_profile_y(p2_r0, "barrel_HDM_r0", args.eta_min, args.eta_max)
    h_r1 = project_barrel_profile_y(p2_r1, "barrel_HDM_r1", args.eta_min, args.eta_max)
    h_rn = project_barrel_profile_y(p2_rn, "barrel_HDM_rn", args.eta_min, args.eta_max)
    h_ru = project_barrel_profile_y(p2_ru, "barrel_HDM_ru", args.eta_min, args.eta_max)

    # HDM definitions
    h_hdm_mpf = build_hdm_hist(
        "barrel_HDM_MPF_branch",
        h_r0,
        h_r1,
        h_rn,
        h_ru,
        args.rn,
        args.ru,
        branch="mpf",
    )

    h_hdm_db = build_hdm_hist(
        "barrel_HDM_DB_branch",
        h_r0,
        h_r1,
        h_rn,
        h_ru,
        args.rn,
        args.ru,
        branch="db",
    )

    # Ratios/checks
    h_r0_over_mpf = clone_ratio("barrel_HDM_r0_over_MPF", h_r0, h_mpf)
    h_r1_over_db = clone_ratio("barrel_HDM_r1_over_DB", h_r1, h_db)

    h_hdm_mpf_over_mpf = clone_ratio("barrel_HDM_MPF_over_MPF", h_hdm_mpf, h_mpf)
    h_hdm_db_over_db = clone_ratio("barrel_HDM_DB_over_DB", h_hdm_db, h_db)
    h_hdm_mpf_over_hdm_db = clone_ratio("barrel_HDM_MPF_over_HDM_DB", h_hdm_mpf, h_hdm_db)

    # Closure check
    h_closure = h_r0.Clone("barrel_HDM_closure_r0_minus_r1_minus_rn_minus_ru")
    h_closure.SetDirectory(0)
    h_closure.Reset()
    h_closure.SetStats(False)
    format_x_axis(h_closure)

    for b in range(1, h_closure.GetNbinsX() + 1):
        value = (
            h_r0.GetBinContent(b)
            - h_r1.GetBinContent(b)
            - h_rn.GetBinContent(b)
            - h_ru.GetBinContent(b)
        )

        # Conservative uncorrelated uncertainty
        error = math.sqrt(
            h_r0.GetBinError(b) ** 2
            + h_r1.GetBinError(b) ** 2
            + h_rn.GetBinError(b) ** 2
            + h_ru.GetBinError(b) ** 2
        )

        h_closure.SetBinContent(b, value)
        h_closure.SetBinError(b, error)

    # Plots
    eta_label = f"{args.eta_min:.3f} < #eta < {args.eta_max:.3f}"

    # draw_overlay(
    #     os.path.join(args.plot_dir, "barrel_response_mpf_db_hdm.png"),
    #     [h_mpf, h_db, h_hdm_mpf, h_hdm_db],
    #     ["MPF", "DB", "HDM MPF", "HDM DB"],
    #     f"Data barrel response, {eta_label}",
    #     "Response",
    #     y_range=hist_y_range([h_mpf, h_db, h_hdm_mpf, h_hdm_db], default=(0.5, 1.5)),
    # )

    draw_overlay(
        os.path.join(args.plot_dir, "barrel_response_mpf_db_hdm.png"),
        [h_mpf, h_db, h_hdm_mpf],
        ["MPF", "DB", "HDM"],
        f"Data barrel response, {eta_label}",
        "Response",
        y_range=hist_y_range([h_mpf, h_db, h_hdm_mpf], default=(0.5, 1.5)),
    )

    draw_overlay(
        os.path.join(args.plot_dir, "barrel_hdm_over_baseline.png"),
        [h_hdm_mpf_over_mpf, h_hdm_db_over_db, h_hdm_mpf_over_hdm_db],
        ["HDM MPF / MPF", "HDM DB / DB", "HDM MPF / HDM DB"],
        f"HDM ratio checks, {eta_label}",
        "Ratio",
        y_range=hist_y_range(
            [h_hdm_mpf_over_mpf, h_hdm_db_over_db, h_hdm_mpf_over_hdm_db],
            default=(0.8, 1.2),
        ),
    )

    draw_overlay(
        os.path.join(args.plot_dir, "barrel_hdm_inputs_vs_existing.png"),
        [h_r0_over_mpf, h_r1_over_db],
        ["HDM r0 / MPF", "HDM r1 / DB"],
        f"Input consistency checks, {eta_label}",
        "Ratio",
        y_range=(0.8, 1.2),
    )

    draw_overlay(
        os.path.join(args.plot_dir, "barrel_hdm_components.png"),
        [h_r0, h_r1, h_rn, h_ru],
        ["r0", "r1", "rn", "ru"],
        f"HDM components, {eta_label}",
        "Projected component",
        y_range=hist_y_range([h_r0, h_r1, h_rn, h_ru], default=(-0.2, 1.2)),
        draw_line_at_one=False,
    )

    draw_overlay(
        os.path.join(args.plot_dir, "barrel_hdm_closure.png"),
        [h_closure],
        ["r0 - r1 - rn - ru"],
        f"HDM closure, {eta_label}",
        "Closure",
        y_range=hist_y_range([h_closure], default=(-0.05, 0.05)),
        draw_line_at_one=False,
    )

    # Output ROOT file
    fout = ROOT.TFile.Open(args.out, "RECREATE")

    write_all(
        fout,
        [
            h_mpf,
            h_db,
            h_r0,
            h_r1,
            h_rn,
            h_ru,
            h_hdm_mpf,
            h_hdm_db,
            h_r0_over_mpf,
            h_r1_over_db,
            h_hdm_mpf_over_mpf,
            h_hdm_db_over_db,
            h_hdm_mpf_over_hdm_db,
            h_closure,
        ],
    )

    fout.Close()
    f_data.Close()

    print("Done.")
    print(f"Input file: {args.data}")
    print(f"Output ROOT file: {args.out}")
    print(f"Plots written to: {args.plot_dir}")
    print(f"Barrel projection: {args.eta_min} < eta < {args.eta_max}")
    print(f"Fixed HDM nuisance inputs: Rn={args.rn}, Ru={args.ru}")


if __name__ == "__main__":
    main()