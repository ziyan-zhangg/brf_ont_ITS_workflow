"""Fixture test for qc_report.py and the demux_qc.py Q histogram (ITS).

Ported from the 16S workflow's test (v6.00). Builds a small fake ITS run:
per-sample demultiplexed FASTQs profiled by the real demux_qc.run_demux_qc,
and read_counts_summary.txt written by the production writer, in post and in
pre mode. Then runs the report end to end with the repo's qc_config.yaml and
checks the Q bins and median_q, the metrics, the flags (Q and modal_dev_bp
never flagged), the parser (incl. pre-mode NA handling), run-level relative
mode, controls, the old-config warnings, client safety and the
self-contained qc_report.html, plus the ITS-specific config handling:
the ITS config loads, a core window below 1 kb is accepted, an info_only
metric never sets a flag, and the HTML title / labels come from the config.

Needs PyYAML (the config is YAML). matplotlib is optional; without it the
plot checks are skipped.

Run:  python tests/test_qc_report.py
"""
import contextlib
import getpass
import gzip
import importlib.util
import io
import json
import math
import re
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Must never appear in an HTML report: facility paths, project, username
# (the username only when long enough not to match by accident).
INTERNAL = ("/g/data", "/scratch", "/home", "vz35",
            *([getpass.getuser()] if len(getpass.getuser()) >= 4 else []))

sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location("qc_report", ROOT / "qc_report.py")
qc = importlib.util.module_from_spec(spec)
sys.modules["qc_report"] = qc
spec.loader.exec_module(qc)
core = qc.core
from cutadapt_2pass import CutadaptStats  # noqa: E402
from demux_qc import run_demux_qc  # noqa: E402

Q40, Q20, Q17, Q10 = "I", "5", "2", "+"

# Production ITS windows: Chopper 400-2000 bp (run_script.qsub), which the
# config mirrors as minknow_min_len / chopper_max_len.
MIN_LEN, MAX_LEN = 400, 2000
PEAK, SHOULDER = 790, 650          # ITS1F/ITS4 main peak and shoulder


def outside_core(n):
    """n reads spread over 1,000-1,175 bp: inside the delivered window but
    above the 450-1000 core window, and spread so none becomes modal."""
    return [(1, 1000 + 25 * (i % 8), Q40) for i in range(n)]


# sample -> list of (n_reads, length, qual_char)
DESIGN = {
    **{f"S{i:02d}": [(300, PEAK, Q40)] for i in range(1, 10)},
    "S10": [(300, PEAK, Q20)],                          # exactly Q20: counts
    "Mixed": [(200, PEAK, Q40), (100, SHOULDER, Q40)],  # peak + shoulder: PASS
    "BP_core": [(140, PEAK, Q40), *outside_core(160)],  # core 47%
    "BP_long": [(270, PEAK, Q40), (30, 2600, Q40)],     # 10% > 2 kb
    "BP_q": [(150, PEAK, Q40), (150, PEAK, Q17)],       # median Q17: not flagged
    "ModalOff": [(300, 550, Q40)],                      # modal 240 bp off: info only
    "BP_primer": [(300, PEAK, Q40)],                    # 70% primer found
    "NTC1": [(20, PEAK, Q40)],
    "Tiny": [(100, PEAK, Q40)],                         # n < 200
    "POS_B": [(300, PEAK, Q40)],                        # ClientB's own control
    "unk": [(50, PEAK, Q10)],
}
SAMPLES = [s for s in DESIGN if s != "unk"]
# Sample sheet: ClientB owns BP_core, ModalOff and the control POS_B; NTC1 has
# no client, so it is a shared control; everything else is ClientA.
CLIENT = {s: "ClientA" for s in SAMPLES} | {
    "BP_core": "ClientB", "ModalOff": "ClientB", "POS_B": "ClientB", "NTC1": ""}
TYPE = {"NTC1": "NTC", "POS_B": "POS_CTRL"}
COMMENT = {"Tiny": "low DNA input", "ModalOff": "clientB private note"}
BARCODE = {s: f"ITS{i + 1:02d}" for i, s in enumerate(SAMPLES)}
PASS1 = {s: 97.0 for s in SAMPLES} | {"BP_primer": 70.0}
DELIVERED = {s: 20000 for s in SAMPLES} | {"NTC1": 40, "Tiny": 100}


def median_low(xs):
    xs = sorted(xs)
    return xs[(len(xs) - 1) // 2] if xs else None


def mean_q(qual: str) -> float:
    err = math.fsum(10 ** (-(ord(c) - 33) / 10) for c in qual) / len(qual)
    return -10 * math.log10(err)


def reads_of(sid):
    for n, length, q in DESIGN[sid]:
        for _ in range(n):
            yield length, q


def write_demux_fastqs(d: Path) -> None:
    d.mkdir(parents=True)
    for sid in DESIGN:
        with (d / f"sample_{sid}.fastq").open("w") as fh:
            for k, (length, q) in enumerate(reads_of(sid)):
                fh.write(f"@r{k}\n{'A' * length}\n+\n{q * length}\n")


def chopper_counts(sid):
    """(raw, kept) as production Chopper (-q 15, 400-2000 bp) would."""
    raw = kept = 0
    for length, q in reads_of(sid):
        raw += 1
        kept += MIN_LEN <= length <= MAX_LEN and mean_q(q * 10) >= 15
    return raw, kept


def build_run(tmp: Path, name: str, stage: str, long_name=None) -> Path:
    """Run dir with demux_qc/ and run_log_*/read_counts_summary.txt."""
    out = tmp / name
    raw_demux = tmp / f"{name}_demux"
    write_demux_fastqs(raw_demux)
    with contextlib.redirect_stdout(io.StringIO()):
        qc_by = run_demux_qc(raw_demux, out / "demux_qc", MIN_LEN, MAX_LEN,
                             filter_stage=stage)

    integ = out / "integrated_demultiplexing"
    (integ / "ClientA").mkdir(parents=True)
    delivered = dict(DELIVERED)
    if long_name:
        delivered[long_name] = 12345
    for sid, n in delivered.items():
        with gzip.open(integ / "ClientA" / f"{sid}.fastq.gz", "wt") as fh:
            fh.write("@r\nA\n+\nI\n" * n)
    (integ / "sample_unk.fastq").write_text("@r\nA\n+\nI\n" * 7)

    pass1 = dict(PASS1)
    if long_name:
        pass1[long_name] = 91.2
    cut = {sid: CutadaptStats(sid, input_reads=1000,
                              pass1_kept=int(10 * pct))
           for sid, pct in pass1.items()}
    chop = ({sid: chopper_counts(sid) for sid in SAMPLES}
            if stage == "post" else None)
    log_dir = out / "run_log_20261013"
    log_dir.mkdir(parents=True)
    with contextlib.redirect_stdout(io.StringIO()):
        core.write_run_summary(integ, log_dir / "read_counts_summary.txt",
                               out / "raw", None, out, total_input=10**6,
                               filter_stage=stage, qc_by_sid=qc_by,
                               chopper_by_sid=chop, cutadapt_by_sid=cut)
    return out


def write_sheet(path: Path) -> None:
    path.write_text("Client,Sample_ID,Location,Barcode,Sample_Type,Comment\n" + "".join(
        f"{CLIENT[s]},{s},A{i + 1},{BARCODE[s]},{TYPE.get(s, '')},"
        f"{COMMENT.get(s, '')}\n"
        for i, s in enumerate(SAMPLES)))


def expect_exit(fn, *args) -> str:
    try:
        fn(*args)
    except SystemExit as e:
        return str(e)
    return ""


def table(path: Path, key: int = 1) -> dict:
    rows = [ln.split("\t") for ln in path.read_text().splitlines()
            if not ln.startswith("#")]
    head = rows[0]
    return {r[key]: dict(zip(head, r)) for r in rows[1:]}


def run_report(out: Path, sheet: Path, *extra, config=None) -> int:
    with contextlib.redirect_stdout(io.StringIO()):
        return qc.main(["--output-dir", str(out), "--samplesheet", str(sheet),
                        "--config", str(config or ROOT / "qc_config.yaml"),
                        *extra])


def html_blob(path: Path) -> dict:
    m = re.search(r'<script id="qc-data" type="application/json">(.*?)</script>',
                  path.read_text(encoding="utf-8"), re.S)
    return json.loads(m.group(1)) if m else {}


def says_2kb(text: str) -> bool:
    """True if the text says '2 kb' ('1.2 kb' does not count). Used on a page
    built from a config whose chopper_max_len is not 2000, so any '2 kb' left
    would be a hard-coded label."""
    return re.search(r"(?<![\d.])2 kb", text) is not None


def variant_config(tmp: Path, name: str, **subs) -> Path:
    """Copy of the repo config with exact-text substitutions applied."""
    text = (ROOT / "qc_config.yaml").read_text()
    for old, new in subs.values():
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    path = tmp / name
    path.write_text(text)
    return path


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="qc_report_test_"))
    failures = []

    def check(ok, msg_pass, msg_fail):
        if ok:
            print(f"PASS  {msg_pass}")
        else:
            failures.append(msg_fail)

    sheet = tmp / "sheet.csv"
    write_sheet(sheet)

    # ---- 0. ITS config -----------------------------------------------------
    cfg = qc.load_config(ROOT / "qc_config.yaml")
    check(cfg["config_version"] == "1.1-its" and cfg["assay"] == "ITS"
          and cfg["title"] == "ITS QC report"
          and cfg["filter"] == {"minknow_min_len": MIN_LEN, "minknow_min_q": 15,
                                "chopper_max_len": MAX_LEN}
          and cfg["core_window"] == {"min": 450, "max": 1000}
          and cfg["plot"] == {"len_min": 400, "len_max": 3000}
          and cfg["metrics"]["modal_dev_bp"]["info_only"]
          and qc.flagged_metrics(cfg) == ["delivered_reads", "primer_found",
                                          "long_frac", "core_frac"]
          and not cfg["_warnings"],
          "ITS config loads: 1.1-its, 400-2000 bp, core 450-1000, plot "
          "400-3000, modal_dev_bp info_only, no warnings",
          f"ITS config: {cfg}")

    narrow = variant_config(tmp, "narrow.yaml", cw=(
        "core_window: {min: 450, max: 1000}", "core_window: {min: 600, max: 900}"))
    c2 = qc.load_config(narrow)
    check(c2["core_window"] == {"min": 600, "max": 900},
          "core window entirely below 1 kb (600-900) accepted",
          f"narrow core window: {c2['core_window']}")
    for label, bad_cw, want in (
            ("below the delivered window", "{min: 350, max: 1000}", "must lie inside"),
            ("above the delivered window", "{min: 450, max: 2050}", "must lie inside"),
            ("off the 50 bp grid", "{min: 425, max: 1000}", "multiples of the 50")):
        p = variant_config(tmp, "bad_cw.yaml", cw=(
            "core_window: {min: 450, max: 1000}", f"core_window: {bad_cw}"))
        msg = expect_exit(qc.load_config, p)
        check(want in msg, f"core window {label} rejected",
              f"core window {label}: {msg!r}")
    for label, bad_plot in (("plot ending inside the delivered window",
                             "{len_min: 400, len_max: 2000}"),
                            ("plot starting above minknow_min_len",
                             "{len_min: 450, len_max: 3000}")):
        p = variant_config(tmp, "bad_plot.yaml", plot=(
            "plot: {len_min: 400, len_max: 3000}", f"plot: {bad_plot}"))
        msg = expect_exit(qc.load_config, p)
        check("plot range" in msg, f"{label} rejected", f"{label}: {msg!r}")

    # A config written before assay/title/plot/info_only existed (the 16S
    # v2.2 shape) still loads, with the 16S defaults.
    legacy = tmp / "legacy_16s.yaml"
    legacy.write_text(
        "config_version: 2.2\n"
        "filter: {minknow_min_len: 1000, minknow_min_q: 15, chopper_max_len: 2000}\n"
        "core_window: {min: 1300, max: 1700}\n"
        "min_reads_for_stats: 200\n"
        "metrics:\n"
        "  delivered_reads: {direction: higher_better, check: 15000, fail: 5000, absolute_only: true}\n"
        "  primer_found:    {direction: higher_better, check: 0.90, fail: 0.75}\n"
        "  long_frac:       {direction: lower_better,  check: 0.05, fail: 0.15}\n"
        "  core_frac:       {direction: higher_better, check: 0.75, fail: 0.50}\n"
        "  modal_dev_bp:    {direction: lower_better,  check: 50,   fail: 100, absolute_only: true}\n"
        "run_level: {outlier_k_check: 3, outlier_k_fail: 5, mad_floor: 0.02, min_samples: 8}\n"
        "controls: {types: [NTC, BLANK, POS_CTRL]}\n")
    c16 = qc.load_config(legacy)
    check(c16["title"] == "16S QC report" and qc.long_label(c16) == "> 2 kb"
          and c16["plot"] == {"len_min": 1000, "len_max": 4000}
          and qc.flagged_metrics(c16) == list(qc.METRICS),
          "legacy 16S config: title '16S QC report', '> 2 kb', 1-4 kb axis, "
          "all five metrics flagged",
          f"legacy config: {c16['title']} {qc.long_label(c16)} {c16['plot']}")

    # ---- 1. demux_qc Q bins ------------------------------------------------
    out = build_run(tmp, "run_post", "post")
    dq = out / "demux_qc"
    bad = []
    for sid in DESIGN:
        p = qc.read_demux_profile(dq, sid)
        n = sum(1 for _ in reads_of(sid))
        n_win = sum(1 for L, _ in reads_of(sid) if MIN_LEN <= L <= MAX_LEN)
        direct = median_low(int(mean_q(q * L) + 1e-6) for L, q in reads_of(sid)
                            if MIN_LEN <= L <= MAX_LEN)
        if (sum(p.q_all.values()) != n or sum(p.q_win.values()) != n_win
                or qc.median_q(p) != direct):
            bad.append((sid, sum(p.q_all.values()), n, qc.median_q(p), direct))
    check(not bad, f"Q bins sum to reads and median_q matches a direct "
                   f"calculation for {len(DESIGN)} bins",
          f"Q bins / median_q mismatch: {bad}")
    summ = table(dq / "demux_qc_summary.tsv", key=0)
    check(summ["S10"]["median_q"] == "20" and summ["BP_q"]["median_q"] == "17"
          and summ["S01"]["window_min"] == "400"
          and summ["S01"]["window_max"] == "2000",
          "demux_qc_summary.tsv median_q column and 400-2000 window",
          f"summary: S10={summ['S10']} BP_q={summ['BP_q'].get('median_q')}")

    # ---- 2. post-mode end to end -------------------------------------------
    rc = run_report(out, sheet)
    qd = out / "qc_report"
    check(rc == 0, "post-mode end-to-end run completed", f"qc.main returned {rc}")
    t = table(qd / "per_sample_qc.tsv")
    want = {"S01": "PASS", "S10": "PASS", "Mixed": "PASS", "BP_core": "FAIL",
            "BP_long": "CHECK", "BP_q": "PASS", "ModalOff": "PASS",
            "BP_primer": "FAIL", "NTC1": "CTRL", "Tiny": "FAIL"}
    got = {s: t[s]["flag"] for s in want}
    check(got == want, f"flags {got}", f"flags: want {want}, got {got}")
    expect_reasons = {
        "BP_core": "FAIL: LOW_CORE(47%)",
        "BP_long": "CHECK: LONG(10%)",
        "BP_q": "",
        "ModalOff": "",
        "Mixed": "",
        "BP_primer": "FAIL: PRIMER(70%)",
        "Tiny": "FAIL: READS(100)",
    }
    wrong = {s: t[s]["reasons"] for s, r in expect_reasons.items()
             if t[s]["reasons"] != r}
    check(not wrong, "reason codes as designed", f"reasons: {wrong}")
    s01 = t["S01"]
    check(s01["raw_demuxed"] == "300" and s01["delivered_reads"] == "20000"
          and s01["primer_found"] == "0.9700" and s01["long_frac"] == "0.0000"
          and s01["core_frac"] == "1.0000" and s01["modal_length"] == str(PEAK)
          and s01["modal_dev_bp"] == "0" and s01["median_q"] == "40",
          "S01 metrics exact, median_q informational",
          f"S01 row: {s01}")
    check(t["Mixed"]["core_frac"] == "1.0000"
          and t["Mixed"]["modal_length"] == str(PEAK),
          "peak + 650 bp shoulder: both inside the core window, modal = peak",
          f"Mixed row: {t['Mixed']}")
    check(t["BP_q"]["median_q"] == "17" and t["Tiny"]["median_q"] == "NA (n<200)",
          "median_q reported but never flagged; NA below min reads",
          f"median_q: BP_q={t['BP_q']} Tiny={t['Tiny']['median_q']}")
    check(t["BP_long"]["core_frac"] == "1.0000",
          "core_frac counts only reads inside 400-2000 bp",
          f"BP_long core_frac: {t['BP_long']['core_frac']}")
    check(t["Tiny"]["core_frac"] == "NA (n<200)"
          and t["Tiny"]["primer_found"] == "NA (n<200)",
          "n<200: stats NA, delivered_reads still flagged",
          f"Tiny row: {t['Tiny']}")
    check(t["NTC1"]["reasons"] == "" and t["NTC1"]["flag"] == "CTRL",
          "control flagged CTRL with no reasons", f"NTC1 row: {t['NTC1']}")

    # info_only: the value is reported, but never a level, flag or reason.
    check(t["ModalOff"]["modal_dev_bp"] == str(PEAK - 550)
          and t["ModalOff"]["modal_length"] == "550",
          f"info_only modal_dev_bp still reported ({PEAK - 550} bp)",
          f"ModalOff row: {t['ModalOff']}")

    head = (qd / "per_sample_qc.tsv").read_text()
    check("# config_version\t1.1-its" in head and "# filter_stage\tpost" in head
          and "# assay\tITS" in head
          and "# modal_dev_bp\tinformational, not flagged" in head
          and "> 2 kb reads visible" in head,
          "TSV header records assay, config version, filter stage, info_only",
          "TSV header lacks assay / config version / filter stage / info_only")
    check((qd / "per_client" / "ClientA_qc.tsv").is_file(),
          "per-client table written", "per_client/ClientA_qc.tsv missing")
    with gzip.open(qd / "qc_bins.tsv.gz", "rt") as fh:
        kinds = {ln.split("\t")[1] for ln in fh.read().splitlines()[1:]}
    check(kinds == {"len", "q", "q_all"}, f"qc_bins.tsv.gz kinds {sorted(kinds)}",
          f"qc_bins.tsv.gz kinds: {kinds}")
    run_qc = (qd / "run_qc.txt").read_text()
    check("Filter stage: post" in run_qc and "unk" in run_qc
          and "NOTE: pre-mode" not in run_qc and "(ITS)" in run_qc
          and re.search(r"modal_dev_bp .* info only", run_qc)
          and "MODAL" not in run_qc,
          "run_qc.txt: post stage, catch-all bins, modal_dev_bp 'info only'",
          f"run_qc.txt:\n{run_qc}")

    # ---- 2b. qc_report.html -----------------------------------------------
    page = (qd / "qc_report.html").read_text(encoding="utf-8")
    blob = html_blob(qd / "qc_report.html")
    ids = {s["id"] for s in blob.get("samples", [])}
    check(ids == set(SAMPLES), f"qc_report.html JSON parses, {len(ids)} samples",
          f"qc_report.html samples: {sorted(ids)}")
    check("http://" not in page and "https://" not in page
          and "<script src" not in page and "<link" not in page,
          "qc_report.html is self-contained (no network references)",
          "qc_report.html references something external")
    hs = {s["id"]: s for s in blob["samples"]}
    mby = {m["key"]: m for m in blob["metrics"]}
    check(hs["BP_core"]["flag"] == "FAIL" and hs["BP_core"]["levels"]["core_frac"] == 2
          and sum(hs["S01"]["len"]) == 300 and sum(hs["S01"]["q"]) == 300
          and [m["key"] for m in blob["metrics"]] == list(qc.METRICS)
          and mby["modal_dev_bp"]["info_only"] and not mby["core_frac"]["info_only"]
          and hs["BP_q"]["values"]["median_q"] == 17,
          "qc_report.html data: flags, per-metric levels, bins, info_only marked",
          f"qc_report.html data: {hs['BP_core']['levels']} {blob['metrics']}")
    check(all("modal_dev_bp" not in s["levels"] for s in blob["samples"])
          and hs["ModalOff"]["values"]["modal_dev_bp"] == PEAK - 550
          and hs["ModalOff"]["flag"] == "PASS",
          "info_only metric has a value but no level in any sample",
          f"modal_dev_bp levels: {[s['levels'] for s in blob['samples']]}")
    check(page.count("__QC_DATA__") == 0 and page.count('id="qc-data"') == 1
          and page.count("__TITLE__") == 0,
          "qc_report.html template placeholders filled once",
          "qc_report.html placeholder left in or duplicated")

    # Title and labels come from the config, never the 16S hard-codes.
    a = blob["axes"]["len"]
    check("<title>ITS QC report - run_post</title>" in page
          and blob["run"]["title"] == "ITS QC report"
          and blob["run"]["long_label"] == "> 2 kb"
          and mby["long_frac"]["label"] == "Long reads (> 2 kb)"
          and "16S" not in page and "1.2 kb" not in page
          and (a["min"], a["max"], a["step"]) == (400, 3000, 50),
          "HTML title, long-read label and length axis from the ITS config; "
          "no '16S' / old '1.2 kb' anywhere",
          f"title/labels: {blob['run'].get('title')} {blob['run'].get('long_label')} "
          f"axes {a}; '16S' in page: {'16S' in page}; '1.2 kb' in page: {'1.2 kb' in page}")
    custom = variant_config(tmp, "custom.yaml",
                            title=("title: ITS QC report", "title: Fungal ITS2 QC"),
                            maxlen=("chopper_max_len: 2000", "chopper_max_len: 1500"),
                            plot=("plot: {len_min: 400, len_max: 3000}",
                                  "plot: {len_min: 0, len_max: 2500}"))
    run_report(out, sheet, config=custom)
    cpage = (qd / "qc_report.html").read_text(encoding="utf-8")
    cblob = html_blob(qd / "qc_report.html")
    check("<title>Fungal ITS2 QC - run_post</title>" in cpage
          and cblob["run"]["long_label"] == "> 1.5 kb"
          and not says_2kb(cpage)
          and (cblob["axes"]["len"]["min"], cblob["axes"]["len"]["max"]) == (0, 2500),
          "changing title / chopper_max_len / plot in the config changes the page",
          f"custom config page: {cblob['run'].get('title')} "
          f"{cblob['run'].get('long_label')} {cblob['axes']['len']}")
    run_report(out, sheet)                      # restore the default report

    try:
        import matplotlib  # noqa: F401
        pngs = {p.name for p in qd.glob("*.png")}
        check(pngs == {"length_heatmap.png", "q_heatmap.png"}
              and not (qd / "flagged").exists(),
              f"plots: {sorted(pngs)}, no flagged/ by default",
              f"plots: {pngs}, flagged/ exists: {(qd / 'flagged').exists()}")
        run_report(out, sheet, "--png-drilldown")
        flagged = {p.stem for p in (qd / "flagged").glob("*.png")}
        want_f = {s for s, f in want.items() if f in ("CHECK", "FAIL")}
        check(flagged == want_f, f"--png-drilldown: flagged {sorted(flagged)}",
              f"--png-drilldown flagged: {flagged}")
    except ImportError:
        print("SKIP  plots (matplotlib not installed)")

    # HTML must not depend on matplotlib.
    saved = {k: sys.modules[k] for k in list(sys.modules)
             if k == "matplotlib" or k.startswith("matplotlib.")}
    sys.modules["matplotlib"] = None          # makes 'import matplotlib' fail
    try:
        (qd / "qc_report.html").unlink()
        run_report(out, sheet)
    finally:
        del sys.modules["matplotlib"]
        sys.modules.update(saved)
    check((qd / "qc_report.html").is_file()
          and "matplotlib not available" in (qd / "run_qc.txt").read_text(),
          "qc_report.html still written without matplotlib",
          "qc_report.html missing when matplotlib is unavailable")

    # ---- 2c. old config entries: pct_q20, thresholds on an info_only metric -
    old_cfg = tmp / "old_config.yaml"
    old_cfg.write_text((ROOT / "qc_config.yaml").read_text().replace(
        "metrics:\n", "metrics:\n  pct_q20: {direction: higher_better, "
        "check: 0.60, fail: 0.40}\n", 1).replace(
        "modal_dev_bp:    {direction: lower_better,  info_only: true}",
        "modal_dev_bp:    {direction: lower_better,  check: 50, fail: 100, "
        "info_only: true}"))
    rc = run_report(out, sheet, config=old_cfg)
    rq = (qd / "run_qc.txt").read_text()
    t3 = table(qd / "per_sample_qc.tsv")
    check(rc == 0 and "metrics.pct_q20 ignored: Q is no longer flagged" in rq
          and "metrics.modal_dev_bp is info_only; its check/fail thresholds "
              "are ignored" in rq
          and t3["BP_q"]["flag"] == "PASS" and t3["ModalOff"]["flag"] == "PASS",
          "old config: pct_q20 and thresholds on an info_only metric warn, "
          "nothing is flagged by them",
          f"old config: rc={rc} BP_q={t3['BP_q']['flag']} "
          f"ModalOff={t3['ModalOff']['flag']}")

    # ---- 2d. qc_report.html is client-safe ---------------------------------
    # A run under a path holding every facility marker, with an old config
    # and a second counts file beside it so warnings carry absolute paths.
    leak = build_run(tmp, "g/data/vz35/scratch/home/ONT_ITS_LEAK", "post")
    shutil.copy(leak / "run_log_20261013" / "read_counts_summary.txt", leak)
    leak_cfg = leak / "old_config.yaml"
    leak_cfg.write_text(old_cfg.read_text())
    run_report(leak, sheet, config=leak_cfg)
    lpage = (leak / "qc_report" / "qc_report.html").read_text(encoding="utf-8")
    flat = lpage.replace("\\\\", "/").replace("\\", "/")
    found = [s for s in INTERNAL if s in flat]
    lblob = html_blob(leak / "qc_report" / "qc_report.html")
    check(not found and "output_dir" not in lblob["run"]
          and lblob["run"]["name"] == "ONT_ITS_LEAK"
          and any("old_config.yaml: metrics.pct_q20" in w
                  for w in lblob["run"]["warnings"]),
          "qc_report.html has no facility paths (/g/data, /scratch, /home, "
          "vz35); run name kept, warnings cut to file names",
          f"qc_report.html leaks {found}; run={lblob.get('run')}")
    ls = {s["id"]: s for s in lblob["samples"]}
    check(ls["Tiny"]["comment"] == "low DNA input" and ls["S01"]["comment"] == ""
          and ls["NTC1"]["type"] == "NTC" and ls["NTC1"]["flag"] == "CTRL",
          "qc_report.html data keeps Comment and Sample_Type (controls -> CTRL)",
          f"comment/type: {ls['Tiny']} {ls['NTC1']}")
    check('{ key: "comment", label: "Comment"' in lpage
          and '{ key: "type", label: "Type"' not in lpage,
          "sample table shows Comment instead of Type",
          "sample table columns not as expected")
    check('key: "reasons"' not in lpage
          and ls["Tiny"]["reasons"] == "FAIL: READS(100)"
          and "<b>Reasons:</b>" in lpage and "esc(s.reasons)" in lpage
          and table(leak / "qc_report" / "per_sample_qc.tsv")["Tiny"]["reasons"]
          == "FAIL: READS(100)",
          "no Reasons column in the HTML table; reasons kept in the embedded "
          "data (tooltip / drill-down) and in per_sample_qc.tsv",
          "Reasons column still in the HTML table, or reasons missing elsewhere")

    leak_pages = sorted((leak / "integrated_demultiplexing").glob("*/*.html"))
    check(len(leak_pages) == 2, "per-client HTML present in the leak run",
          f"leak run per-client HTML: {leak_pages}")
    for f in leak_pages:
        fp = f.read_text(encoding="utf-8").replace("\\\\", "/").replace("\\", "/")
        found = [s for s in INTERNAL if s in fp]
        check(not found, f"{f.name} has no facility paths (/g/data, /scratch, "
                         f"/home, vz35)",
              f"{f.name} leaks {found}")

    # No integrated_demultiplexing/: per-client HTML falls back to per_client/.
    integ_l = leak / "integrated_demultiplexing"
    integ_l.rename(leak / "moved_away")
    try:
        run_report(leak, sheet)
        fb = sorted(p.name for p in (leak / "qc_report" / "per_client").glob("*.html"))
        rq = (leak / "qc_report" / "run_qc.txt").read_text()
    finally:
        (leak / "moved_away").rename(integ_l)
    check(fb == ["ClientA_qc_report.html", "ClientB_qc_report.html"]
          and "written to qc_report/per_client/ instead" in rq,
          "no integrated_demultiplexing/: per-client HTML in per_client/, with a warning",
          f"fallback: {fb}")

    # ---- 2e. per-client reports --------------------------------------------
    pc = qd / "per_client"
    integ = out / "integrated_demultiplexing"
    run_flags = {s["id"]: s["flag"] for s in html_blob(qd / "qc_report.html")["samples"]}
    placed = sorted(p.relative_to(integ).as_posix() for p in integ.glob("*/*.html"))
    check(placed == ["ClientA/ClientA_qc_report.html", "ClientB/ClientB_qc_report.html"]
          and not list(pc.glob("*.html")) and not list(integ.glob("*.html")),
          "<Client>_qc_report.html written into integrated_demultiplexing/<Client>/, "
          "none for shared controls, none left in per_client/",
          f"per-client HTML: {placed}, per_client/: {sorted(p.name for p in pc.glob('*'))}")
    for me, other in (("ClientA", "ClientB"), ("ClientB", "ClientA")):
        path = integ / me / f"{me}_qc_report.html"
        page = path.read_text(encoding="utf-8")
        blob = html_blob(path)
        ids = {s["id"] for s in blob["samples"]}
        mine = {s for s in SAMPLES if CLIENT[s] == me} | {"NTC1"}
        theirs = [s for s in SAMPLES if CLIENT[s] == other]
        leaked = [x for x in (other, *(f'"{s}"' for s in theirs),
                              *(f'"{BARCODE[s]}"' for s in theirs),
                              *(COMMENT[s] for s in theirs if s in COMMENT))
                  if x in page]
        check(ids == mine and not leaked,
              f"{me} report: all {me} samples + shared control, nothing of {other}",
              f"{me} report: ids {sorted(ids)} (want {sorted(mine)}), leaked {leaked}")
        diff = {s["id"]: (s["flag"], run_flags[s["id"]]) for s in blob["samples"]
                if s["flag"] != run_flags[s["id"]]}
        counts = blob["run"]["flag_counts"]
        check(not diff and sum(counts.values()) == len(mine)
              and counts["CTRL"] == len({"NTC1", "POS_B"} & mine)
              and blob["run"]["client"] == me and not blob["run"]["warnings"]
              and not blob["run"]["catch_all"]
              and blob["run"]["title"] == "ITS QC report",
              f"{me} report: flags identical to the run report, counts "
              f"recomputed ({counts}), no warnings / catch-all",
              f"{me} report: flag diffs {diff}, counts {counts}, run {blob['run']}")
    check((pc / "ClientA_qc.tsv").is_file() and "NTC1" in table(pc / "ClientA_qc.tsv")
          and "POS_B" not in table(pc / "ClientA_qc.tsv"),
          "per-client TSV: shared control included, other client's control not",
          f"ClientA_qc.tsv rows: {sorted(table(pc / 'ClientA_qc.tsv'))}")
    tp = (qd / "per_sample_qc.tsv")
    check(table(tp)["NTC1"]["client"] == "" and table(tp)["NTC1"]["flag"] == "CTRL",
          "control with empty Client kept in the run table as CTRL",
          f"NTC1 row: {table(tp).get('NTC1')}")

    # ---- 3. pre mode -------------------------------------------------------
    out_pre = build_run(tmp, "run_pre", "pre")
    c = qc.parse_read_counts(out_pre / "run_log_20261013" / "read_counts_summary.txt")
    check(c.filter_stage == "pre" and c.rows["S01"].raw_demuxed is None
          and c.rows["S01"].pct_retained is None
          and c.rows["S01"].post_chopper == 300,
          "parser: pre-mode table with blank raw_demuxed / pct_retained",
          f"parser pre-mode: {c.filter_stage} {c.rows.get('S01')}")
    rc = run_report(out_pre, sheet)
    tp = table(out_pre / "qc_report" / "per_sample_qc.tsv")
    rq = (out_pre / "qc_report" / "run_qc.txt").read_text()
    check(rc == 0 and tp["S01"]["long_frac"] == "NA (pre mode)"
          and tp["BP_long"]["long_frac"] == "NA (pre mode)"
          and tp["S01"]["core_frac"] == "1.0000" and "NOTE: pre-mode run" in rq
          and "reads > 2 kb were never demultiplexed" in rq,
          "pre mode: long_frac NA (pre mode), noted at top of run_qc.txt "
          "with the ITS long-read label",
          f"pre-mode S01: {tp.get('S01')}\n{rq}")
    ppage = (out_pre / "qc_report" / "qc_report.html").read_text(encoding="utf-8")
    pblob = html_blob(out_pre / "qc_report" / "qc_report.html")
    check(pblob["run"]["filter_stage"] == "pre"
          and pblob["run"]["long_label"] == "> 2 kb" and "1.2 kb" not in ppage,
          "pre-mode HTML takes its long-read label from the config",
          f"pre-mode HTML label: {pblob['run'].get('long_label')}")

    # ---- 4. counts parser --------------------------------------------------
    cp = qc.parse_read_counts(out / "run_log_20261013" / "read_counts_summary.txt")
    r = cp.rows["BP_long"]
    check(cp.filter_stage == "post" and r.raw_demuxed == 300
          and r.post_chopper == 270 and r.pct_retained == 90.0
          and r.pct_pass1_kept == 97.0,
          "parser: post-mode raw_demuxed / post_chopper / pct_retained",
          f"parser post-mode: {cp.filter_stage} {r}")
    check("TOTAL DEMULTIPLEXED" not in cp.rows and "unk" in cp.rows,
          "parser: TOTAL rows ignored", f"parser rows: {sorted(cp.rows)}")

    long_name = "X" * 55
    out_long = build_run(tmp, "run_long", "post", long_name=long_name)
    rl = qc.parse_read_counts(out_long / "run_log_20261013" / "read_counts_summary.txt")
    check(rl.rows[long_name].reads == 12345
          and rl.rows[long_name].pct_pass1_kept == 91.2
          and rl.rows["BP_long"].pct_retained == 90.0,
          "parser: over-long sample name shifts columns correctly",
          f"parser long name: {rl.rows.get(long_name)}")

    counts_txt = (out / "run_log_20261013" / "read_counts_summary.txt").read_text()
    broken = tmp / "broken.txt"
    broken.write_text(counts_txt.replace("pct_pass1_kept", "pct_kept"))
    msg = expect_exit(qc.parse_read_counts, broken)
    check("pct_pass1_kept" in msg, "parser: missing column fails loudly by name",
          f"parser missing column message: {msg!r}")

    old = tmp / "old.txt"
    old.write_text("Sample                                        Reads % of Input\n"
                   "----------------------------------------\n"
                   "sample_S01                                    45249      4.00%\n")
    msg = expect_exit(qc.parse_read_counts, old)
    check("missing column" in msg, "parser: legacy 3-column table rejected",
          f"parser legacy message: {msg!r}")

    # ---- 5. demux_qc without Q histogram (older runs) ----------------------
    for q in (out / "demux_qc").glob("*.q.tsv"):
        q.unlink()
    run_report(out, sheet)
    t2 = table(qd / "per_sample_qc.tsv")
    rq = (qd / "run_qc.txt").read_text()
    check(t2["S01"]["median_q"] == "NA" and t2["BP_q"]["flag"] == "PASS"
          and "median_q is NA" in rq,
          "no .q.tsv: median_q NA with a warning, other metrics unaffected",
          f"no-Q run: {t2['S01']} / {t2['BP_q']['flag']}")
    msg = expect_exit(run_report, out, sheet, "--demux-qc-dir", str(tmp / "nope"))
    check("demux QC dir not found" in msg, "missing demux_qc/ fails loudly",
          f"missing demux_qc message: {msg!r}")

    # ---- 6. run-level relative mode ----------------------------------------
    sheet_rows = [qc.SheetSample("C", f"R{i:02d}", order=i) for i in range(12)]

    def synthetic(primer_pct, reads=30000, modal=lambda i: PEAK):
        rows = {s.sample_id: qc.CountsRow(s.sample_id, reads, 1000, 1000, 100.0,
                                          970, primer_pct(i))
                for i, s in enumerate(sheet_rows)}
        profs = {s.sample_id: qc.DemuxProfile(
                     s.sample_id, reads=1000, filter_stage="post",
                     modal_length=modal(i), len_bins={750: 1000}, has_q=True,
                     q_win={17: 1000})
                 for i, s in enumerate(sheet_rows)}
        return qc.Counts(Path("x"), "post", rows), profs

    # Whole run at ~70% primer found (FAIL), one sample far below the rest.
    counts, profs = synthetic(lambda i: 30.0 if i == 0 else 68.0 + i % 5)
    ev = qc.evaluate(sheet_rows, counts, profs, cfg, [])
    p_run = next(m for m in ev.runs if m.metric == "primer_found")
    flags = {r.sample.sample_id: r.flag for r in ev.results}
    check(p_run.run_level == 2 and p_run.relative
          and any(ln.startswith("RUN: LOW_PRIMER (median ") and ", FAIL)" in ln
                  for ln in qc.run_flag_lines(ev.runs))
          and not any("Q20" in ln for ln in qc.run_flag_lines(ev.runs)),
          "run flag RUN: LOW_PRIMER raised and metric switched to relative mode",
          f"run-level primer_found: {p_run}")
    check(flags["R00"] == "FAIL"
          and all(f == "PASS" for s, f in flags.items() if s != "R00"),
          "relative mode flags only the outlier, not the whole run",
          f"relative flags: {flags}")
    check("z=" in ev.results[0].reasons, f"outlier reason: {ev.results[0].reasons}",
          f"outlier reason lacks z: {ev.results[0].reasons}")

    # delivered_reads stays absolute even when the whole run is low
    counts, profs = synthetic(lambda i: 97.0)
    for c in counts.rows.values():
        c.reads = 9000
    ev = qc.evaluate(sheet_rows, counts, profs, cfg, [])
    rd = next(m for m in ev.runs if m.metric == "delivered_reads")
    check(rd.run_level == 1 and not rd.relative
          and all(r.levels["delivered_reads"] == 1 for r in ev.results),
          "delivered_reads: run flag raised but samples stay absolute",
          f"delivered_reads run: {rd}")

    # info_only never flags, even when the modal lengths are wildly spread
    # (a run median far off for every sample, and one extreme outlier) --
    # under the same config with info_only removed the same data does flag.
    spread = lambda i: 400 + 60 * i if i else 2000   # noqa: E731
    counts, profs = synthetic(lambda i: 97.0, modal=spread)
    ev = qc.evaluate(sheet_rows, counts, profs, cfg, [])
    md = next(m for m in ev.runs if m.metric == "modal_dev_bp")
    check(md.median is not None and md.run_level == 0 and not md.relative
          and all("modal_dev_bp" not in r.levels for r in ev.results)
          and all(r.flag == "PASS" and "MODAL" not in r.reasons for r in ev.results)
          and not any("MODAL" in ln for ln in qc.run_flag_lines(ev.runs)),
          f"info_only modal_dev_bp: median {md.median:.0f} bp reported, no "
          f"level, flag, reason or run flag for any sample",
          f"info_only modal_dev_bp: {md} {[(r.flag, r.reasons) for r in ev.results]}")
    flagging = variant_config(tmp, "modal_flagged.yaml", m=(
        "modal_dev_bp:    {direction: lower_better,  info_only: true}",
        "modal_dev_bp:    {direction: lower_better,  check: 50, fail: 100, "
        "absolute_only: true}"))
    ev2 = qc.evaluate(sheet_rows, counts, profs, qc.load_config(flagging), [])
    check(any("MODAL" in r.reasons for r in ev2.results),
          "same modal spread does flag once info_only is removed (control)",
          f"modal flagged config: {[(r.flag, r.reasons) for r in ev2.results]}")

    # controls are left out of run medians
    counts, profs = synthetic(lambda i: 97.0)
    ctrl = qc.SheetSample("C", "POS", sample_type="POS_CTRL", order=99)
    counts.rows["POS"] = qc.CountsRow("POS", 30000, 1000, 1000, 100.0, 970, 97.0)
    profs["POS"] = qc.DemuxProfile("POS", reads=1000, modal_length=1150,
                                   len_bins={750: 1000}, has_q=True,
                                   q_win={17: 1000})
    ev = qc.evaluate([*sheet_rows, ctrl], counts, profs, cfg, [])
    pos = ev.results[-1]
    check(ev.modal_median == PEAK and pos.flag == "CTRL"
          and pos.values["modal_dev_bp"] == 1150 - PEAK,
          "controls: flag CTRL, metrics reported, excluded from run medians",
          f"control: median {ev.modal_median}, {pos}")

    shutil.rmtree(tmp, ignore_errors=True)
    print()
    if failures:
        for f in failures:
            print("FAIL  " + f)
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
