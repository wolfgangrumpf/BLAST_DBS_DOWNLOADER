#!/usr/bin/env python3
"""
download_blastdb.py - list and download pre-formatted NCBI BLAST databases.

Run with NO arguments:
    Fetches the list of the latest-version BLAST databases from the NCBI
    FTP site and prints them (name, type, download size, last update,
    description). If databases are already installed in the output directory
    and you're at an interactive terminal, it then offers to update them.

Run with one or more DATABASE NAMES:
    For each database it looks up the volume archives, downloads every
    .tar.gz volume with wget, verifies each archive against its NCBI .md5
    checksum, and then extracts (untars) the volumes.

    Use -j/--jobs N to download N volumes at a time (faster for big,
    multi-volume databases like nt). Keep N modest (3-4) to be a good
    citizen of NCBI's servers.

Run with -u/--update:
    Detects the databases already present in --outdir and refreshes them to
    the latest version, re-downloading only the volumes whose checksum
    changed upstream. Non-interactive, so it's well suited to cron.

Examples:
    ./download_blastdb.py                    # list dbs (+ offer to update local)
    ./download_blastdb.py nt                 # download the 'nt' database
    ./download_blastdb.py 16S_ribosomal_RNA  # download a small db
    ./download_blastdb.py nt -j 4 -o blastdb # 4 volumes at once, into ./blastdb
    ./download_blastdb.py nt nr -j 3 --keep
    ./download_blastdb.py --update -o /data/blastdb -j 4 -J 2   # cron refresh
"""

import argparse
import glob
import hashlib
import json
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import time

NCBI_BASE = "https://ftp.ncbi.nlm.nih.gov/blast/db"
AGGREGATE_METADATA_URL = f"{NCBI_BASE}/blastdb-metadata-1-1.json"
# Molecule-type suffixes used in the per-database metadata file names.
MOLECULE_TYPES = ("nucl", "prot")

# Serializes console output so parallel workers don't garble each other's lines.
_PRINT_LOCK = threading.Lock()


def log(msg):
    with _PRINT_LOCK:
        print(msg, flush=True)


class DownloadError(Exception):
    """Raised when a volume cannot be downloaded or verified."""


def fmt_dur(seconds):
    """Format a duration as e.g. '45s', '9m32s', '1h03m'."""
    seconds = int(round(max(0.0, seconds)))
    if seconds < 60:
        return f"{seconds}s"
    minutes, sec = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes}m{sec:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h{minutes:02d}m"


def timestamp():
    """Local wall-clock timestamp, e.g. '2026-06-26 14:32:07 EDT'."""
    return time.strftime("%Y-%m-%d %H:%M:%S %Z").strip()


def mb_per_sec(num_bytes, seconds):
    """Megabytes per second (decimal: 1 MB = 1e6 bytes) for a transfer."""
    if seconds <= 0 or num_bytes <= 0:
        return 0.0
    return num_bytes / 1e6 / seconds


class ProgressTracker:
    """Thread-safe progress/timing/ETA for a download -> extract pipeline.

    The two stages are accounted separately: download bytes/time feed the
    transfer-speed figures, while extraction time is tracked on its own. The
    ETA is based on the download stage (usually the bottleneck): with J download
    workers, finishing R not-yet-downloaded volumes takes about
    R * (avg download time) / J. Skipped and failed volumes are excluded from
    the averages so a restart full of instant skips doesn't skew the estimate.
    """

    def __init__(self, total, jobs):
        self.total = total
        self.jobs = max(1, jobs)
        self.lock = threading.Lock()
        self.done = 0          # fully resolved volumes (printed lines)
        self.skipped = 0       # up-to-date / no-op volumes (no download)
        self.dl_count = 0      # volumes actually downloaded + verified
        self.dl_time = 0.0     # summed wall time of those downloads
        self.dl_bytes = 0      # summed bytes pulled over the wire
        self.ext_count = 0     # volumes extracted
        self.ext_time = 0.0    # summed extraction time
        self.t0 = time.monotonic()

    def note_download(self, dl_seconds, dl_bytes):
        """Record a completed download+verify (called by a download worker)."""
        with self.lock:
            self.dl_count += 1
            self.dl_time += dl_seconds
            self.dl_bytes += dl_bytes

    def complete(self, basename, status, *, skipped=False, ok=True,
                 dl_seconds=0.0, dl_bytes=0, ext_seconds=0.0):
        """Emit the final line for a volume and update completion stats."""
        with self.lock:
            self.done += 1
            n = self.done
            if skipped:
                self.skipped += 1
            if ext_seconds > 0:
                self.ext_count += 1
                self.ext_time += ext_seconds
            remaining = self.total - n
            remaining_dl = max(0, self.total - self.dl_count - self.skipped)
            avg_dl = (self.dl_time / self.dl_count) if self.dl_count else None
            avg_speed = mb_per_sec(self.dl_bytes, self.dl_time) if self.dl_time else None

        if skipped or not ok:
            parts = ["skipped" if skipped else "failed"]
        else:
            seg = []
            if dl_bytes > 0 and dl_seconds > 0:
                seg.append(f"dl {mb_per_sec(dl_bytes, dl_seconds):.1f} MB/s")
            if ext_seconds > 0:
                seg.append(f"extract {fmt_dur(ext_seconds)}")
            parts = [", ".join(seg) if seg else "ok"]

        if avg_speed:
            parts.append(f"avg {avg_speed:.1f} MB/s")
        if remaining > 0 and avg_dl:
            parts.append(f"ETA {fmt_dur(remaining_dl * avg_dl / self.jobs)}")
        if remaining == 0:
            parts.append(f"total {fmt_dur(time.monotonic() - self.t0)}")

        log(f"{timestamp()}  [{n:>3}/{self.total}] {basename}: {status}  "
            f"[{' | '.join(parts)}]")

    def elapsed(self):
        return time.monotonic() - self.t0


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def die(msg, code=1):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(code)


def require_wget():
    if shutil.which("wget") is None:
        die("wget is not installed or not on PATH; please install wget first.")


def wget_fetch_text(url):
    """Download a (small) URL and return its body as text, or None on failure."""
    proc = subprocess.run(
        ["wget", "-q", "--tries=3", "--timeout=60", "-O", "-", url],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    if proc.returncode != 0 or not proc.stdout:
        return None
    return proc.stdout.decode("utf-8", errors="replace")


def wget_download(url, outdir, quiet=False):
    """
    Download `url` into `outdir` with wget (resuming partial files).
    Returns the path to the downloaded file on success, else None.
    When quiet, wget is silent (used for parallel downloads so the live
    progress bars of several files don't collide); otherwise it shows a bar.
    """
    basename = url.rstrip("/").rsplit("/", 1)[-1]
    dest = os.path.join(outdir, basename)
    cmd = ["wget", "-c", "--tries=3", "--timeout=60"]
    cmd.append("-q" if quiet else "--show-progress")
    cmd.append(url)
    proc = subprocess.run(cmd, cwd=outdir)
    if proc.returncode != 0:
        return None
    return dest


def md5_of_file(path, chunk_size=1024 * 1024):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_md5_file(text):
    """NCBI .md5 files look like:  '<hexdigest>  <filename>'."""
    if not text:
        return None
    token = text.split()[0].strip().lower()
    return token if len(token) == 32 and all(c in "0123456789abcdef" for c in token) else None


# --------------------------------------------------------------------------- #
# Metadata handling
# --------------------------------------------------------------------------- #
def extract_db_records(obj):
    """
    Normalize the aggregate metadata JSON into a list of per-database dicts,
    being tolerant of the exact top-level shape.
    """
    if isinstance(obj, list):
        return [r for r in obj if isinstance(r, dict) and "dbname" in r]

    if isinstance(obj, dict):
        for value in obj.values():
            if isinstance(value, list) and any(
                isinstance(r, dict) and "dbname" in r for r in value
            ):
                return [r for r in value if isinstance(r, dict) and "dbname" in r]
        records = []
        for key, value in obj.items():
            if isinstance(value, dict):
                rec = dict(value)
                rec.setdefault("dbname", key)
                records.append(rec)
        if records:
            return records

    return []


def fetch_database_list():
    text = wget_fetch_text(AGGREGATE_METADATA_URL)
    if text is None:
        die(
            "could not download the database list from NCBI "
            f"({AGGREGATE_METADATA_URL}). Check your network connection."
        )
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        die(f"could not parse the NCBI metadata file: {exc}")
    records = extract_db_records(data)
    if not records:
        die("the NCBI metadata file did not contain any recognizable databases.")
    return records


def human_gb(num_bytes):
    try:
        return f"{float(num_bytes) / 1e9:6.2f}"
    except (TypeError, ValueError):
        return "   ?  "


def render_database_table(records):
    """Print the aligned NAME / TYPE / SIZE / UPDATED / DESCRIPTION table."""
    records = sorted(records, key=lambda r: r.get("dbname", "").lower())

    name_w = max([len("NAME")] + [len(r.get("dbname", "")) for r in records])
    name_w = min(name_w, 36)

    header = (
        f"{'NAME':<{name_w}}  {'TYPE':<10}  {'SIZE/GB':>7}  "
        f"{'UPDATED':<10}  DESCRIPTION"
    )
    print(header)
    print("-" * len(header))

    for r in records:
        name = r.get("dbname", "")[:name_w]
        dbtype = (r.get("dbtype") or "")[:10]
        size = human_gb(r.get("bytes-total-compressed", r.get("bytes-total")))
        updated = (r.get("last-updated") or "")[:10]
        desc = r.get("description") or ""
        if len(desc) > 70:
            desc = desc[:67] + "..."
        print(f"{name:<{name_w}}  {dbtype:<10}  {size:>7}  {updated:<10}  {desc}")


def print_database_list(records):
    render_database_table(records)
    print()
    print(f"{len(records)} databases available.")
    print("Download one or more with:  download_blastdb.py NAME [NAME ...]")


# --------------------------------------------------------------------------- #
# Resolving which archive files belong to a database
# --------------------------------------------------------------------------- #
def archive_url_from_entry(entry):
    """
    Turn a metadata 'files' entry into a concrete NCBI HTTPS URL.
    Entries may be bare names or ftp://, https://, s3:// or gs:// URLs;
    we always pull the file from the NCBI FTP HTTPS endpoint.
    """
    basename = str(entry).rstrip("/").rsplit("/", 1)[-1]
    return f"{NCBI_BASE}/{basename}"


def resolve_archives(db_name, dbtype_hint=None):
    """
    Return the list of .tar.gz archive URLs for `db_name` by reading its
    per-database metadata JSON. Tries the hinted molecule type first.
    Returns [] if the database cannot be found.
    """
    order = list(MOLECULE_TYPES)
    if dbtype_hint:
        hint = "nucl" if dbtype_hint.lower().startswith("nuc") else "prot"
        order = [hint] + [m for m in MOLECULE_TYPES if m != hint]

    for mol in order:
        meta_url = f"{NCBI_BASE}/{db_name}-{mol}-metadata.json"
        text = wget_fetch_text(meta_url)
        if not text:
            continue
        try:
            meta = json.loads(text)
        except json.JSONDecodeError:
            continue
        files = meta.get("files") or []
        archives = [
            archive_url_from_entry(f) for f in files if str(f).endswith(".tar.gz")
        ]
        if archives:
            return archives
    return []


# --------------------------------------------------------------------------- #
# Download + verify + extract
# --------------------------------------------------------------------------- #
def fetch_expected_md5(archive_url):
    text = wget_fetch_text(archive_url + ".md5")
    return parse_md5_file(text)


# Sentinel meaning "no checksum supplied; fetch it" (distinct from None = no
# .md5 exists on the server).
_FETCH_MD5 = object()


def md5_marker_path(archive_path):
    """Path of the on-disk .md5 marker kept beside a volume (NCBI's convention)."""
    return archive_path + ".md5"


def read_stored_md5(archive_path):
    """Return the md5 recorded in the local .md5 marker, or None if absent."""
    try:
        with open(md5_marker_path(archive_path)) as fh:
            return parse_md5_file(fh.read())
    except OSError:
        return None


def write_md5_marker(archive_path, md5hex):
    """Persist the verified md5 next to the volume as a version marker."""
    basename = os.path.basename(archive_path)
    try:
        with open(md5_marker_path(archive_path), "w") as fh:
            fh.write(f"{md5hex}  {basename}\n")
    except OSError:
        pass


def already_extracted(stem, outdir):
    """
    True if files extracted from a volume already exist (e.g. 'nt.12.nhr'
    for the stem 'nt.12'), ignoring the archive and its .md5.
    """
    pattern = os.path.join(outdir, glob.escape(stem) + ".*")
    for p in glob.glob(pattern):
        if p.endswith(".tar.gz") or p.endswith(".tar.gz.md5"):
            continue
        return True
    return False


# Trailing volume number on an archive/file stem, e.g. the '.12' in 'nt.12'.
_VOLUME_SUFFIX = re.compile(r"\.\d{2,4}$")


def installed_databases(outdir):
    """Best-effort list of BLAST database names already present in `outdir`.

    Detected from (in order of reliability): our own '<db>.NN.tar.gz.md5'
    version markers, BLAST alias files ('<db>.nal'/'.pal' for multi-volume
    databases), and single-volume index files ('<db>.nin'/'.pin').
    """
    names = set()

    for p in glob.glob(os.path.join(outdir, "*.tar.gz.md5")):
        stem = os.path.basename(p)[: -len(".tar.gz.md5")]   # '<db>' or '<db>.NN'
        names.add(_VOLUME_SUFFIX.sub("", stem))

    for ext in (".nal", ".pal"):                            # multi-volume aliases
        for p in glob.glob(os.path.join(outdir, "*" + ext)):
            names.add(os.path.basename(p)[: -len(ext)])

    for ext in (".nin", ".pin"):                            # single-volume indexes
        for p in glob.glob(os.path.join(outdir, "*" + ext)):
            stem = os.path.basename(p)[: -len(ext)]
            if not _VOLUME_SUFFIX.search(stem):             # skip volume members
                names.add(stem)

    return sorted(n for n in names if n)


def confirm(question, default=False):
    """Ask a yes/no question on an interactive terminal. Non-TTY -> default."""
    if not sys.stdin.isatty():
        return default
    suffix = " [Y/n] " if default else " [y/N] "
    try:
        answer = input(question + suffix).strip().lower()
    except (EOFError, KeyboardInterrupt):
        print()
        return default
    if not answer:
        return default
    return answer in ("y", "yes")


def ask_int(question, default, minimum=1, maximum=None):
    """Ask for a whole number on an interactive terminal. Non-TTY -> default.

    Re-prompts on invalid input; empty input accepts the default.
    """
    if not sys.stdin.isatty():
        return default
    while True:
        try:
            raw = input(f"{question} [{default}] ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return default
        if not raw:
            return default
        try:
            value = int(raw)
        except ValueError:
            print("  please enter a whole number")
            continue
        if value < minimum:
            print(f"  must be at least {minimum}")
            continue
        if maximum is not None and value > maximum:
            print(f"  must be at most {maximum}")
            continue
        return value


def download_and_verify_archive(archive_url, outdir, verbose=True,
                                expected=_FETCH_MD5):
    """Download one .tar.gz volume and verify it against its .md5.

    Returns (path, xfer_seconds, xfer_bytes). xfer_* describe the bytes
    actually pulled over the wire and how long that took (resume-aware), so
    callers can report transmission speed. Both are 0 if nothing was fetched.

    On successful verification the checksum is persisted as a '<archive>.md5'
    marker so a later run can tell whether an already-extracted volume is still
    the current upstream version.

    `expected` may be passed in to reuse an already-fetched checksum; by default
    it is fetched here. (None means the server has no .md5 for this volume.)

    Raises DownloadError on unrecoverable failure (safe to call from threads).
    """
    basename = archive_url.rsplit("/", 1)[-1]
    dest = os.path.join(outdir, basename)

    if expected is _FETCH_MD5:
        expected = fetch_expected_md5(archive_url)
    if expected is None and verbose:
        log(f"  ! warning: no .md5 found for {basename}; cannot verify checksum")

    if expected and os.path.exists(dest) and md5_of_file(dest) == expected:
        if verbose:
            log(f"  = {basename} already present and verified; skipping download")
        write_md5_marker(dest, expected)
        return dest, 0.0, 0

    xfer_seconds = 0.0
    xfer_bytes = 0
    for attempt in (1, 2):
        if verbose:
            note = "" if attempt == 1 else " (re-downloading from scratch)"
            log(f"  > downloading {basename}{note}")
        size_before = os.path.getsize(dest) if os.path.exists(dest) else 0
        t0 = time.monotonic()
        path = wget_download(archive_url, outdir, quiet=not verbose)
        xfer_seconds += time.monotonic() - t0
        if path is None:
            raise DownloadError(f"wget failed to download {archive_url}")
        # Bytes added this attempt = final size minus whatever resume started from.
        xfer_bytes += max(0, os.path.getsize(path) - size_before)

        if expected is None:
            return path, xfer_seconds, xfer_bytes

        actual = md5_of_file(path)
        if actual == expected:
            if verbose:
                log(f"  + checksum OK ({actual})")
            write_md5_marker(path, expected)
            return path, xfer_seconds, xfer_bytes

        if verbose:
            log(
                f"  ! checksum MISMATCH for {basename}: "
                f"expected {expected}, got {actual}"
            )
        try:
            os.remove(path)
        except OSError:
            pass

    raise DownloadError(f"checksum verification failed for {basename}")


def extract_archive(path, outdir):
    import tarfile

    with tarfile.open(path, "r:gz") as tar:
        try:  # safe 'data' filter where available (Python 3.12+)
            tar.extractall(path=outdir, filter="data")
        except TypeError:
            tar.extractall(path=outdir)


def download_stage(archive_url, outdir, keep, extract, verbose, tracker):
    """Download + verify one volume (the network-bound stage).

    Returns one of:
      ("skip", message)                  -> nothing to do (already current)
      ("extract", path, dl_s, dl_b)      -> verified archive ready to extract
      ("done", status, dl_s, dl_b)       -> verified, no extraction requested
      ("fail", message)                  -> could not download/verify
    Does NOT touch the tracker; the caller reports, so download and extract
    milestones are recorded in the right stage.
    """
    basename = archive_url.rsplit("/", 1)[-1]
    stem = basename[:-7] if basename.endswith(".tar.gz") else basename
    archive_path = os.path.join(outdir, basename)
    prefetched_md5 = _FETCH_MD5

    try:
        # Restart currency check: archive gone but extracted files present.
        if extract and not keep and not os.path.exists(archive_path) \
                and already_extracted(stem, outdir):
            current = fetch_expected_md5(archive_url)
            stored = read_stored_md5(archive_path)
            if current is None:
                return ("skip", "already extracted; skipping (no upstream checksum)")
            if stored == current:
                return ("skip", "already extracted, up to date; skipping")
            if verbose:
                reason = ("upstream version changed" if stored is not None
                          else "no local version marker")
                log(f"  i {basename}: {reason}; re-downloading to refresh")
            prefetched_md5 = current

        path, dl_s, dl_b = download_and_verify_archive(
            archive_url, outdir, verbose=verbose, expected=prefetched_md5
        )
        if extract:
            return ("extract", path, dl_s, dl_b)
        return ("done", "ok (verified)", dl_s, dl_b)

    except DownloadError as exc:
        return ("fail", str(exc))


def extract_stage(path, outdir, keep):
    """Decompress one verified archive (the CPU-bound stage). Returns seconds."""
    t0 = time.monotonic()
    extract_archive(path, outdir)
    if not keep:
        try:
            os.remove(path)  # the .md5 marker stays as the version record
        except OSError:
            pass
    return time.monotonic() - t0


def download_databases(names, outdir, keep_archives, extract, jobs, extract_jobs):
    os.makedirs(outdir, exist_ok=True)

    # One aggregate fetch lets us hint the molecule type and catch typos early.
    records = []
    agg_text = wget_fetch_text(AGGREGATE_METADATA_URL)
    if agg_text:
        try:
            records = extract_db_records(json.loads(agg_text))
        except json.JSONDecodeError:
            records = []
    dbtype_by_name = {
        r.get("dbname"): r.get("dbtype") for r in records if r.get("dbname")
    }

    # Resolve every requested database into a flat list of volume tasks.
    tasks = []                 # list of (archive_url, db_name)
    volumes_by_db = {}         # db_name -> count of volumes
    failures = []
    for name in names:
        archives = resolve_archives(name, dbtype_by_name.get(name))
        if not archives:
            log(f"! could not find a database named '{name}' on the NCBI site")
            if records:
                hits = [
                    r["dbname"] for r in records
                    if name.lower() in r.get("dbname", "").lower()
                ][:6]
                if hits:
                    log(f"    did you mean: {', '.join(hits)} ?")
            failures.append(name)
            continue
        volumes_by_db[name] = len(archives)
        for url in archives:
            tasks.append((url, name))

    if not tasks:
        die(f"nothing to download; unresolved: {', '.join(failures) or '(none)'}")

    total = len(tasks)
    jobs = max(1, jobs)
    extract_jobs = max(1, extract_jobs) if extract else 0
    # Show wget's live progress bar only when nothing else prints concurrently.
    verbose = jobs == 1 and not extract

    pipe_note = ""
    if extract:
        pipe_note = f", extracting {extract_jobs} at a time"
    log(
        f"{timestamp()}  starting: {total} volume(s) across "
        f"{len(volumes_by_db)} database(s); downloading {jobs} at a time"
        f"{pipe_note} into '{outdir}'\n"
    )

    tracker = ProgressTracker(total, jobs)

    work_q = queue.Queue()
    for task in tasks:
        work_q.put(task)
    # Bounded handoff queue = backpressure: downloads may run ahead of
    # extraction, but only by a few volumes, so archives can't pile up and
    # exhaust disk on a big database like nt.
    handoff = queue.Queue(maxsize=jobs + extract_jobs) if extract else None

    results = []  # (db_name, ok, message)
    results_lock = threading.Lock()

    def record(db_name, ok, message):
        with results_lock:
            results.append((db_name, ok, message))

    def downloader():
        while True:
            try:
                url, db_name = work_q.get_nowait()
            except queue.Empty:
                return
            basename = url.rsplit("/", 1)[-1]
            try:
                outcome = download_stage(
                    url, outdir, keep_archives, extract, verbose, tracker
                )
                kind = outcome[0]
                if kind == "skip":
                    tracker.complete(basename, outcome[1], skipped=True)
                    record(db_name, True, None)
                elif kind == "fail":
                    tracker.complete(basename, f"FAILED -- {outcome[1]}", ok=False)
                    record(db_name, False, outcome[1])
                elif kind == "done":  # verified, no extraction
                    _, status, dl_s, dl_b = outcome
                    tracker.note_download(dl_s, dl_b)
                    tracker.complete(basename, status, dl_seconds=dl_s, dl_bytes=dl_b)
                    record(db_name, True, None)
                else:  # ("extract", path, dl_s, dl_b)
                    _, path, dl_s, dl_b = outcome
                    tracker.note_download(dl_s, dl_b)
                    handoff.put((db_name, basename, path, dl_s, dl_b))
            except Exception as exc:  # noqa: BLE001 - never let a worker die silently
                tracker.complete(basename, f"FAILED -- {exc}", ok=False)
                record(db_name, False, str(exc))
            finally:
                work_q.task_done()

    def extractor():
        while True:
            item = handoff.get()
            try:
                if item is None:  # shutdown sentinel
                    return
                db_name, basename, path, dl_s, dl_b = item
                try:
                    ext_s = extract_stage(path, outdir, keep_archives)
                    tracker.complete(basename, "ok (verified, extracted)",
                                     dl_seconds=dl_s, dl_bytes=dl_b, ext_seconds=ext_s)
                    record(db_name, True, None)
                except Exception as exc:  # noqa: BLE001 - keep the pipeline alive
                    tracker.complete(basename, f"FAILED -- extract: {exc}", ok=False)
                    record(db_name, False, str(exc))
            finally:
                handoff.task_done()

    # Start the download workers (and extract workers if extracting).
    dl_threads = [threading.Thread(target=downloader, name=f"dl-{i}")
                  for i in range(jobs)]
    ex_threads = [threading.Thread(target=extractor, name=f"ex-{i}")
                  for i in range(extract_jobs)]
    for t in ex_threads:
        t.start()
    for t in dl_threads:
        t.start()

    for t in dl_threads:           # all downloads finished/queued for extract
        t.join()
    if extract:
        for _ in ex_threads:       # tell extractors to stop once drained
            handoff.put(None)
        for t in ex_threads:
            t.join()

    # Per-database summary: a database is done only if all its volumes succeeded.
    bad_volumes = {}
    for name, ok, _msg in results:
        if not ok:
            bad_volumes[name] = bad_volumes.get(name, 0) + 1

    # Safety net: every queued volume must have produced exactly one result.
    unaccounted = total - len(results)
    if unaccounted > 0:
        log(f"! WARNING: {unaccounted} volume(s) were not accounted for "
            f"(unexpected worker error)")

    log("\n--- summary ---")
    for name in names:
        if name in failures:
            log(f"  {name}: NOT FOUND")
        elif name in bad_volumes:
            log(f"  {name}: FAILED ({bad_volumes[name]}/{volumes_by_db[name]} volume(s) bad)")
        else:
            log(f"  {name}: ok ({volumes_by_db[name]} volume(s))")

    elapsed = tracker.elapsed()
    if tracker.dl_count:
        gb = tracker.dl_bytes / 1e9
        avg_dl = tracker.dl_time / tracker.dl_count
        per_stream = mb_per_sec(tracker.dl_bytes, tracker.dl_time)
        aggregate = mb_per_sec(tracker.dl_bytes, elapsed)
        log(
            f"  ---\n  {tracker.dl_count} volume(s), {gb:.2f} GB downloaded "
            f"in {fmt_dur(elapsed)} (avg {fmt_dur(avg_dl)} download/volume at -j {jobs})"
        )
        speed_line = f"  download speed: {per_stream:.1f} MB/s per stream"
        if jobs > 1:
            speed_line += f"; {aggregate:.1f} MB/s effective (all streams)"
        log(speed_line)
        if tracker.ext_count:
            avg_ext = tracker.ext_time / tracker.ext_count
            log(f"  extraction: {tracker.ext_count} volume(s), "
                f"avg {fmt_dur(avg_ext)}/volume across {extract_jobs} worker(s)")
    else:
        log(f"  ---\n  no new volumes downloaded (all already present); "
            f"elapsed {fmt_dur(elapsed)}")
    log(f"  finished {timestamp()}")

    if failures or bad_volumes or unaccounted > 0:
        die("one or more databases did not download cleanly (see summary above)")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser():
    parser = argparse.ArgumentParser(
        description="List or download pre-formatted NCBI BLAST databases.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="With no database names, the available databases are listed.",
    )
    parser.add_argument(
        "databases", nargs="*",
        help="one or more BLAST database names to download (e.g. nt nr pdbaa)",
    )
    parser.add_argument(
        "-l", "--list", action="store_true",
        help="list available databases (the default when no names are given)",
    )
    parser.add_argument(
        "-o", "--outdir", default=".",
        help="directory to download/extract into (default: current directory)",
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=1, metavar="N",
        help="download N volumes at a time (default: 1; 3-4 is a good, "
             "server-friendly value for large databases like nt)",
    )
    parser.add_argument(
        "-J", "--extract-jobs", type=int, default=None, metavar="N",
        help="decompress N volumes at a time in parallel with downloading "
             "(default: same as --jobs). Downloads never wait on extraction; "
             "raise this if your CPU can't keep up with the download rate.",
    )
    parser.add_argument(
        "-u", "--update", action="store_true",
        help="refresh databases already present in --outdir to the latest "
             "version (only volumes whose checksum changed are downloaded). "
             "Non-interactive, so it's suitable for cron.",
    )
    parser.add_argument(
        "--keep", action="store_true",
        help="keep the downloaded .tar.gz archives after extracting them",
    )
    parser.add_argument(
        "--no-extract", action="store_true",
        help="download and verify only; do not untar the archives",
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    require_wget()

    if args.jobs < 1:
        die("--jobs must be 1 or greater")
    if args.extract_jobs is not None and args.extract_jobs < 1:
        die("--extract-jobs must be 1 or greater")
    extract_jobs = args.extract_jobs if args.extract_jobs is not None else args.jobs

    def run(names, jobs=args.jobs, extract_jobs=extract_jobs):
        download_databases(
            names=names,
            outdir=args.outdir,
            keep_archives=args.keep,
            extract=not args.no_extract,
            jobs=jobs,
            extract_jobs=extract_jobs,
        )

    # Explicit, non-interactive update of whatever is installed in --outdir.
    if args.update and not args.databases:
        names = installed_databases(args.outdir)
        if not names:
            log(f"{timestamp()}  update: no installed databases found in "
                f"'{args.outdir}'; nothing to do.")
            return
        log(f"{timestamp()}  update: refreshing {len(names)} database(s) in "
            f"'{args.outdir}': {', '.join(names)}")
        run(names)
        return

    # Explicit database names (download, or refresh just those).
    if args.databases:
        run(args.databases)
        return

    # No names given: show the catalog.
    records = fetch_database_list()
    print_database_list(records)
    if args.list:
        return

    # Then, if databases are installed here, show them in the same table format
    # and offer to update them.
    installed = installed_databases(args.outdir)
    if not installed:
        return
    by_name = {r.get("dbname"): r for r in records}
    installed_records = [
        by_name.get(name, {"dbname": name,
                            "description": "(not listed in current NCBI catalog)"})
        for name in installed
    ]
    where = "the current directory" if args.outdir == "." else f"'{args.outdir}'"
    noun = "database" if len(installed) == 1 else "databases"
    print(f"\nInstalled in {where} ({len(installed)} {noun}):")
    render_database_table(installed_records)
    print()

    if confirm(f"Update {'this' if len(installed) == 1 else 'these'} "
               f"{len(installed)} {noun} to the latest version now?"):
        chosen_jobs = ask_int(
            "How many volumes to download at once? (3-4 is plenty; "
            "high values may get throttled by NCBI)",
            default=args.jobs, minimum=1,
        )
        chosen_ej = (args.extract_jobs if args.extract_jobs is not None
                     else chosen_jobs)
        run(installed, jobs=chosen_jobs, extract_jobs=chosen_ej)
    elif not sys.stdin.isatty():
        log(f"({len(installed)} {noun} installed; run with --update to "
            f"refresh them, e.g. from cron.)")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        die("interrupted", code=130)