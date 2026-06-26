#!/usr/bin/env python3
"""
download_blastdb.py - list and download pre-formatted NCBI BLAST databases.

Run with NO arguments:
    Fetches the list of the latest-version BLAST databases from the NCBI
    FTP site and prints them (name, type, download size, last update,
    description).

Run with one or more DATABASE NAMES:
    For each database it looks up the volume archives, downloads every
    .tar.gz volume with wget, verifies each archive against its NCBI .md5
    checksum, and then extracts (untars) the volumes.

    Use -j/--jobs N to download N volumes at a time (faster for big,
    multi-volume databases like nt). Keep N modest (3-4) to be a good
    citizen of NCBI's servers.

Examples:
    ./download_blastdb.py                    # show available databases
    ./download_blastdb.py nt                 # download the 'nt' database
    ./download_blastdb.py 16S_ribosomal_RNA  # download a small db
    ./download_blastdb.py nt -j 4 -o blastdb # 4 volumes at once, into ./blastdb
    ./download_blastdb.py nt nr -j 3 --keep
"""

import argparse
import concurrent.futures
import glob
import hashlib
import json
import os
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
    """Thread-safe per-volume progress, timing, and ETA reporting.

    Timing stats use only volumes that did real work (downloaded). Skipped or
    failed volumes are counted toward the [done/total] position but excluded
    from the average and ETA, so a restart full of instant skips doesn't make
    the ETA wildly optimistic.

    ETA is concurrency-aware: with J jobs running, the wall-clock time to finish
    R remaining volumes is roughly R * (avg work per volume) / min(J, R).
    """

    def __init__(self, total, jobs):
        self.total = total
        self.jobs = max(1, jobs)
        self.lock = threading.Lock()
        self.done = 0
        self.work_count = 0
        self.work_time = 0.0
        self.xfer_bytes = 0
        self.xfer_time = 0.0
        self.t0 = time.monotonic()

    def complete(self, basename, status, duration, counts_as_work,
                 xfer_bytes=0, xfer_seconds=0.0):
        with self.lock:
            self.done += 1
            n = self.done
            if counts_as_work:
                self.work_count += 1
                self.work_time += duration
                self.xfer_bytes += xfer_bytes
                self.xfer_time += xfer_seconds
            remaining = self.total - n
            avg = (self.work_time / self.work_count) if self.work_count else None
            avg_speed = mb_per_sec(self.xfer_bytes, self.xfer_time) if self.xfer_time else None

        if counts_as_work:
            head = fmt_dur(duration)
            if xfer_seconds > 0 and xfer_bytes > 0:
                head += f" @ {mb_per_sec(xfer_bytes, xfer_seconds):.1f} MB/s"
            parts = [head]
        else:
            parts = ["skipped"]

        if avg is not None:
            avg_part = f"avg {fmt_dur(avg)}/vol"
            if avg_speed:
                avg_part += f", {avg_speed:.1f} MB/s"
            parts.append(avg_part)
            if remaining > 0:
                conc = min(self.jobs, remaining)
                parts.append(f"ETA {fmt_dur(remaining * avg / conc)}")
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


def print_database_list(records):
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


def process_volume(archive_url, outdir, keep, extract, verbose, tracker):
    """Fetch, verify and (optionally) extract a single volume. Thread-safe.

    Returns (archive_url, ok: bool, message: str|None).
    """
    basename = archive_url.rsplit("/", 1)[-1]
    stem = basename[:-7] if basename.endswith(".tar.gz") else basename
    archive_path = os.path.join(outdir, basename)
    t_start = time.monotonic()
    prefetched_md5 = _FETCH_MD5

    try:
        # Restart currency check: in delete-after-extract mode a finished volume
        # has no archive left but its extracted files are present. Only skip it
        # if our stored .md5 marker matches NCBI's CURRENT .md5 (i.e. it's still
        # the latest version); otherwise fall through and re-download.
        if extract and not keep and not os.path.exists(archive_path) \
                and already_extracted(stem, outdir):
            current = fetch_expected_md5(archive_url)
            stored = read_stored_md5(archive_path)

            if current is None:
                # Server has no checksum to compare against; trust presence.
                tracker.complete(
                    basename, "already extracted; skipping (no upstream checksum)",
                    time.monotonic() - t_start, counts_as_work=False)
                return (archive_url, True, None)

            if stored == current:
                tracker.complete(
                    basename, "already extracted, up to date; skipping",
                    time.monotonic() - t_start, counts_as_work=False)
                return (archive_url, True, None)

            # Stale (marker differs) or unverifiable (no marker) -> refresh,
            # reusing the checksum we just fetched.
            if verbose:
                reason = ("upstream version changed" if stored is not None
                          else "no local version marker")
                log(f"  i {basename}: {reason}; re-downloading to refresh")
            prefetched_md5 = current

        path, xfer_seconds, xfer_bytes = download_and_verify_archive(
            archive_url, outdir, verbose=verbose, expected=prefetched_md5
        )

        if extract:
            extract_archive(path, outdir)
            if not keep:
                try:
                    os.remove(path)
                except OSError:
                    pass
            status = "ok (verified, extracted)"
        else:
            status = "ok (verified)"

        tracker.complete(basename, status, time.monotonic() - t_start,
                         counts_as_work=True, xfer_bytes=xfer_bytes,
                         xfer_seconds=xfer_seconds)
        return (archive_url, True, None)

    except DownloadError as exc:
        # A fast failure shouldn't drag the average down -> not counted as work.
        tracker.complete(basename, f"FAILED -- {exc}",
                         time.monotonic() - t_start, counts_as_work=False)
        return (archive_url, False, str(exc))


def download_databases(names, outdir, keep_archives, extract, jobs):
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
    verbose = jobs == 1  # show wget's live progress bar only when sequential

    log(
        f"{timestamp()}  starting: {total} volume(s) across "
        f"{len(volumes_by_db)} database(s); downloading {jobs} at a time "
        f"into '{outdir}'\n"
    )

    tracker = ProgressTracker(total, jobs)
    results = []  # (archive_url, db_name, ok, message)

    if jobs == 1:
        for url, name in tasks:
            log(f"=== {name}: {url.rsplit('/', 1)[-1]} ===")
            url_, ok, msg = process_volume(
                url, outdir, keep_archives, extract, verbose, tracker
            )
            results.append((url_, name, ok, msg))
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
            future_to_db = {
                pool.submit(
                    process_volume, url, outdir, keep_archives, extract,
                    verbose, tracker,
                ): name
                for url, name in tasks
            }
            for fut in concurrent.futures.as_completed(future_to_db):
                name = future_to_db[fut]
                url_, ok, msg = fut.result()
                results.append((url_, name, ok, msg))

    # Per-database summary: a database is done only if all its volumes succeeded.
    bad_volumes = {}
    for _url, name, ok, _msg in results:
        if not ok:
            bad_volumes.setdefault(name, 0)
            bad_volumes[name] += 1

    log("\n--- summary ---")
    for name in names:
        if name in failures:
            log(f"  {name}: NOT FOUND")
        elif name in bad_volumes:
            log(f"  {name}: FAILED ({bad_volumes[name]}/{volumes_by_db[name]} volume(s) bad)")
        else:
            log(f"  {name}: ok ({volumes_by_db[name]} volume(s))")

    elapsed = tracker.elapsed()
    if tracker.work_count:
        avg = tracker.work_time / tracker.work_count
        gb = tracker.xfer_bytes / 1e9
        per_stream = mb_per_sec(tracker.xfer_bytes, tracker.xfer_time)
        aggregate = mb_per_sec(tracker.xfer_bytes, elapsed)
        log(
            f"  ---\n  {tracker.work_count} volume(s), {gb:.2f} GB downloaded "
            f"in {fmt_dur(elapsed)} (avg {fmt_dur(avg)}/volume at -j {jobs})"
        )
        if tracker.xfer_time > 0:
            speed_line = f"  average speed: {per_stream:.1f} MB/s per stream"
            if jobs > 1:
                speed_line += f"; {aggregate:.1f} MB/s effective (all streams)"
            log(speed_line)
    else:
        log(f"  ---\n  no new volumes downloaded (all already present); "
            f"elapsed {fmt_dur(elapsed)}")
    log(f"  finished {timestamp()}")

    if failures or bad_volumes:
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

    if args.list or not args.databases:
        print_database_list(fetch_database_list())
        return

    download_databases(
        names=args.databases,
        outdir=args.outdir,
        keep_archives=args.keep,
        extract=not args.no_extract,
        jobs=args.jobs,
    )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        die("interrupted", code=130)
