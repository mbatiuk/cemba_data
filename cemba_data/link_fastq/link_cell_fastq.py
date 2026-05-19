import re
import logging
from pathlib import Path

log = logging.getLogger(__name__)


def link_cell_fastq(in_fq_dir, out_fq_dir,
                    plate_pattern,
                    multiplex_group_pattern,
                    well_pattern,
                    read_type_pattern=r'(R[12])',
                    recursive=False):
    """
    Create symlinks for already-demultiplexed cell-level FASTQ files, renaming to
    the canonical format expected by mapping_cell_fastq:

        {plate}-{multiplex_group}-{well}-R1.fq.gz
        {plate}-{multiplex_group}-{well}-R2.fq.gz

    All four fields are extracted via user-supplied regexes.

    Parameters
    ----------
    in_fq_dir : str or Path
        Directory containing raw cell-level FASTQ files.
    out_fq_dir : str or Path
        Output directory for renamed symlinks.
    plate_pattern : str
        Regex with one capture group matching the plate ID.
    multiplex_group_pattern : str
        Regex with one capture group matching the multiplex group.
    well_pattern : str
        Regex with one capture group matching the well ID.
    read_type_pattern : str
        Regex with one capture group matching read type (R1/R2).
    recursive : bool
        Search for FASTQ files recursively under in_fq_dir.
    """
    in_fq_dir = Path(in_fq_dir).absolute()
    out_fq_dir = Path(out_fq_dir).absolute()
    out_fq_dir.mkdir(parents=True, exist_ok=True)

    try:
        plate_re = re.compile(plate_pattern)
        group_re = re.compile(multiplex_group_pattern)
        well_re = re.compile(well_pattern)
        read_type_re = re.compile(read_type_pattern)
    except re.error as e:
        log.error(f"Regex compilation error: {e}")
        return

    search_pattern = '**/*.f*q.gz' if recursive else '*.f*q.gz'
    all_files = list(in_fq_dir.glob(search_pattern))
    log.info(f"Found {len(all_files)} FASTQ files in {in_fq_dir}")

    created = 0
    conflict = 0
    skipped = 0

    for fpath in all_files:
        fname = fpath.name

        # 1. Plate
        m = plate_re.search(fname)
        if not m:
            log.warning(f"No plate match in '{fname}', skipping.")
            skipped += 1
            continue
        try:
            plate = m.group(1)
        except IndexError:
            plate = m.group(0)

        # 2. Multiplex group
        m = group_re.search(fname)
        if not m:
            log.warning(f"No multiplex group match in '{fname}', skipping.")
            skipped += 1
            continue
        try:
            multiplex_group = m.group(1)
        except IndexError:
            multiplex_group = m.group(0)

        # 3. Well
        m = well_re.search(fname)
        if not m:
            log.warning(f"No well match in '{fname}', skipping.")
            skipped += 1
            continue
        try:
            well = m.group(1)
        except IndexError:
            well = m.group(0)

        # 4. Read type
        m = read_type_re.search(fname)
        if not m:
            log.warning(f"No read type match in '{fname}', skipping.")
            skipped += 1
            continue
        try:
            read_type_val = m.group(1)
        except IndexError:
            read_type_val = m.group(0)

        if '1' in read_type_val:
            read_type = 'R1'
        elif '2' in read_type_val:
            read_type = 'R2'
        else:
            log.warning(f"Could not normalize read type '{read_type_val}' in '{fname}', skipping.")
            skipped += 1
            continue

        dest_name = f"{plate}-{multiplex_group}-{well}-{read_type}.fq.gz"
        dest = out_fq_dir / dest_name

        if dest.exists() or dest.is_symlink():
            if dest.is_symlink() and dest.resolve() == fpath.resolve():
                log.debug(f"Symlink already correct, skipping: {dest_name}")
            else:
                log.warning(f"'{dest_name}' exists and points elsewhere, skipping.")
                conflict += 1
            continue

        try:
            dest.symlink_to(fpath)
            created += 1
        except Exception as e:
            log.error(f"Failed to create symlink '{dest_name}': {e}")

    log.info(f"Done: {created} links created, {conflict} conflicts, {skipped} skipped (no match).")
