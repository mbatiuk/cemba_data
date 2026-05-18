import re
import logging
from pathlib import Path

log = logging.getLogger(__name__)


def link_fastq(in_fq_dir, out_fq_dir,
               plate_pattern,
               read_type_pattern=r'(R[12])',
               multiplex_group_pattern=None,
               lane_pattern=None,
               lane=None,
               recursive=False):
    """
    Standardize FASTQ file names using symlinks.
    Target naming: {plate}-{multiplex_group}_{lane}_{read_type}.fastq.gz

    Lane resolution priority (highest to lowest):
      1. --lane  (manual override, applied to all files in in_fq_dir)
      2. --lane_pattern  (extracted from each filename via regex)
      3. L001  (hardcoded fallback when neither is provided)

    Multiplex group fallback: '1' (when --multiplex_group_pattern is not provided)
    """
    in_fq_dir = Path(in_fq_dir).absolute()
    out_fq_dir = Path(out_fq_dir).absolute()
    out_fq_dir.mkdir(parents=True, exist_ok=True)

    try:
        plate_re = re.compile(plate_pattern)
        read_type_re = re.compile(read_type_pattern)
        group_re = re.compile(multiplex_group_pattern) if multiplex_group_pattern else None
        lane_re = re.compile(lane_pattern) if lane_pattern else None
    except re.error as e:
        log.error(f"Regex compilation error: {e}")
        return

    if lane and lane_pattern:
        log.warning("Both --lane and --lane_pattern provided. --lane takes precedence; --lane_pattern will be ignored.")

    search_pattern = '**/*.f*q.gz' if recursive else '*.f*q.gz'
    all_files = list(in_fq_dir.glob(search_pattern))
    log.info(f"Found {len(all_files)} FASTQ files in {in_fq_dir}")

    # (plate, multiplex_group, lane) -> {read_type: path}
    groups = {}

    for fpath in all_files:
        fname = fpath.name

        # 1. Plate
        m = plate_re.search(fname)
        if not m:
            log.warning(f"No plate match in '{fname}', skipping.")
            continue
        try:
            plate = m.group(1)
        except IndexError:
            plate = m.group(0)

        # 2. Read type
        m = read_type_re.search(fname)
        if not m:
            log.warning(f"No read type match in '{fname}', skipping.")
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
            continue

        # 3. Multiplex group
        multiplex_group = '1'
        if group_re:
            m = group_re.search(fname)
            if m:
                try:
                    multiplex_group = m.group(1)
                except IndexError:
                    multiplex_group = m.group(0)

        # 4. Lane: manual > pattern > L001
        if lane:
            resolved_lane = lane
        elif lane_re:
            m = lane_re.search(fname)
            if m:
                try:
                    resolved_lane = m.group(1)
                except IndexError:
                    # No capturing group — use the whole match as the lane label
                    resolved_lane = m.group(0)
            else:
                resolved_lane = 'L001'
        else:
            resolved_lane = 'L001'

        key = (plate, multiplex_group, resolved_lane)
        if key not in groups:
            groups[key] = {}

        if read_type in groups[key]:
            log.warning(f"Collision for {key} read {read_type}: "
                        f"skipping '{fname}' (already have '{groups[key][read_type].name}')")
            continue

        groups[key][read_type] = fpath

    # Warn about unpaired files
    for (plate, group, resolved_lane), files in groups.items():
        if len(files) < 2:
            log.warning(f"No mate for {plate}-{group}_{resolved_lane}: only {list(files.keys())} found.")

    # Create symlinks for all files, paired or not
    created = 0
    conflict = 0

    for (plate, group, resolved_lane), files in groups.items():
        for rt, src in files.items():
            dest_name = f"{plate}-{group}_{resolved_lane}_{rt}.fastq.gz"
            dest = out_fq_dir / dest_name

            if dest.exists() or dest.is_symlink():
                if dest.is_symlink() and dest.resolve() == src.resolve():
                    log.debug(f"Symlink already correct, skipping: {dest_name}")
                else:
                    log.warning(f"'{dest_name}' exists and points elsewhere, skipping.")
                    conflict += 1
                continue

            try:
                dest.symlink_to(src)
                created += 1
            except Exception as e:
                log.error(f"Failed to create symlink '{dest_name}': {e}")

    log.info(f"Done: {created} links created, {conflict} skipped (conflict).")
