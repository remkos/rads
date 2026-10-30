#!/usr/bin/env bash
#
# Make cycle index out of pass index (read from standard input).

set -euo pipefail
gdate=$(rads-config --gdate)

convdate_sec() {
    local sec="${1:-0}"
    $gdate -u -d "1985-01-01 00:00:00 UTC + $sec seconds" +%y%m%d%H%M%S.%3N
}

declare -A starttime endtime number npass phases passes

while read -r phase cycle pass startsec endsec eqtime eqlon recs; do
    [[ -n "${phase:-}" ]] || continue

    key="$cycle,$pass"

    # Track each distinct pass once.
    if [[ ! -v "starttime[$key]" ]]; then
        passes["$cycle"]+="$pass"$'\n'
    fi

    starttime["$key"]="$startsec"
    endtime["$key"]="$endsec"
    phases["$cycle"]="$phase"

    # Treat record counts as decimal, including leading zeros.
    number["$cycle"]=$(( ${number[$cycle]:-0} + 10#$recs ))
    npass["$cycle"]=$(( ${npass[$cycle]:-0} + 1 ))
done

((${#passes[@]})) || exit 0

while IFS= read -r cycle; do
    # Preserve the Perl script's lexicographic pass ordering.
    mapfile -t sorted_passes < <(
        printf '%s' "${passes[$cycle]}" | LC_ALL=C sort
    )

    firstpass="${sorted_passes[0]}"
    lastpass="${sorted_passes[${#sorted_passes[@]}-1]}"

    startdate=$(convdate_sec "${starttime[$cycle,$firstpass]}")
    enddate=$(convdate_sec "${endtime[$cycle,$lastpass]}")

    printf '%s %s %s %s %s %s %4i %7i\n' \
        "${phases[$cycle]}" "$cycle" "$firstpass" "$lastpass" \
        "$startdate" "$enddate" "${npass[$cycle]}" "${number[$cycle]}"
done < <(printf '%s\n' "${!passes[@]}" | sort -n)
