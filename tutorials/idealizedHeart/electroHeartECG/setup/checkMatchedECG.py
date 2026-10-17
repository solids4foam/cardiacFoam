#!/usr/bin/env python3
"""Check common-source recordings and report their finite-grid differences."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def series(path):
    lines = path.read_text().splitlines()
    header = next(line[1:].split() for line in lines
                  if line.startswith('#') and 'time' in line.split())
    rows = [list(map(float, line.split())) for line in lines
            if line.strip() and not line.startswith('#')]
    assert rows and all(len(row) == len(header) for row in rows), path
    assert all(math.isfinite(v) for row in rows for v in row), path
    assert all(a[0] < b[0] for a, b in zip(rows, rows[1:])), path
    return header, rows


def hashes(case):
    return {str(p.relative_to(case)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in case.glob('*/Vm') if p.parent.name.replace('.', '', 1).isdigit()}


def common_row(rows, time):
    result = min(rows, key=lambda row: abs(row[0]-time))
    assert abs(result[0]-time) < 1e-10, (time, result[0])
    return result


def check(case, online):
    output = case/'postProcessing'
    h, rows = series(output/'matchedLeadFieldECG.dat')
    th, tr = series(output/'matchedTorsoPotentials.dat')
    assert len(tr) == len(rows)
    assert all(r[h.index('sourceUnchanged')] == 1 for r in rows)
    assert all(r[h.index('correctors')] >= 2 for r in rows)
    assert all(0 <= r[h.index('linearResidual')] <= 1e-11 for r in rows)
    saved = json.loads((output/'matchedVmInputHashes.json').read_text())
    assert saved['sha256'] == hashes(case), 'Stored source files changed'
    expected_times = sorted(float(Path(name).parent.name) for name in saved['sha256']
                            if float(Path(name).parent.name) > 0)
    assert len(rows) == len(expected_times)
    assert all(abs(row[0]-time) < 1e-10 for row, time in zip(rows, expected_times))
    if online:
        oh, records = series(output/'leadFieldECG.dat')
    result = {'samples': len(rows), 'source_hash_capture_phase': saved['phase'],
              'source_files_unchanged': True, 'channels': {},
              'scope': 'Finite-grid comparison, not exact-reference error or physiological validation.'}
    total_error = total_reference = 0.
    for lead in ['V'+str(i) for i in range(1, 7)]:
        square = reference = maximum = online_difference = 0.
        for row in rows:
            sampled = common_row(tr, row[0])
            wilson = sum(sampled[th.index(e)] for e in ['RA', 'LA', 'LL'])/3
            forward = sampled[th.index('E'+lead[1:])]-wilson
            assert abs(forward-row[h.index(lead+'_torso')]) < 1e-13
            projected = row[h.index(lead+'_leadField')]
            difference = projected-forward
            assert abs(difference-row[h.index(lead+'_difference')]) < 1e-13
            square += difference*difference
            reference += forward*forward
            maximum = max(maximum, abs(difference))
            if online:
                written = common_row(records, row[0])[oh.index(lead)]
                delta = abs(projected-written)
                assert delta < 1e-12+5e-8*abs(projected), (lead, row[0], delta)
                online_difference = max(online_difference, delta)
        result['channels'][lead] = {
            'RMS_difference_V': math.sqrt(square/len(rows)),
            'RMS_forward_V': math.sqrt(reference/len(rows)),
            'relative_RMS_difference': math.sqrt(square/reference) if reference > 1e-24 else None,
            'maximum_difference_V': maximum,
            'max_online_replay_difference_V': online_difference if online else None}
        total_error += square
        total_reference += reference
    assert total_reference > 1e-24, 'No observable nonzero forward cohort'
    result['cohort_relative_RMS_difference'] = math.sqrt(total_error/total_reference)
    if online:
        ph, pr = series(output/'pseudoECG.dat')
        weighted = ['# time V1 V2 V3 V4 V5 V6']
        for row in rows:
            sample = common_row(pr, row[0])
            wilson = sum(sample[ph.index(e)] for e in ['RA', 'LA', 'LL'])/3
            values = [sample[ph.index('E'+str(i))]-wilson for i in range(1, 7)]
            weighted.append(' '.join(f'{x:.17g}' for x in [row[0]]+values))
        (output/'pseudoWilsonECG.dat').write_text('\n'.join(weighted)+'\n')
    (output/'matchedECGAssessment.json').write_text(json.dumps(result, indent=2)+'\n')
    print('PASS: frozen source, Wilson weighting, finite/converged records'
          + (' and online/replay agreement' if online else ''))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path, nargs='?', default=Path.cwd())
    parser.add_argument('--online', action='store_true')
    parser.add_argument('--capture-source', action='store_true')
    parser.add_argument('--phase', default='before replay')
    args = parser.parse_args()
    if args.capture_source:
        destination = args.case/'postProcessing/matchedVmInputHashes.json'
        captured = hashes(args.case)
        assert captured, 'No stored Vm files'
        destination.parent.mkdir(exist_ok=True)
        if destination.exists():
            assert json.loads(destination.read_text())['sha256'] == captured, 'Stored source files changed'
            print('Reusing the preserved source hash record')
        else:
            destination.write_text(json.dumps({'phase': args.phase, 'sha256': captured}, indent=2)+'\n')
    else:
        check(args.case, args.online)
