"""Read completed checkpoint totals without examining association results."""
from pathlib import Path
import json
root=Path(__file__).resolve().parents[1]
def main():
    rows=[json.loads(p.read_text()) for p in (root/'results/chunks').glob('cell_*/*.json')]
    cfg=json.loads((root/'protocol/design.json').read_text())
    cal=sum(r['stop']-r['start'] for r in rows if r['phase']=='calibration')
    evaluation=sum(r['stop']-r['start'] for r in rows if r['phase']=='evaluation')
    print(json.dumps(dict(calibration_completed=cal,calibration_planned=4*cfg['calibration_null'],
        evaluation_completed=evaluation,evaluation_planned=4*(cfg['evaluation_null']+8*cfg['evaluation_per_alternative']),
        recorded_method_failures=sum(r['failures'] for r in rows),completed_chunks=len(rows),
        thresholds_frozen=(root/'results/thresholds.json').exists()),indent=2))
if __name__=='__main__':main()
