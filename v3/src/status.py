"""Read progress without changing simulations or inspecting inferential results."""
from pathlib import Path
import json
V=Path(__file__).resolve().parents[1]
totals={'cal':0,'eval':0};chunks=0;seconds=0
for p in (V/'results/chunks').rglob('*.json'):
    d=json.loads(p.read_text());phase=p.stem.split('_')[0]
    totals[phase]+=d['rows'];chunks+=1;seconds+=d['seconds']
print(json.dumps(dict(completed_datasets=totals,target_calibration=12000,target_evaluation=19200,completed_chunks=chunks,worker_seconds=round(seconds,1),thresholds_frozen=(V/'results/thresholds.json').exists()),indent=2))
