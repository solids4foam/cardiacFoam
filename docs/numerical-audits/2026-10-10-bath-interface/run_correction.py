from pathlib import Path
import subprocess,json,time
root=Path(__file__).resolve().parent
for item in json.loads((root/'plan-correction.json').read_text()):
 case=root/'cases'/item['name']
 start=time.monotonic()
 with (case/'audit-run.log').open('w') as log:
  proc=subprocess.run(['/bin/bash','./Allrun',*item['args']],cwd=case,stdout=log,stderr=subprocess.STDOUT)
 data=list((case/'postProcessing').glob('*cells.dat')) if (case/'postProcessing').exists() else []
 print(item['name'],'exit',proc.returncode,'seconds',round(time.monotonic()-start,2),flush=True)
 if proc.returncode:
  for p in [case/'audit-run.log',case/'log.setManufacturedFields',case/'log.cardiacFoam']:
   if p.exists():print(p.name,p.read_text(errors='replace')[-2200:],flush=True)
  break
 rows={}
 for p in data:
  for line in p.read_text().splitlines():
   cols=line.split()
   if len(cols)==4 and not line.startswith('#'):
    try:rows[cols[0]]=list(map(float,cols[1:]))
    except ValueError:pass
 item['errors']=rows
 (case/'audit-summary.json').write_text(json.dumps(item,indent=2))
 print('  '+' '.join(f'{k}={rows[k][1]:.5g}' for k in ('Vm','phiE_heart','phiE_bath','interfaceFlux','gradTangentialError','gradHeartTraceError') if k in rows),flush=True)
