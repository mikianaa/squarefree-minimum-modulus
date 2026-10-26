"""Batch price all candidate projections, then independently audit them.

Build scores are provisional until the separate int64 core audit and final
transport/completeness audit pass. Candidate LP objectives are never trusted.
"""
import argparse,contextlib,io,json,time
from fractions import Fraction
from pathlib import Path
from batch_core2310_power import build
from audit_batch_core2310_power import audit as audit_batch


def read(path):return json.loads(Path(path).read_text())


def run(mode,shard=0,shards=1):
    first=read('data/core2310-common-tail-basis.json');second=read('data/audit-core2310-common-tail-basis.json')
    G=[max(Fraction(int(a),int(first['scale'])),Fraction(int(b),int(second['scale']))) for a,b in zip(first['upper_integers'],second['upper_integers'])]
    allcores=[];scores=[];started=time.time()
    for batch in range(30):
        if mode=='audit' and batch%shards!=shard:continue
        source=f'data/core2310-allproj-batch-{batch}-input.json';manifest=f'data/core2310-allproj-batch-{batch}-manifest.json'
        if mode=='build':
            while True:
                try:rows=read('data/core2310-all-projections-summary.json')
                except (FileNotFoundError,json.JSONDecodeError):time.sleep(1);continue
                if len(rows)>=(batch+1)*60:break
                time.sleep(1)
            subset=rows[batch*60:(batch+1)*60];Path(source).write_text(json.dumps(subset,indent=2))
            with contextlib.redirect_stdout(io.StringIO()):build(source,manifest)
            cores=read(manifest)
            for row,entry in zip(subset,cores):
                core=read(entry['core_file']);w=core['weights_numerators'];den=core['weights_denominator'];free=core['initial_free_moduli']
                extra=Fraction(sum(max(sum(w[x] for x in range(r,2310,d)) for r in range(d)) for d in free),den)
                total=extra+sum(Fraction(v,core['power_seed_denominator'])*g for v,g in zip(core['power_seed_numerators'],G))
                scores.append(dict(projections=row['projections'],core_file=entry['core_file'],total_upper=str(total),value=float(total),
                                   scope='Exact primary core score; independent core and transport audits pending'))
            allcores.extend(cores);Path('data/core2310-allproj-manifest.json').write_text(json.dumps(allcores,indent=2))
            Path('data/core2310-allproj-primary-scores.json').write_text(json.dumps(scores,indent=2))
            print('PRIMARY BATCH',batch,'cases',len(scores),'worst',max(r['value'] for r in scores),'failures',sum(r['value']>=1 for r in scores),'seconds',round(time.time()-started,1),flush=True)
        else:
            with contextlib.redirect_stdout(io.StringIO()):audit_batch(manifest)
            print('INDEPENDENT BATCH',batch,'seconds',round(time.time()-started,1),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['build','audit']);ap.add_argument('--shard',type=int,default=0);ap.add_argument('--shards',type=int,default=1);a=ap.parse_args();run(a.mode,a.shard,a.shards)
