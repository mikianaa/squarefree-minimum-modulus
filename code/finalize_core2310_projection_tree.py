"""Assemble and audit the complete adaptive projection tree for one core."""
import argparse,contextlib,io,json
from pathlib import Path
from audit_batch_core2310_power import audit as audit_batch
from core2310_density30_transport import audit as audit_transport


def read(path):return json.loads(Path(path).read_text())


def prepare():
    state=read('data/core2310-refinement-state.json');base=read('data/core2310-allproj-manifest.json');branch=read('data/core2310-branch165-manifest.json');rows=read('data/core2310-branch165-summary.json')
    replacement=read('data/core2310-last-branch-manifest.json')[0]
    for i,row in enumerate(rows):
        if row['projections']==[0,1,1,4,6] and row['residue165']==13:branch[i]=replacement
    final=[dict(core_file=r['core_file']) for r in state['accepted']]+branch
    assert len(state['accepted'])==1783 and len(branch)==255 and len(final)==2038
    Path('data/core2310-complete-projection-tree.json').write_text(json.dumps(final,indent=2))
    already={r['core_file'] for r in base};extra=[r for r in final if r['core_file'] not in already]
    manifests=[]
    for lo in range(0,len(extra),60):
        name=f'data/core2310-extra-audit-{lo//60}.json';Path(name).write_text(json.dumps(extra[lo:lo+60],indent=2));manifests.append(name)
    Path('data/core2310-extra-audit-manifests.json').write_text(json.dumps(manifests,indent=2));print('LEAVES',len(final),'EXTRA AUDITS',len(extra),'BATCHES',len(manifests),flush=True)


def extra(shard,shards):
    for i,path in enumerate(read('data/core2310-extra-audit-manifests.json')):
        if i%shards!=shard:continue
        with contextlib.redirect_stdout(io.StringIO()):audit_batch(path)
        print('EXTRA INTEGER BATCH',i,'VERIFIED',flush=True)


def final():
    audit_transport([r['core_file'] for r in read('data/core2310-complete-projection-tree.json')],
                    released=(22,33,55,66,77,110,165),adaptive165=True,
                    output='data/audit-core2310-complete-projection-tree.json')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['prepare','extra','audit']);ap.add_argument('--shard',type=int,default=0);ap.add_argument('--shards',type=int,default=1);a=ap.parse_args()
    {'prepare':prepare,'extra':lambda:extra(a.shard,a.shards),'audit':final}[a.mode]()
