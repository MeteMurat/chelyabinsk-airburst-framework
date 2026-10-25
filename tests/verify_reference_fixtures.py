"""Verify preserved fixture identities and downstream arithmetic, NOT MATLAB replay.
Run: python tests/verify_reference_fixtures.py
Python standard library only. No network, writes, or original-archive search.
"""
from __future__ import annotations
import csv, hashlib, json, math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def read(name,base='expected'):
    with (ROOT/base/name).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))

def quantile(x,p):
    x=sorted(x); i=(len(x)-1)*p; lo=math.floor(i); hi=math.ceil(i)
    return x[lo]+(i-lo)*(x[hi]-x[lo])

def verify():
    checks=[]
    def need(ok,name):
        checks.append({'name':name,'passed':bool(ok)})
        if not ok:raise AssertionError(name)
    man=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
    for row in man['files']:
        p=(ROOT/row['path']).resolve()
        need(ROOT in p.parents,'manifest_path:'+row['path'])
        need(hashlib.sha256(p.read_bytes()).hexdigest().upper()==row['sha256'],'sha256:'+row['path'])
    orig=json.loads((ROOT/'docs/source_origins.json').read_text())
    for row in orig:
        if row['transformation']=='NONE_BYTE_IDENTICAL':
            need(hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest().upper()==row['sha256'],
                 'historical_identity:'+row['path'])
    cfg=json.loads((ROOT/'inputs/replay_config.json').read_text())
    design=read('mc_design.csv','inputs'); metrics=read('mc_metrics.csv'); pressure=read('pressure.csv'); pop=read('population.csv')
    need([int(r['run_id']) for r in design]==list(range(1,1001)),'frozen_design_1000_unique_ids')
    need(len({r['fragmentation_seed'] for r in design})==1000,'distinct_fragmentation_seeds')
    need([int(r['run_id']) for r in metrics]==list(range(1,1001)),'metrics_1000_ids')
    need([int(r['run_id']) for r in pressure]==list(range(1,1001)),'pressure_1000_ids')
    need(len(pop)==4000 and len({(r['run_id'],r['threshold_kpa']) for r in pop})==4000,'exposure_4000_unique_id_threshold_pairs')
    cols=['R_3_45_kPa_km','R_1_00_kPa_km','R_0_50_kPa_km','R_0_20_kPa_km']
    radius_error=0.0
    for m,p in zip(metrics,pressure):
        need(m['run_id']==p['run_id'] and m['fragmentation_seed']==p['fragmentation_seed'],'source_pressure_row_identity:'+m['run_id'])
        need(float(m['publish_main_airburst_h_km'])==float(p['h_burst_km']) and
             float(m['publish_E_burst_coupled_kton'])==float(p['E_c_kton']),'source_term_join:'+m['run_id'])
        for t,c in zip(cfg['thresholds_kPa'],cols):
            s=(cfg['pressure_coefficient_kPa_km3_per_kt']*float(p['E_c_kton'])/t)**(1/3)
            radius=math.sqrt(max(s*s-float(p['h_burst_km'])**2,0))
            radius_error=max(radius_error,abs(radius-float(p[c])))
    need(radius_error<1e-9,'all_4000_radius_formula_checks')
    pm={int(p['run_id']):p for p in pressure}
    for p in pop:
        t=float(p['threshold_kpa']); col=cols[cfg['thresholds_kPa'].index(t)]
        need(abs(float(p['radius_km'])-float(pm[int(p['run_id'])][col]))<1e-10,'radius_exposure_join:'+p['run_id']+':'+str(t))
    need(all(float(r['publish_cap_hit_event_count'])==0 for r in read('mc_diagnostics.csv')),'zero_archived_cap_hits')
    need(quantile([float(r['publish_main_airburst_h_km']) for r in metrics],.5)==22.75,'archived_altitude_median')
    need(abs(quantile([float(r['publish_E_burst_coupled_kton']) for r in metrics],.5)-330.8671)<.0001,'archived_energy_median')
    return {'status':'PASS_ARCHIVED_FIXTURE_IDENTITIES_AND_ARITHMETIC',
            'check_count':len(checks),'all_checks_pass':True,'pressure_radius_comparisons':4000,
            'max_radius_error_km':radius_error,'matlab_forward_replay_executed':False,
            'full_reproduction_verified':False,
            'interpretation':'Historical fixture consistency, not a new MATLAB forward calculation.'}

if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
