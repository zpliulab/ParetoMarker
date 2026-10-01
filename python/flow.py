"""Path-induced modules and maximum-flow/minimum-cut candidate selection.

Module construction follows recovered path.py: union bounded simple-path nodes,
then use the induced directed subgraph and virtual target-to-sink edges.
The flow solver implements residual reverse edges, absent from the recovered
path_flow.py; this corrected release does not claim historical equivalence.
"""
import argparse
from collections import deque
import json
from pathlib import Path
import math
import pandas as pd


def maximum_flow(capacities, source, sink, tolerance=1e-12):
    nodes = set(capacities) | {v for edges in capacities.values() for v in edges}
    if source == sink or source not in nodes or sink not in nodes:
        raise ValueError('Distinct source and sink must exist')
    residual = {u: {} for u in nodes}
    for u, edges in capacities.items():
        for v, c in edges.items():
            if not math.isfinite(c) or c < 0 or u == v:
                raise ValueError('Capacities must be finite/nonnegative; no self loops')
            residual[u][v] = c
            residual[v].setdefault(u, 0.)
    value = 0.
    while True:
        parent = {source: None}
        queue = deque([source])
        while queue and sink not in parent:
            u = queue.popleft()
            for v in sorted(residual[u]):
                if residual[u][v] > tolerance and v not in parent:
                    parent[v] = u
                    queue.append(v)
        if sink not in parent:
            reachable = set(parent)
            break
        delta = math.inf
        v = sink
        while v != source:
            u = parent[v]; delta = min(delta, residual[u][v]); v = u
        v = sink
        while v != source:
            u = parent[v]
            residual[u][v] -= delta
            residual[v][u] = residual[v].get(u, 0.) + delta
            v = u
        value += delta
    # For antiparallel original edges retain the nonnegative net-flow form.
    flows = {(u,v): max(0., c-residual[u].get(v,0.))
             for u, edges in capacities.items() for v,c in edges.items()}
    cut = [(u,v,c) for u,edges in capacities.items() for v,c in edges.items()
           if u in reachable and v not in reachable]
    if not math.isclose(value, sum(c for _,_,c in cut), rel_tol=1e-8, abs_tol=1e-9):
        raise ArithmeticError('Flow/cut equality failed')
    return value, flows, cut


def bounded_nodes(graph, source, target, limit):
    """Union of nodes on simple paths with at most limit nodes (not hops)."""
    selected = set()
    def visit(u, path):
        if u == target:
            selected.update(path)
            return
        if len(path) >= limit:
            return
        for v in sorted(graph.get(u, {})):
            if v not in path:
                visit(v, path+[v])
    visit(source, [source])
    return selected


def run(cfg):
    root = Path(cfg['input_dir']); out = Path(cfg['output_dir'])
    if out.exists():
        raise FileExistsError('Choose a fresh output directory')
    limit = cfg['max_path_nodes']; terminal_capacity = float(cfg['sink_capacity'])
    if not isinstance(limit, int) or limit < 2 or not math.isfinite(terminal_capacity) or terminal_capacity <= 0:
        raise ValueError('Invalid path bound or sink capacity')
    edges = pd.read_csv(root/'weighted_network.csv', dtype={'source':str,'target':str})
    if edges[['source','target','capacity']].isna().any().any() or edges.duplicated(['source','target']).any():
        raise ValueError('Missing values or duplicate directed edges')
    graph = {}
    for e in edges.itertuples(index=False):
        c = float(e.capacity)
        if not math.isfinite(c) or c < 0:
            raise ValueError('Invalid edge capacity')
        if e.source != e.target:
            graph.setdefault(e.source, {})[e.target] = c
            graph.setdefault(e.target, {})
    sources = pd.read_csv(root/'sources.csv',dtype=str).gene.dropna().drop_duplicates().tolist()
    targets = pd.read_csv(root/'targets.csv',dtype=str).gene.dropna().drop_duplicates().tolist()
    sink = '__PARETOMARKER_VIRTUAL_SINK__'
    if sink in graph:
        raise ValueError('Reserved virtual sink identifier is present')
    allgenes = set(); modules = []; allflows = []; allcuts = []
    for s in sources:
        nodes = set(); reached = []
        for t in targets:
            if s == t:
                continue
            found = bounded_nodes(graph, s, t, limit)
            if found:
                nodes.update(found); reached.append(t)
        if not reached:
            modules.append(dict(source=s,flow=0.,genes=[]));continue
        cap = {u:{v:c for v,c in graph[u].items() if v in nodes} for u in sorted(nodes)}
        cap[sink] = {}
        for t in reached:
            cap[t][sink] = terminal_capacity
        value, flows, cuts = maximum_flow(cap,s,sink)
        genes = set()
        for (u,v),f in flows.items():
            if f > 1e-12:
                genes.update([u,v]);allflows.append(dict(module_source=s,source=u,target=v,flow=f))
        genes.discard(sink);allgenes.update(genes)
        allcuts.extend(dict(module_source=s,source=u,target=v,capacity=c) for u,v,c in cuts)
        modules.append(dict(source=s,flow=value,genes=sorted(genes)))
    expected_modules = cfg.get('expected_modules')
    if expected_modules is not None:
        if not isinstance(expected_modules, int) or expected_modules < 1:
            raise ValueError('expected_modules must be a positive integer when provided')
        if len(modules) != expected_modules:
            raise ValueError(f'Expected {expected_modules} source-centered modules, found {len(modules)}')
        empty = [m['source'] for m in modules if not m['genes']]
        if empty:
            raise ValueError(f'Expected non-empty initial modules; empty modules for sources: {empty}')
    if len(allgenes)<2:
        raise ValueError('Fewer than two candidate genes; inspect inputs without outcome-driven threshold changes')
    expression = pd.read_csv(root/'expression.csv',dtype={'sample_id':str})
    labels = pd.read_csv(root/'labels.csv',dtype={'sample_id':str})
    if expression.sample_id.duplicated().any() or labels.sample_id.duplicated().any():
        raise ValueError('Duplicate sample IDs')
    genes = sorted(allgenes)
    if set(genes)-set(expression.columns):
        raise ValueError('Selected gene is absent from expression matrix')
    if set(expression.sample_id)-set(labels.sample_id):
        raise ValueError('Missing labels')
    out.mkdir(parents=True); dest=out/'matlab_input';dest.mkdir()
    expression[['sample_id']+genes].to_csv(dest/'expression.csv',index=False)
    labels.set_index('sample_id').loc[expression.sample_id].reset_index().to_csv(dest/'labels.csv',index=False)
    edges.loc[edges.source.isin(genes)&edges.target.isin(genes),['source','target']].to_csv(dest/'network.csv',index=False)

    # Export one binary mask per source-centered flow module.  This is the
    # initial population consumed by the MATLAB evolutionary search.
    # Keeping the masks in the same candidate-gene order makes the Python ->
    # MATLAB boundary explicit and avoids reconstructing modules inside MATLAB.
    initial_rows = []
    for module in modules:
        module_genes = set(module['genes'])
        row = {'module_source': module['source']}
        row.update({g: int(g in module_genes) for g in genes})
        initial_rows.append(row)
    initial = pd.DataFrame(initial_rows, columns=['module_source'] + genes)
    initial.to_csv(dest/'initial_population.csv', index=False)

    pd.DataFrame({'gene':genes}).to_csv(out/'candidate_genes.csv',index=False)
    pd.DataFrame(allflows).to_csv(out/'flow_edges.csv',index=False)
    pd.DataFrame(allcuts).to_csv(out/'cut_edges.csv',index=False)
    (out/'modules.json').write_text(json.dumps(modules,indent=2))
    (out/'config_used.json').write_text(json.dumps(cfg,indent=2))
    print(f'{len(genes)} candidate genes; {len(modules)} source-centered modules; MATLAB inputs: {dest}')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);a=p.parse_args()
    run(json.loads(Path(a.config).read_text()))
