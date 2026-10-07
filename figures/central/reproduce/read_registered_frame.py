"""Read registered saved frames, without starting or calibrating a simulator.

Only synthetic agents are read. Output names distinguish counts, sampled death
records and an uncalibrated geometric area proxy from experimental observations.
"""
from pathlib import Path
import argparse
import json
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy.io import loadmat


def read_frame(path, registry_path):
    root = ET.parse(path).getroot()
    nodes = list(root.iter('simplified_data'))
    if len(nodes) != 1:
        raise ValueError('Exactly one cell matrix schema is required')
    node = nodes[0]
    indexed = {}
    for label in node.iter('label'):
        start, size = int(label.get('index')), int(label.get('size', '1'))
        for k in range(size):
            if start+k in indexed:
                raise ValueError('Overlapping row labels')
            indexed[start+k] = label.text if size == 1 else f'{label.text}_{k}'
    if sorted(indexed) != list(range(len(indexed))):
        raise ValueError('Incomplete row schema')
    mat_name = node.findtext('filename')
    if not mat_name or Path(mat_name).name != mat_name:
        raise ValueError('Expected a local matrix basename')
    matrices = {k: v for k, v in loadmat(path.parent / mat_name).items()
                if not k.startswith('__')}
    if len(matrices) != 1:
        raise ValueError('Ambiguous matrix file')
    matrix = next(iter(matrices.values()))
    if matrix.ndim != 2 or matrix.shape[0] != len(indexed):
        raise ValueError('Matrix shape does not match explicit row schema')
    d = pd.DataFrame(matrix.T, columns=[indexed[k] for k in range(len(indexed))])
    for name in ['ID', 'dead', 'cell_type', 'total_volume', 'position_0', 'position_1', 'position_2']:
        if name not in d or not np.isfinite(d[name]).all():
            raise ValueError('Missing or nonfinite required data: '+name)
    if d.ID.duplicated().any() or not d.dead.isin([0, 1]).all():
        raise ValueError('Duplicate IDs or invalid death flags')
    types = {int(n.get('ID')): n.text for n in node.findall('./cell_types/type')}
    registry = json.loads(Path(registry_path).read_text())['states']
    expected = {r['id']: r['name'] for r in registry}
    if types != expected or not set(d.cell_type).issubset(types):
        raise ValueError('Frame type map does not match the explicit target-state registry')
    if len(expected) != len(registry):
        raise ValueError('Duplicate registered type ID')
    clock = root.find('.//current_time')
    if clock is None or clock.get('units') != 'min':
        raise ValueError('Time unit must be minutes')
    return float(clock.text) / 60, d, types


