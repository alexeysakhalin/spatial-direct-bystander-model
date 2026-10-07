from array import array
from collections import Counter
import hashlib,math,struct,sys,xml.etree.ElementTree as ET
from pathlib import Path

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def read_cells(xml, registry):
    root = ET.parse(xml).getroot()
    simplified = root.find('.//simplified_data')
    filename = simplified.findtext('filename')
    if Path(filename).name != filename:
        raise ValueError('Cell filename must be a local basename')
    labels = simplified.findall('labels/label')
    width = max(int(label.attrib['index']) + int(label.attrib['size']) for label in labels)
    rows = {}
    for label in labels:
        rows.setdefault(label.text, []).append((int(label.attrib['index']), int(label.attrib['size'])))
    for required in ['ID', 'cell_type', 'dead']:
        if len(rows.get(required, [])) != 1 or rows[required][0][1] != 1:
            raise ValueError('Missing or ambiguous scalar: ' + required)
    with (xml.parent / filename).open('rb') as stream:
        header = stream.read(20)
        if len(header) != 20:
            raise ValueError('Truncated MAT header')
        little = struct.unpack('<5i', header)
        big = struct.unpack('>5i', header)
        if little[0] == 0:
            endian, (_, nrows, ncols, imaginary, name_length) = 'little', little
        elif big[0] == 1000:
            endian, (_, nrows, ncols, imaginary, name_length) = 'big', big
        else:
            raise ValueError('Only native dense real-double v4 cell exports are supported')
        if nrows != width or ncols < 0 or imaginary != 0 or not 1 <= name_length <= 128:
            raise ValueError('Unexpected MAT shape, name length or complex array')
        variable = stream.read(name_length).rstrip(b'\0').decode('ascii')
        if variable != 'cells':
            raise ValueError('Unexpected MAT variable: ' + variable)
        values = array('d')
        values.fromfile(stream, nrows * ncols)
        if sys.byteorder != endian:
            values.byteswap()
        if stream.read(1):
            raise ValueError('Unexpected trailing MAT payload')
    if not all(map(math.isfinite, values)):
        raise ValueError('Nonfinite cell data')
    def scalar(name):
        return values[rows[name][0][0]::nrows]
    ids, types, dead = scalar('ID'), scalar('cell_type'), scalar('dead')
    if len(ids) != len(set(ids)) or any(x != int(x) for x in ids):
        raise ValueError('Invalid or duplicate IDs')
    if any(x not in (0, 1) for x in dead):
        raise ValueError('Invalid death flag')
    names = {int(el.attrib['ID']): el.text for el in simplified.findall('cell_types/type')}
    for x in types:
        if x != int(x) or int(x) not in registry or names[int(x)] != registry[int(x)]['name']:
            raise ValueError('Unknown type or registry-name mismatch')
    living = Counter(int(t) for t, d in zip(types, dead) if d == 0)
    result = {'time_h': float(root.findtext('metadata/current_time')) / 60,
              'live_targets': 0, 'live_effectors': 0, 'retained_dead': int(sum(dead)),
              'live_targets_by_type': {}}
    if root.find('metadata/current_time').attrib['units'] != 'min':
        raise ValueError('Unexpected time units')
    for typ, count in living.items():
        record = registry[typ]
        if record['role'] == 'target':
            result['live_targets'] += count
            result['live_targets_by_type'][record['name']] = count
        elif record['role'] == 'effector':
            result['live_effectors'] += count
        else:
            raise ValueError('Unknown cell role')
    if sum(result[k] for k in ('live_targets', 'live_effectors', 'retained_dead')) != ncols:
        raise ValueError('Cell-count balance failed')
    return result
