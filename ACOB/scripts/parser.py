"""
API to store the anatomical properties of an embryo in an .xlsx file from an .xml file.
.xml-specific global variables are described in the parser_settings util script.
_______________

    + data is structured around time points and named cells only instead of cell snapshots
    + cell count, contact surfaces and volume properties are computed
    + embryo-, tissue- and cell-resolution properties are generated
_______________

@ ASLOUDJ Yanis
07/06/2022
"""


from utils.parser_settings import *

import ast
import pandas as pd
import xml.etree.ElementTree as ET


def get_node_key(node):
    """
    # Description
    ---
    Returns a string usable as a unique node dict key.

    # Argument(s)
    ---
        `node` (ET Element): an ElementTree node.
    """
    if node.tag == cell_snapshot_tag:
        return ast.literal_eval(node.attrib[cell_snapshot_attrib])
    return node.tag

def node_to_dict(node, data):
    """
    # Description
    ---
    Converts an ElementTree node into a dict, and add it to an existing dict.
    Recursively calls itself to process paired (e.g. 'cell_name') and nested (e.g. 'cell_contact_surface') nodes.

    # Argument(s)
    ---
        `node` (ET Element): an ElementTree node.
        `data` (dict): a dict to update with the nodes' data.
    """
    key = get_node_key(node)
    data[key] = {}
    children = node.getchildren()
    if len(children) == 0:
        data[key] = ast.literal_eval(node.text)
    else:
        for child in children:
            data[key] = node_to_dict(child, data[key])
    return data

def xml_to_dict(filename):
    """
    # Description
    ---
    Parses an .xml file into a dict.

    # Argument(s)
    ---
        `filename` (str): an Astec .xml file name.

    # Usage
    ---
    >>> data = xml_to_dict('Astec-pm1.recomputed.xml')
    >>> print(data[name_tag][1330300], len(data[contact_tag][1340307]))
    ... b10.0017* 24
    """
    data = {}
    tree = ET.parse(filename)
    root = tree.getroot()
    for node in root:
        data = node_to_dict(node, data)
    return data

def get_timepoint(cell_snapshot_id):
    """
    # Description
    ---
    Given a cell snapshot id, returns the corresponding timepoint.

    # Argument(s)
    ---
        `cell_snapshot_id` (int): a cell snapshot id.

    # Usage
    ---
    >>> print(get_timepoint(1330300))
    ... 133
    """
    tp = 0
    if len(str(cell_snapshot_id)) > cell_snapshot_k_digits:
        tp = int(str(cell_snapshot_id)[:-cell_snapshot_k_digits])
    return tp

def sort_by_timepoint(data):
    """
    # Description
    ---
    Adds an intermediate level of dict, between the property and the cell snapshot id, corresponding to the timepoint.

    # Argument(s)
    ---
        `data` (dict): data structure generated using xml_to_dict().

    # Usage
    ---
    >>> data = xml_to_dict('Astec-pm3.recomputed.xml')
    >>> print(len(data[contact_tag]))
    ... 13484
    >>> sort_by_timepoint(data)
    >>> print(len(data[contact_tag]))
    ... 90
    """
    for td_prop in time_dependent_properties:
        sorted_prop = {}
        for cell_snapshot_id, value in data[td_prop].items():
            tp = get_timepoint(cell_snapshot_id)
            try:
                sorted_prop[tp][cell_snapshot_id] = value
            except KeyError:
                sorted_prop[tp] = {cell_snapshot_id:value}
        data[td_prop] = sorted_prop

def is_exterior(cell_snapshot_id):
    """
    # Description
    ---
    Tests if a cell snapshot corresponds to the exterior. In this case, its id is 1+tp*10**k.

    # Argument(s)
    ---
        `cell_snapshot_id` (int): a cell snapshot id.

    # Usage
    ---
    >>> print(is_exterior(40001), is_exterior(40012))
    ... True False
    """
    return cell_snapshot_id == 1 + get_timepoint(cell_snapshot_id) * 10 ** cell_snapshot_k_digits

def get_cell_tissue(data, cell_snapshot_id):
    """
    # Description
    ---
    Returns the tissue (or 'exterior') inside which a cell snapshot belongs. It is determined using the tissue fates of the cell.
    If the cell has a unique tissue fate, returns it. If it has more than one fate, returns 'undifferentiated'. 
    Else, returns an integer 0.

    # Argument(s)
    ---
        `data` (dict): data structure generated using xml_to_dict().
        `cell_snapshot_id` (int): a cell snapshot id.

    # Usage
    ---
    >>> data = xml_to_dict('Astec-pm1.recomputed.xml')
    >>> print([get_cell_tissue(data, 40001), get_cell_tissue(data, 10049), get_cell_tissue(data, 1400100)])
    ... ['exterior', 'undifferentiated', 'head epidermis']
    """
    try:
        cell_fate = data[fate_tag][cell_snapshot_id]
        if isinstance(cell_fate, list):
            return 'undifferentiated'
        return cell_fate.lower()
    except KeyError:
        if is_exterior(cell_snapshot_id):
            return 'exterior'
        return 0

def rebrand_cell_name(cell_name):
    """
    # Description
    ---
    Returns a cell name with a human-readable syntax.

    # Argument(s)
    ---
        `cell_name` (str): a .xml cell name.

    # Usage
    ---
    >>> print(rebrand_cell_name('a10.0004*'), rebrand_cell_name('b7.0007_'))
    ... A10.4* B7.7_
    """
    try:
        [lineage, position] = cell_name.split('.')
    except ValueError:
        return cell_name

    [progenitor, generation] = lineage[0], int(lineage[1:])
    [proximity, symmetry] = int(position[:-1]), position[-1]

    mother_generation, mother_proximity = generation, proximity
    while mother_generation > 4:
        if mother_proximity % 2 != 0:
            mother_proximity += 1
        mother_proximity /= 2
        mother_generation -= 1
    # verify if the 8-cell stage progenitor is A4.1, B4.1, a4.2 or b4.2

    if mother_proximity == 1:
        progenitor = progenitor.upper()

    return f'{progenitor}{generation}.{proximity}{symmetry}'

def get_cell_name(data, cell_snapshot_id):
    """
    # Description
    ---
    Returns the cell name (or 'exterior') corresponding to a cell snapshot id. 
    Else, returns the integer cell snapshot id.

    # Argument(s)
    ---
        `data` (dict): data structure generated using xml_to_dict().
        `cell_snapshot_id` (int): a cell snapshot id.

    # Usage
    ---
    >>> data = xml_to_dict('Astec-pm1.recomputed.xml')
    >>> print([get_cell_name(data, 1330300), get_cell_name(data, 1340307), get_cell_name(data, 1340001)])
    ... ['B10.17*', 1340307, 'exterior']
    """
    try:
        return rebrand_cell_name(data[name_tag][cell_snapshot_id])
    except KeyError:
        if is_exterior(cell_snapshot_id):
            return 'exterior'
        return 0

def get_named_dict_and_counter(data, cell_level_dict, name_fun):
    """
    # Description
    ---
    Renames the cell snapshots ids into human readable labels based on an input function (i.e. get_cell_name or get_cell_tissue).
    The individual cell values are aggregated for simultaneous cell snapshots belonging to the same tissue.
    A cell counter (dict) is also returned.

    # Argument(s)
    ---
        `data` (dict): data structure generated using xml_to_dict().
        `cell_level_dict` (dict): dict associating a cell snapshot to a value.
        `name_fun` (fun): function applied to rename the cell snapshots ids.

    # Usage
    ---
    >>> data = xml_to_dict('Astec-pm7.recomputed.xml')
    >>> sort_by_timepoint(data)
    >>> cell_level_contacts = data[contact_tag][5][50002]  # contacts of cell 50002 at tp5
    >>> print(cell_level_contacts)
    ... {50001: 11429.941406, 50003: 14455.183594, ... 50046: 5198.400879}

    >>> tissue_level_contacts, tissue_level_cell_counter = get_named_dict_and_counter(data, cell_level_contacts, get_cell_tissue)
    >>> print(tissue_level_contacts)
    ... {'exterior': 11429.941406, '1st lineage, notochord': 27109.943359999997, ... 'head epidermis': 5198.400879}
    >>> print(tissue_level_cell_counter)
    ... {'exterior': 1, '1st lineage, notochord': 2, ... 'head epidermis': 1}

    >>> named_cell_level_contacts, named_cell_level_cell_counter = get_named_dict_and_counter(data, cell_level_contacts, get_cell_name)
    >>> print(named_cell_level_contacts)
    ... {'exterior': 11429.941406, 'A7.3_': 14455.183594, ... 'a6.6*': 5198.400879}
    >>> print(named_cell_level_cell_counter)
    ... {'exterior': 1, 'A7.3_': 1, ... 'a6.6*': 1}
    """
    named_dict, cell_counter = {}, {}
    for cell_label, value in cell_level_dict.items():
        name = name_fun(data, cell_label)
        named_dict[name] = named_dict.get(name, 0) + value
        cell_counter[name] = cell_counter.get(name, 0) + 1
    return named_dict, cell_counter

def frame_sublevel_global(data, name_fun):
    """
    # Description
    ---
    Returns a dataframe associating at each timepoint an object (i.e. a tissue or a cell) with its volume and its cell count.
    The objects volumes are calculated by summing the volumes of every cell in the object at a given tp.

    # Argument(s)
    ---
        `data` (dict): data structure generated using usage_xml_to_dict() then sort_by_timepoint().
        `name_fun` (fun): function applied to rename the cell snapshots ids.

    # Usage
    >>> data = xml_to_dict('Astec-pm8.recomputed.xml')
    >>> sort_by_timepoint(data)
    >>> tissue_global = frame_sublevel_global(data, get_cell_tissue)
    >>> cell_global = frame_sublevel_global(data, get_cell_name)
    >>> print(tissue_global.head())
    ...   tp                          object    volume  cell_count
        0   1          1st lineage, notochord   3186730           4
        1   1          anterior head endoderm   5836171           6
        2   1                undifferentiated  17041369          20
        3   1  posterior ventral neural plate   2031386           2
        4   1              trunk lateral cell   1851587           2
    >>> print(cell_global.head())
    ...   tp object   volume  cell_count
        0   1  A7.3_   845462           1
        1   1  A7.5_   985306           1
        2   1  A7.2_   962170           1
        3   1  A7.8_  1073303           1
        4   1  A7.4_  1052662           1
    """
    rows = []
    for tp, volume_data in data[volume_tag].items():
        sublevel_volume, cell_counter = get_named_dict_and_counter(data, volume_data, name_fun)
        for object, volume in sublevel_volume.items():
            rows.append({'tp': tp, 'object': object, 'volume': volume, 'cell_count': cell_counter[object]})
    df = pd.DataFrame.from_records(rows)
    df.drop(df[df['object']=='exterior'].index, inplace=True)
    
    # remove cell snapshots ids unpaired to a human-readable biological object :
    return df[ df['object'].apply(lambda x: not str(x).isdigit()) ].reset_index(drop=True)