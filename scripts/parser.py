"""
API to store the anatomical properties of an embryo in an .xlsx file from an .xml file.
.xml-specific global variables are described in the script utils.parser_settings.py.
_______________

    + data is structured around time points and named cells only instead of cell snapshots
    + cell count, contact surfaces and volume properties are computed
    + embryo-, tissue- and cell-resolution properties are generated
_______________

@ ASLOUDJ Yanis
07/07/2022
"""


from scripts.utils.parser_settings import *

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

def get_named_dict_and_counter(data, cell_resolution_dict, name_fun):
    """
    # Description
    ---
    Renames the cell snapshots ids into human readable labels based on an input function (i.e. get_cell_name or get_cell_tissue).
    The individual cell values are aggregated for simultaneous cell snapshots belonging to the same tissue.
    A cell counter (dict) is also returned.

    # Argument(s)
    ---
        `data` (dict): data structure generated using xml_to_dict().
        `cell_resolution_dict` (dict): dict associating a cell snapshot to a value.
        `name_fun` (fun): function applied to rename the cell snapshots ids.

    # Usage
    ---
    >>> data = xml_to_dict('Astec-pm7.recomputed.xml')
    >>> sort_by_timepoint(data)
    >>> cell_resolution_contacts = data[contact_tag][5][50002]  # contacts of cell 50002 at tp5
    >>> print(cell_resolution_contacts)
    ... {50001: 11429.941406, 50003: 14455.183594, ... 50046: 5198.400879}

    >>> tissue_resolution_contacts, tissue_resolution_cell_counter = get_named_dict_and_counter(data, cell_resolution_contacts, get_cell_tissue)
    >>> print(tissue_resolution_contacts)
    ... {'exterior': 11429.941406, '1st lineage, notochord': 27109.943359999997, ... 'head epidermis': 5198.400879}
    >>> print(tissue_resolution_cell_counter)
    ... {'exterior': 1, '1st lineage, notochord': 2, ... 'head epidermis': 1}

    >>> named_cell_resolution_contacts, named_cell_resolution_cell_counter = get_named_dict_and_counter(data, cell_resolution_contacts, get_cell_name)
    >>> print(named_cell_resolution_contacts)
    ... {'exterior': 11429.941406, 'A7.3_': 14455.183594, ... 'a6.6*': 5198.400879}
    >>> print(named_cell_resolution_cell_counter)
    ... {'exterior': 1, 'A7.3_': 1, ... 'a6.6*': 1}
    """
    named_dict, cell_counter = {}, {}
    for cell_label, value in cell_resolution_dict.items():
        name = name_fun(data, cell_label)
        named_dict[name] = named_dict.get(name, 0) + value
        cell_counter[name] = cell_counter.get(name, 0) + 1
    return named_dict, cell_counter

def frame_sub_resolution_global(data, name_fun):
    """
    # Description
    ---
    Returns a dataframe associating at each timepoint an object (i.e. a tissue or a cell) with its volume and its cell count.
    The volumes of the objects are calculated by adding the volumes of every cell in the object at a given tp.

    # Argument(s)
    ---
        `data` (dict): data structure generated using usage_xml_to_dict() then sort_by_timepoint().
        `name_fun` (fun): function applied to rename the cell snapshots ids.

    # Usage
    >>> data = xml_to_dict('Astec-pm8.recomputed.xml')
    >>> sort_by_timepoint(data)
    >>> tissue_global = frame_sub_resolution_global(data, get_cell_tissue)
    >>> cell_global = frame_sub_resolution_global(data, get_cell_name)
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
        sub_resolution_volume, cell_counter = get_named_dict_and_counter(data, volume_data, name_fun)
        for object, volume in sub_resolution_volume.items():
            rows.append({'tp': tp, 'object': object, 'volume': volume, 'cell_count': cell_counter[object]})
    df = pd.DataFrame.from_records(rows)
    df.drop(df[df['object']=='exterior'].index, inplace=True)
    
    # remove cell snapshots ids unpaired to a human-readable biological object :
    return df[ df['object'].apply(lambda x: not str(x).isdigit()) ].reset_index(drop=True)

def merge_dicts(labels_and_dicts):
    """
    # Definition
    ---
    Merges dicts together if they are associated to the same label. If a key is shared by the two dicts, the associated values are summed.
    The output is a new dict where each unique label is a key and the corresponding merged dicts are the values.

    # Argument(s)
    ---
        `label_and_dict` (list of tuples): a tuple with a string and a dict respectively.

    # Usage
    ---
    >>> tuple_A = ('black', {'a': 10, 'b': 3, 'd': 9})
    >>> tuple_B = ('black', {'a': 3, 'b': 4, 'c': 5})
    >>> tuple_C = ('white', {'a': 2, 'b': 3})
    >>> tuple_D = ('white', {'a': 1, 'c': 4})
    >>> print(merge_dicts([tuple_A, tuple_B, tuple_C, tuple_D]))
    ... {'black': {'a': 13, 'b': 7, 'd': 9, 'c': 5}, 'white': {'a': 3, 'b': 3, 'c': 4}}
    """
    output_dict = {}
    for (label, dic) in labels_and_dicts:
        label_entry = output_dict.get(label, {})
        for k,v in dic.items():
            label_entry[k] = label_entry.get(k, 0) + v
        output_dict[label] = label_entry
    return output_dict

def frame_sub_resolution_contacts(data, name_fun):
    """
    # Description
    ---
    Returns a dataframe associating at each timepoint two neighbor objects (or an object and the exterior) and their shared surface.
    A surface shared by two objects (i.e. single cell or tissue) corresponds to the sum of the surfaces shared by every cell composing the object.

    # Argument(s)
    ---
        `data` (dict): data structure generated using xml_to_dict() then sort_by_timepoint().
        `name_fun` (fun): function applied to rename the cell labels.

    # Usage
    ---
    >>> data = xml_to_dict('Astec-pm5.recomputed.xml')
    >>> sort_by_timepoint(data)
    >>> tissue_contacts = frame_sub_resolution_contacts(data, get_cell_tissue)
    >>> print(tissue_contacts.head())
    ...   tp                     neighbor_1                      neighbor_2       surface
        0   0  anterior ventral neural plate                        exterior  65730.128906
        1   0  anterior ventral neural plate                undifferentiated  14517.524414
        2   0  anterior ventral neural plate  posterior ventral neural plate  28305.182251
        3   0  anterior ventral neural plate          1st lineage, notochord   7238.465102
        4   0  anterior ventral neural plate                  head epidermis  20601.085083
    >>> cell_contacts = frame_sub_resolution_contacts(data, get_cell_name)
    >>> print(cell_contacts.head())
    ...   tp neighbor_1 neighbor_2       surface
        0   0      a6.5_   exterior  33600.445312
        1   0      a6.5_      a6.7_   6119.218750
        2   0      a6.5_      A7.8_   2885.224609
        3   0      a6.5_      A7.4_  13091.570312
        4   0      a6.5_      A7.7_   1821.151611
    """
    rows = []
    for tp, cell_snapshot_contacts in data[contact_tag].items():
        object_contacts_at_tp = []

        # associates a cell snapshot to an object and contact surfaces
        for cell_snapshot_id, contacts_data in cell_snapshot_contacts.items():
            object_contacts, cell_counter = get_named_dict_and_counter(data, contacts_data, name_fun)
            object_contacts_at_tp.append((name_fun(data, cell_snapshot_id), object_contacts))

        # aggregates the contact surfaces of identical objects
        object_contacts_at_tp = merge_dicts(object_contacts_at_tp)
        for object_1, object1_contacts in object_contacts_at_tp.items():
            for object_2, surface in object1_contacts.items():
                rows.append({'tp': tp, 'neighbor_1': object_1, 'neighbor_2': object_2, 'surface': surface})

    df = pd.DataFrame.from_records(rows)
    df.drop(df[df['neighbor_1']=='exterior'].index, inplace=True)
    df.drop(df[df['neighbor_1']==df['neighbor_2']].index, inplace=True)

    # remove cell snapshots ids unpaired to a human-readable biological object :
    return df[ df['neighbor_1'].apply(lambda x: not str(x).isdigit()) ].reset_index(drop=True)

def frame_embryo_resolution_global(sub_resolution_global, sub_resolution_contacts):
    """
    # Description
    ---
    Returns a dataframe associating at each timepoint the embryo with its volume, its cell count and its total surface.
    The embryo volume is calculated by summing the volumes of every tissue at a given tp.
    The embryo surface is calculated by summing the surfaces shared by each tissue and the exterior at a given tp.

    # Argument(s)
    ---
        `sub_resolution_global` (pandas df): df associating an object with a volume at each tp.
        `tissue_contacts` (pandas df): df associating an object with its contact surfaces at each tp.

    # Usage
    >>> data = xml_to_dict('Astec-pm1.recomputed.xml')
    >>> sort_by_timepoint(data)
    >>> tissue_global = frame_sub_resolution_global(data, get_cell_tissue)
    >>> tissue_contacts = frame_sub_resolution_contacts(data, get_cell_tissue)
    >>> embryo_global = frame_embryo_resolution_global(tissue_global, tissue_contacts)
    >>> print(embryo_global.head())
    ...   tp  object    volume  total_surface  cell_count
        0   1  embryo  57879877  767962.160155          64
        1   2  embryo  57689015  763912.465820          64
        2   3  embryo  57891788  767863.582646          64
        3   4  embryo  57883863  772252.382447          64
        4   5  embryo  57900541  776648.596679          64
    """
    rows = []
    embryo_surfaces = sub_resolution_contacts[sub_resolution_contacts['neighbor_2'] == 'exterior'].groupby(['tp']).sum()['surface']
    embryo_global = sub_resolution_global.groupby(['tp']).sum()
    for tp, row in embryo_global.iterrows():
        rows.append({'tp': tp, 'object': 'embryo', 'volume': row['volume'], 'total_surface': embryo_surfaces[tp], 'cell_count': row['cell_count']})
    return pd.DataFrame.from_records(rows)

def update_surface_properties(sub_resolution_global, sub_resolution_contacts):
    ###################################################################################################
    # the current implementation of this function is responsible for the long duration of the parsing #
    #   + use a dict approach similar to the other update properties functions instead of a mask ?    #
    ###################################################################################################
    """
    Description
    ---
    Adds a surface property to the global df corresponding to the sum of all the contacts surfaces of the object.
    Afterwards, adds a surface_ratio property to the contacts df. 
    It corresponds to the contact surface between two neighbors divided by the total surface of the first neighbor.

    # Argument(s)
    ---
        `sub_resolution_global` (pandas df): dataframe associating at each tp an object with its volume and cell count.
        `sub_resolution_contacts` (pandas df): dataframe associating at each tp an object with its neighbors and their shared surface.

    # Usage
    ---
    >>> data = xml_to_dict('Astec-pm3.recomputed.xml')
    >>> sort_by_timepoint(data)
    >>> cell_global = frame_sub_resolution_global(data, get_cell_name)
    >>> cell_contacts = frame_sub_resolution_contacts(data, get_cell_name)
    >>> cell_global, cell_contacts = update_surface_properties(cell_global, cell_contacts)
    >>> print(cell_global.head())
    ...   tp object   volume  cell_count  total_surface
        0   0  B6.3*   693710           1   43068.024123
        1   0  B7.2*  1727691           1   89103.627275
        2   0  B7.3*  1499646           1   90973.076203
        3   0  B6.3_   714346           1   42671.960695
        4   0  B7.2_  1676893           1   94243.194306
    >>> print(cell_contacts.head())
    ...   tp neighbor_1 neighbor_2       surface  surface_ratio
        0   0      B6.3*   exterior  19276.703125       0.447587
        1   0      B6.3*      B7.2*   7857.625000       0.182447
        2   0      B6.3*      B6.3_   2284.056152       0.053034
        3   0      B6.3*      B7.2_    146.656235       0.003405
        4   0      B6.3*      B7.8*   8121.777832       0.188580
    """
    df_global, df_contacts = sub_resolution_global.copy(), sub_resolution_contacts.copy()
    df_global['total_surface'] = 0.0
    df_contacts['surface_ratio'] = 0.0
    objects_total_surfaces = df_contacts.groupby(['tp', 'neighbor_1']).sum()['surface']
    
    for (tp, object) in objects_total_surfaces.index:
        object_snapshot_total_surface = objects_total_surfaces[tp][object]
        global_loc_mask = (df_global['tp']==tp) & (df_global['object']==object)
        df_global.loc[global_loc_mask, 'total_surface'] = object_snapshot_total_surface
        contacts_loc_mask = (df_contacts['tp']==tp) & (df_contacts['neighbor_1']==object)
        df_contacts.loc[contacts_loc_mask, 'surface_ratio'] = df_contacts['surface'] / object_snapshot_total_surface

    return df_global, df_contacts

def update_embryo_cell_count(sub_resolution_df, embryo_global):
    """
    # Description
    ---
    Adds an embryo cell count column to a sub-resolution df. It corresponds to the number of cells in the embryo at a specific time point.

    # Argument(s)
    ---
        `sub_resolution_df` (pandas df): a dataframe at the tissue- or cell-resolution.
        `embryo_global` (pandas df): a global dataframe at the embryo-resolution.

    # Usage
    ---
    >>> data = xml_to_dict('Astec-pm1.recomputed.xml')
    >>> sort_by_timepoint(data)
    >>> cell_global = frame_sub_resolution_global(data, get_cell_name)
    >>> cell_contacts = frame_sub_resolution_contacts(data, get_cell_name)
    >>> embryo_global = frame_embryo_resolution_global(cell_global, cell_contacts)
    >>> update_embryo_cell_count(cell_global, embryo_global)
    >>> print(cell_global.head())
    ...   tp  object  volume  cell_count  embryo_cell_count
        0   1   a7.9_  626879           1                 64
        1   1  a7.10_  604934           1                 64
        2   1  a7.13_  660881           1                 64
        3   1  a7.10*  700065           1                 64
        4   1   A7.4_  907942           1                 64
    """
    tp_and_embryo_cell_count_dict = embryo_global.set_index('tp').to_dict()['cell_count']
    sub_resolution_df['embryo_cell_count'] = sub_resolution_df['tp'].apply(lambda x: tp_and_embryo_cell_count_dict[x])

def update_volume_property(sub_resolution_global, embryo_global):
    """
    # Description
    ---
    Adds a volume_ratio property to a global df.
    It corresponds to the volume of the object divided by the volume of the embryo.

    # Argument(s)
    ---
        `sub_resolution_global` (pandas df): a global dataframe at the tissue- or cell-resolution.
        `embryo_global` (pandas df): a global dataframe at the embryo-resolution.

    # Usage
    ---
    >>> data = xml_to_dict('Astec-pm1.recomputed.xml')
    >>> sort_by_timepoint(data)
    >>> cell_global = frame_sub_resolution_global(data, get_cell_name)
    >>> cell_contacts = frame_sub_resolution_contacts(data, get_cell_name)
    >>> embryo_global = frame_embryo_resolution_global(cell_global, cell_contacts)
    >>> update_volume_property(cell_global, embryo_global)
    >>> print(cell_global.head())
    ...    tp  object  volume  cell_count  volume_ratio
        0   1   a7.9_  626879           1      0.010831
        1   1  a7.10_  604934           1      0.010452
        2   1  a7.13_  660881           1      0.011418
        3   1  a7.10*  700065           1      0.012095
        4   1   A7.4_  907942           1      0.015687
    """
    tp_and_embryo_volume_dict = embryo_global.set_index('tp').to_dict()['volume']
    sub_resolution_global['volume_ratio'] = sub_resolution_global['tp'].apply(lambda x: tp_and_embryo_volume_dict[x])
    sub_resolution_global['volume_ratio'] = sub_resolution_global['volume'] / sub_resolution_global['volume_ratio']

def parse_xml_completely(filename, geometry=False):
    """
    # Description
    ---
    Reads a segmentation .xml file and returns a dict of dataframes.
    The dict has 3 entries corresponding each to a resolution : 'embryo', 'tissue' and 'cell'. 
    The tissue and cell entries are nested dicts, with 2 keys each corresponding to a data type : 'global' and 'contacts'.

    # Argument(s)
    ---
        `filename` (str): an Astec .xml file name.
        `geometry` (bool): if True, updates the geometrical properties of surface and volume (time-consuming).

    # Usage
    ---
    >>> astec = parse_xml_completely('Astec-pm3.recomputed.xml')
    >>> print(astec['embryo'].head())
    ...   tp  object    volume  total_surface  cell_count
        0   0  embryo  67400259  906306.334958          46
        1   1  embryo  67789699  896366.383059          47
        2   2  embryo  67878724  916702.315671          53
        3   3  embryo  68163715  932586.093023          60
        4   4  embryo  68172491  907851.260497          62
    >>> print(astec['tissue']['global'].head())
    ...   tp                    object    volume  cell_count  embryo_cell_count
        0   0          undifferentiated  29154599          18                 46
        1   0   posterior head endoderm   3057113           2                 46
        2   0  1st lineage, tail muscle   5910417           4                 46
        3   0    anterior head endoderm   6926467           6                 46
        4   0                mesenchyme   1251309           2                 46
    >>> print(astec['cell']['contacts'].head())
    ...   tp neighbor_1 neighbor_2       surface  embryo_cell_count
        0   0      B6.3*   exterior  19276.703125                 46
        1   0      B6.3*      B7.2*   7857.625000                 46
        2   0      B6.3*      B6.3_   2284.056152                 46
        3   0      B6.3*      B7.2_    146.656235                 46
        4   0      B6.3*      B7.8*   8121.777832                 46
    """
    output = {}
    data = xml_to_dict(filename)
    sort_by_timepoint(data)

    for (resolution, fun) in [('cell', get_cell_name), ('tissue', get_cell_tissue)]:
        output[resolution] = {'global': frame_sub_resolution_global(data, fun), 'contacts': frame_sub_resolution_contacts(data, fun)}
    output['embryo'] = frame_embryo_resolution_global(*output[resolution].values())

    if geometry:
        for resolution in ['tissue', 'cell']:
            [output[resolution]['global'], output[resolution]['contacts']] = update_surface_properties(*output[resolution].values())
            update_volume_property(output[resolution]['global'], output['embryo'])

    for resolution in ['tissue', 'cell']:
        contacts = output[resolution]['contacts']
        # filter out the neighbors unpaired to a cell name, a tissue fate or the exterior (e.g. 11th generation cells) indicated by 0 :
        output[resolution]['contacts'] = contacts[ contacts['neighbor_2'].apply(lambda x: not str(x).isdigit()) ].reset_index(drop=True)
        for data_type in ['global', 'contacts']:
            update_embryo_cell_count(output[resolution][data_type], output['embryo'])

    return output

def save_multi_resolution_data(astec, filename):
    """
    # Definition
    ---
    Saves a parsed .xml file data into a multi-sheets Excel document.

    # Argument(s)
    ---
        `astec` (dict): data structure generated using parse_xml_completely().
        `filename` (str): name of the saved file.

    # Usage
    >>> astec = parse_xml_completely('Astec-pm3.recomputed.xml')
    >>> save_multi_resolution_data(astec, 'Astec-pm3.xlsx')
    >>> data = pd.read_excel('Astec-pm3.xlsx', engine='openpyxl', sheet_name=None)
    >>> print(data['embryo'].head())
    ...   tp  object    volume  total_surface  cell_count
        0   0  embryo  67400259  906306.334958          46
        1   1  embryo  67789699  896366.383059          47
        2   2  embryo  67878724  916702.315671          53
        3   3  embryo  68163715  932586.093023          60
        4   4  embryo  68172491  907851.260497          62
    """
    with pd.ExcelWriter(filename) as writer:
        astec['embryo'].to_excel(writer, sheet_name='embryo', index=False)
        for resolution in ['tissue', 'cell']:
            astec[resolution]['global'].to_excel(writer, sheet_name=resolution, index=False)
            astec[resolution]['contacts'].to_excel(writer, sheet_name=f'{resolution}_contacts', index=False)