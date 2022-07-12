"""
Global variables called in the script variability.py.
_______________

    + distances computation methods used
    + evenness indexes used
    + pools of tissues (e.g. germ lines) definition
_______________

@ ASLOUDJ Yanis
07/11/2022
"""


    ### distances computation ###

absolute_pairwise_distance_pool = ['canberra', 'cityblock', 'euclidean', 'hamming', 'matching', 'rogerstanimoto', 'sokalmichener', 'sqeuclidean', 'yule']
relative_pairwise_distance_pool = ['braycurtis', 'correlation', 'cosine', 'dice', 'jaccard', 'jensenshannon', 'sokalsneath']

distance_pool = absolute_pairwise_distance_pool + relative_pairwise_distance_pool

"""
The following ssd.pdist distances methods are ignored:

    + chebyshev: is equal to 1 if two vectors have at least one difference, or else 0
        -> the average distance values computable are not continuous.

    + russelrao & kulsinski: the vectors size is a direct parameter in the calculation
        -> 0 are included in the calculation and artificially increase the computed distance.

    + minkowski & mahalanobis: requires to define a supplementary argument
        -> programmatically invalid.
"""

    ### evenness computation ###

evenness_pool = ['simpson_e', 'heip_e', 'pielou_e']
# mcintosh_e is ignored: -> this index has been little used.


    ### tissues definition ###

ectoderm = ['head epidermis', 'lateral tail epidermis', 'medio-lateral tail epidermis', 'midline tail epidermis', 'posterior ventral neural plate', 'anterior ventral neural plate', 'anterior dorsal neural plate', 'posterior lateral neural plate', 'posterior dorsal neural plate']
mesoderm = ['mesenchyme', '1st lineage, tail muscle', '2nd lineage, tail muscle', 'trunk lateral cell', 'trunk ventral cell', '1st lineage, notochord', '2nd lineage, notochord']
endoderm = ['posterior head endoderm', 'anterior head endoderm', '1st endodermal lineage', '2nd endodermal lineage']
germ_lines = {'ectoderm': ectoderm, 'mesoderm': mesoderm, 'endoderm': endoderm}

epidermis = ['head epidermis', 'lateral tail epidermis', 'medio-lateral tail epidermis', 'midline tail epidermis']
neural_plate = ['posterior ventral neural plate', 'anterior ventral neural plate', 'anterior dorsal neural plate', 'posterior lateral neural plate', 'posterior dorsal neural plate']
ectoderm_lines = {'epidermis': epidermis, 'neural plate': neural_plate}

tissue_subgroups = germ_lines