"""
Global variables called in the parser.py script.
_______________

    + informations on the structure of the .xml files
_______________

@ ASLOUDJ Yanis
07/06/2022
"""


    ### .xml anatomical properties ###

contact_tag = 'cell_contact_surface'
volume_tag = 'cell_volume'
time_dependent_properties = [contact_tag, volume_tag]
# time_dependent_properties : different cell snapshots of a unique biological cell have different values.

name_tag = 'cell_name'
fate_tag = 'cell_fate'
lineage_tag = 'cell_lineage'
guignard_tag = 'tissuefate_guignard_2020'
time_independent_properties = [name_tag, fate_tag, lineage_tag, guignard_tag]
# time_independent_properties : different cell snapshots of a unique biological cell have the same value.


    ### .xml cell snapshots  ###

cell_snapshot_tag = 'cell'
cell_snapshot_attrib = 'cell-id'
cell_snapshot_k_digits = 4
# cell_snapshot_k_digts : snapshot ids are concatenated from a time point and a unique combination of k digits.


"""
@ Astec-pm8.recomputed.xml :

        <data>

            <cell_lineage>
                <cell cell-id="530208">[540209]</cell>
            </cell_lineage>

            <cell_name>     
                <cell cell-id="530208">'a9.0048_'</cell>
                <cell cell-id="540209">'a9.0048_'</cell>
            </cell_name>

            <cell_volume>
                <cell cell-id="530208">189776</cell>
                <cell cell-id="540209">191351</cell>
            </cell_volume>

        </data>

"""