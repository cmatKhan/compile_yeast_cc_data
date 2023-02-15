# # stdlib
# from ipaddress import summarize_address_range
# import re
# import os
# from typing import Callable, Literal
# import logging
# import warnings
# # local
# from callingcardstools.BarcodeParser import BarcodeParser
# from callingcardstools.Database.DatabaseApi import DatabaseApi
# from .peak_calling import call_peaks_with_background
# # outside
# import pandas as pd
# import numpy as np
# import edlib

# # TODO switch to SQLalchemy and remove this
# warnings.simplefilter(action='ignore', category=UserWarning)

# logging.getLogger(__name__).addHandler(logging.NullHandler())


# class HopsDb(DatabaseApi):
#     """A database manager for yeast calling cards experiments which includes background hops data"""

#     def __init__(self, *args) -> None:
#         """Create a database structured to perform QC on yeast data. positional arguments are the same as DatabaseAPI"""

#         super().__init__(*args)
#         # create batch summary table
#         self.create_batch_table()
#         # create qc tables
#         self.create_qc_tables()
#         # create trigger
#         self.create_qc_table_triggers()

#     def add_regions(self, regions_df: pd.DataFrame, *args, **kwargs) -> bool:
#         """Add a regions table to the database. note that this function is a convenience function 
#         which wraps DatabaseApi.add_table. It does some checking of the regions dataframe and then 
#         calls DatabaseApi.add_table with the regions_df, "regions" set for the table_type, and 'bed3' 
#         set for the table_format. You can pass additional arguments to DatabaseApi.add_table with either 
#         positional or keyword arguments.

#         Args:
#             regions_df (pd.DataFrame): _description_
#             *args (list, optional): additional arguments to DatabaseApi.add_table
#             **kwargs (dict, optional): additional keyword arguments to DatabaseApi.add_table

#         Raises:
#             AttributeError: _description_
#             ValueError: _description_

#         Returns:
#             bool: _description_
#         """
#         if len(np.setdiff1d(self.required_fields['bed3'], regions_df.columns)) > 0:
#             raise AttributeError(
#                 f"regions df must have at least the headers {self.required_fields['bed3']}")
#         elif len(regions_df) > len(regions_df.groupby(self.required_fields['bed3']).size()):
#             raise ValueError(f"regions_df has duplicate entries in the combination of {self.required_fields['bed6']}. "
#                              f"These regions should be unique as calling cards hop counting is unstranded. "
#                              f"Use an annotation table if you need to annotate the same region with multiple affiliations.")
#         else:
#             self.add_frame(regions_df, 'regions', 'bed3', *args, **kwargs)

#     def add_annotation_fk(self, df: pd.DataFrame, regions_fk_table: str) -> bool:
#         """Add an annotation table. Optionally key to a regions table

#         Args:
#             df (pd.DataFrame): A dataframe with little restriction 
#             to fields. However, the assumption is that you will want these 
#             formatted in such a way that they can be used to annotate your 
#             significant regions. Suggested minimum columns: chr,start,stop, which 
#             can be used to join/overlap with the bed6 and qBed tables. Note that 
#             adding a table with the same tablename_suffix will overwrite the current 
#             table with the same suffix if one exists.
#             tablename_suffix (str): suffix to be appended to annotation_ to form 
#             the new tablename. Eg,
#             regions_fk_table (str): Name of a regions table to which to key. The foreign key 
#             will be <fk_tablename>_id

#         Returns:
#             bool: True if successful
#         """
#         # check that the fk table exists
#         if not regions_fk_table in self.list_tables(self.con):
#             raise AttributeError(f"{regions_fk_table} not in the database!")
#         # create the fk_fieldname and ensure that it exists in the database
#         fk_field = regions_fk_table+"_id"
#         if not fk_field in df.columns:
#             raise KeyError(f"{fk_field} is not a field in the input dataframe")
#         # create the tablename
#         tablename = 'annotation'+f'_{regions_fk_table}'
#         # add this table and the fk to the required_fields dict
#         self.required_fields = {tablename: [fk_field]}  # pylint:disable=W0201
#         # add the frame
#         self.add_frame(df,
#                        table_type='annotation',
#                        table_format='annotation',
#                        tablename=tablename,
#                        fk_table=regions_fk_table,
#                        fk_keys=[fk_field])

#     def create_aggregate_view(self, qbed_tbl: str, regions_tbl: str, clean: bool = False) -> bool:
#         """For each background and experiment table, create a view where the hops are 
#         aggregated over the regions provided in the regions table

#         Args:
#             tbl_list (list): list of table names from which aggregate views will 
#             be created. Eg, a list of all experiment and background tables
#             regions_tbl (str): name of the regions table to use to aggregate the hops

#         Returns:
#             bool: True if successful
#         """
#         # verify qbed_tbl
#         if qbed_tbl not in [x for x in self.list_tables(self.con)
#                             if x.startswith('experiment_') or x.startswith('background_')]:
#             AttributeError(f"No table named {qbed_tbl} "
#                            f"in database with prefix experiment_ or background_")

#         # verify regions_tbl
#         # check if region_tbl is in the regions table list
#         if regions_tbl not in [x for x in self.list_tables(self.con)
#                                if x.startswith('regions_')]:
#             AttributeError(f"No table named {regions_tbl} "
#                            f"in the regions tbls in the database")

#         # create the views
#         cur = self.con.cursor()
#         viewname = regions_tbl + '_' + qbed_tbl

#         if clean:
#             drop_view_sql = f'DROP VIEW IF EXISTS "main"."{viewname}";'
#             self.db_execute(cur, drop_view_sql)

#         create_view_sql = \
#             " ".join([f"CREATE VIEW IF NOT EXISTS {viewname} AS",
#                       "SELECT r.chr as chr,r.start as start,r.end as end, count(*) as hops, batch_id, name",
#                       f"FROM {regions_tbl} as r",
#                       f"LEFT JOIN {qbed_tbl} as x",
#                       "WHERE (x.chr = r.chr AND x.start BETWEEN r.start AND r.end)",
#                       "GROUP BY batch_id, name, r.chr,r.start,r.end"])

#         self.db_execute(cur, create_view_sql)
#         self.con.commit()

#         return True

#     def _regions_background_expr_sql(self, regions_tblname: str, background_tblname: str, experiment_tblname: str) -> str:
#         """_summary_

#         Note that the keyword arguments must remain the same as the expected 
#         keyword arguments in peak_caller

#         Args:
#             regions (str): _description_
#             background (str): _description_

#         Raises:
#             AttributeError: _description_

#         Returns:
#             str: _description_
#         """
#         try:
#             hop_views = {
#                 'background': regions_tblname + '_' + background_tblname,
#                 'experiment': regions_tblname + '_' + experiment_tblname
#             }
#         except KeyError as exp:
#             raise 'regions_background_sql requires that ' +\
#                 'regions,background,experiment ' +\
#                 'be passed as a named argument' from exp

#         # check that the views are in the database
#         tbl_list = self.list_tables(self.con)
#         for t in hop_views.values():
#             if t not in tbl_list:
#                 raise AttributeError(f"No view called {t} exists. "
#                                      f"Try calling region_aggregate() and then resubmit")

#         # finally, I decided that the areas where the experiment is null
#         # actually don't matter -- no reason to save them, anyway.
#         # note: revision 20230207, not actually true since it makes it
#         # difficult to compare across experiments
#         sql = " ".join(["SELECT e.chr AS chr, e.start AS start, e.end AS end,",
#                         "IFNULL(b.hops,0) AS bg_hops, e.hops AS expr_hops,",
#                         "e.batch_id as batch_id",
#                         f"FROM {hop_views['experiment']} as e",
#                         f"LEFT JOIN {hop_views['background']} as b USING(chr,start,end,name)",
#                         "ORDER BY batch_id,chr,start,end"])
#         return sql

#     def create_batch_table(self, **kwargs) -> dict:
#         """_summary_

#         Returns:
#             dict: _description_
#         """
#         tables_dict = {
#             'batch': {
#                 'fk': None,
#                 'fields': {
#                     'batch': 'TEXT NOT NULL',
#                     'tf': 'TEXT NOT NULL',
#                     'replicate': 'TEXT DEFAULT \'none\''}},

#         }

#         return {k: self.new_table(k, v['fields'], v['fk'], **kwargs) for k, v in tables_dict.items()}

#     def create_qc_tables(self, **kwargs) -> dict:
#         """_summary_

#         Returns:
#             dict: _description_
#         """
#         # note DatabaseAPI.pk will be added as primary key -- by default id
#         tables_dict = {
#             'qc_manual_review': {
#                 'fk': ['batch'],
#                 'fields': {
#                     'batch_id': 'INTEGER',
#                     'rank_recall': '"rank_recall"	TEXT DEFAULT \'unreviewed\' CHECK("rank_recall" IN (\'pass\', \'fail\', \'unreviewed\', \'note\'))',
#                     'chip_better': 'TEXT DEFAULT \'unreviewed\' CHECK("chip_better" IN (\'yes\', \'no\', \'unreviewed\', \'note\'))',
#                     'data_usable': 'TEXT DEFAULT \'unreviewed\' CHECK("data_usable" IN (\'yes\', \'no\', \'unreviewed\', \'note\'))',
#                     'passing_replicate': 'TEXT DEFAULT \'unreviewed\' CHECK("passing_replicate" IN (\'yes\', \'no\', \'unreviewed\', \'note\'))',
#                     'notes': 'TEXT DEFAULT \'unreviewed\''
#                 }
#             },
#             'qc_r1_to_r2_tf': {
#                 'fk': ['batch'],
#                 'fields': {
#                     'batch_id': 'INTEGER',
#                     'edit_dist': 'INTEGER DEFAULT -1 NOT NULL',
#                     'tally': 'INTEGER DEFAULT -1 NOT NULL'
#                 }},
#             'qc_r2_to_r1_tf': {
#                 'fk': ['batch'],
#                 'fields': {
#                     'batch_id': 'INTEGER',
#                     'edit_dist': 'INTEGER DEFAULT -1 NOT NULL',
#                     'tally': 'INTEGER DEFAULT -1 NOT NULL'
#                 }},
#             'qc_tf_to_transposon': {
#                 'fk': ['batch'],
#                 'fields': {
#                     'batch_id': 'INTEGER',
#                     'edit_dist': 'INTEGER DEFAULT -1 NOT NULL',
#                     'tally': 'INTEGER DEFAULT -1 NOT NULL'
#                 }},
#             'qc_alignment': {
#                 'fk': ['batch'],
#                 'fields': {
#                     'batch_id': 'INTEGER',
#                     'total': 'INTEGER DEFAULT -1 NOT NULL',
#                     'mapped': 'INTEGER DEFAULT -1 NOT NULL',
#                     'multimap': 'INTEGER DEFAULT -1 NOT NULL',
#                 }},
#             'qc_hops': {
#                 'fk': ['batch'],
#                 'fields': {
#                     'batch_id': 'INTEGER',
#                     'total': 'INTEGER DEFAULT -1 NOT NULL',
#                     'transpositions': 'INTEGER DEFAULT -1 NOT NULL',
#                     'plasmid_transpositions': 'INTEGER DEFAULT -1 NOT NULL'
#                 }}
#         }

#         return {k: self.new_table(k, v['fields'], v['fk'], **kwargs) for k, v in tables_dict.items()}

#     def create_qc_table_triggers(self) -> None:

#         # TODO address hard coding here -- turn into property?
#         batch_tablename = 'batch'
#         batchtable_fk = 'batch_id'

#         repeat_record_qc_tbls = ['qc_r1_to_r2_tf',
#                                  'qc_r2_to_r1_tf', 'qc_tf_to_transposon']

#         repeat_sql_prefix = f"CREATE TRIGGER IF NOT EXISTS repeat_qc_records AFTER INSERT ON {batch_tablename} BEGIN "
#         insert_stmts = []
#         for tbl in repeat_record_qc_tbls:
#             insert_stmts.append("; ".join([f'INSERT INTO {tbl} ({batchtable_fk},edit_dist) VALUES (new.id, {x})'
#                                            for x in range(5)]))
#         repeat_entry_sql = repeat_sql_prefix + \
#             ";".join(insert_stmts) + '; END;'

#         cur = self.con.cursor()
#         self.db_execute(cur, repeat_entry_sql)
#         self.con.commit()

#         single_entry_qc_tbls = ['qc_manual_review', 'qc_alignment', 'qc_hops']
#         single_entry_sql = f"CREATE TRIGGER IF NOT EXISTS single_qc_records AFTER INSERT ON {batch_tablename} BEGIN " +\
#             "; ".join([f'INSERT INTO {x} ({batchtable_fk}) VALUES (new.id)'
#                        for x in single_entry_qc_tbls]) +\
#             '; END;'

#         cur = self.con.cursor()
#         self.db_execute(cur, single_entry_sql)
#         self.con.commit()

#     def add_batch_qc(self, barcode_details: BarcodeParser, id_to_bc_map_path: str = None, split_char="x") -> int:
#         """_summary_

#         Args:
#             barcode_details (BarcodeDetails): _description_
#             id_to_bc_map_path (str, optional): _description_. Defaults to None.
#             split_char (str, optional): _description_. Defaults to "x".

#         Returns:
#             int: _description_
#         """
#         # extract data from the barcode_details
#         tfs = list(
#             barcode_details.barcode_dict['components']['tf']['map'].values())
#         batch = barcode_details.barcode_dict['batch']

#         parsed_tfs = []
#         replicates = []
#         for x in tfs:
#             try:
#                 tf = x.split(split_char)[0].replace(')', '')
#                 rep = x.split(split_char)[1].replace(')', '')
#             except IndexError:
#                 tf = x
#                 rep = 'none'
#             parsed_tfs.append(tf)
#             replicates.append(rep)

#         parsed_tfs.append('undetermined')
#         replicates.append('none')

#         # create df with nrow == len(tfs). batch is repeated for each record.
#         # records is either split from the tf name, or 'none'
#         df = pd.DataFrame(
#             {'batch': [batch for x in parsed_tfs],
#              'tf': parsed_tfs,
#              'replicate': replicates})

#         # by row, check if the entry already exists. if it does, then insert.
#         # else, skip
#         add_id_map = False
#         insert_sql = "INSERT INTO batch (batch,tf,replicate) VALUES ('%s','%s','%s')"
#         for idx, row in df.iterrows():
#             try:
#                 self.get_batch_id(row['batch'], row['tf'], row['replicate'])
#             except IndexError:
#                 add_id_map = True
#                 cur = self.con.cursor()
#                 self.db_execute(cur, insert_sql %
#                                 (row['batch'], row['tf'], row['replicate']))
#                 self.con.commit()

#         # note that a trigger creates which creates default entries in all other
#         # qc tables

#         if id_to_bc_map_path and add_id_map:
#             self.add_read_qc(barcode_details, id_to_bc_map_path, split_char)

#     def _update_qc_table(self, tablename: str, id_col: str, id_value: str, update_dict: dict) -> None:
#         """_summary_

#         Args:
#             tablename (str): _description_
#             id_col (str): _description_
#             id_value (str): _description_
#             update_dict (dict): _description_
#         """
#         for record in update_dict:
#             sql = f"UPDATE {tablename} SET tally = {record.get('tally')} "\
#                 f"WHERE {id_col} = {id_value} AND "\
#                 f"edit_dist = {record.get('edit_dist')}"
#             cur = self.con.cursor()
#             self.db_execute(cur, sql)
#             self.con.commit()

#     def _summarize_r1_primer(self, id_to_bc_df: pd.DataFrame, r1_primer: str, r2_transposon: str) -> pd.DataFrame:
#         """group by a given r1_primer and return the number of r2_transposon seq within 0 to 4 edit distance 
#         of the expected r2_transposon seq

#         Args:
#             id_to_bc_df (pd.DataFrame): _description_
#             r1_primer (str): _description_
#             r2_transposon (str): _description_

#         Returns:
#             pd.DataFrame: dataframe with 5 rows
#         """

#         primer_df = id_to_bc_df[id_to_bc_df.r1_primer == r1_primer]
#         summarised_dict = {}
#         if not primer_df.shape[0] == 0:
#             primer_df = primer_df\
#                 .assign(edit_dist=primer_df
#                         .apply(lambda x:
#                                edlib.align(
#                                    x['r2_transposon'],
#                                    r2_transposon)['editDistance'],
#                                axis=1))

#             summarised_dict = primer_df\
#                 .sort_values('edit_dist')\
#                 .groupby('edit_dist')[['edit_dist']]\
#                 .count()\
#                 .rename(columns={'edit_dist': 'tally'})\
#                 .loc[:10, :]\
#                 .to_dict(orient='index')

#         output_dict = {'edit_dist': [], 'tally': []}
#         for edit_dist in range(5):
#             output_dict['edit_dist'].append(edit_dist)
#             output_dict['tally'].append(
#                 summarised_dict.get(edit_dist, {}).get('tally', 0))

#         return pd.DataFrame(output_dict)

#     def _summarize_r2_transposon(self, id_to_bc_df: pd.DataFrame, r1_primer: str, r2_transposon: str) -> pd.DataFrame:
#         """group by a given r2_transposon seq and return the number of r1_primer seqs within 0 to 4 edit distance of the 
#         expected r1_primer seq

#         Args:
#             id_to_bc_df (pd.DataFrame): _description_
#             r1_primer (str): _description_
#             r2_transposon (str): _description_

#         Returns:
#             pd.DataFrame: dataframe with 5 rows
#         """

#         trans_df = id_to_bc_df[id_to_bc_df.r2_transposon == r2_transposon]
#         summarised_dict = {}
#         if not trans_df.shape[0] == 0:
#             trans_df = trans_df\
#                 .assign(edit_dist=trans_df
#                         .apply(lambda x:
#                                edlib.align(
#                                    x['r1_primer'],
#                                    r1_primer)['editDistance'],
#                                axis=1))

#             summarised_dict = trans_df\
#                 .sort_values('edit_dist')\
#                 .groupby('edit_dist')[['edit_dist']]\
#                 .count()\
#                 .rename(columns={'edit_dist': 'tally'})\
#                 .loc[:10, :]\
#                 .to_dict(orient='index')

#         output_dict = {'edit_dist': [], 'tally': []}
#         for edit_dist in range(5):
#             output_dict['edit_dist'].append(edit_dist)
#             output_dict['tally'].append(
#                 summarised_dict.get(edit_dist, {}).get('tally', 0))

#         return pd.DataFrame(output_dict)

#     def _summarize_tf_to_r1_transposon(self, id_to_bc_df: pd.DataFrame, r1_primer: str, r2_transposon: str, r1_transposon: str) -> pd.DataFrame:
#         """_summary_

#         Args:
#             id_to_bc_df (pd.DataFrame): _description_
#             r1_primer (str): _description_
#             r2_transposon (str): _description_
#             r1_transposon (str): _description_

#         Returns:
#             pd.DataFrame: _description_
#         """

#         r1_trans_df = id_to_bc_df[id_to_bc_df.r1_primer == r1_primer]
#         r1_trans_df = r1_trans_df[r1_trans_df.r2_transposon == r2_transposon]

#         # this is here to handle case in which there are no r1, r2 expected
#         # matches
#         summarised_dict = {}
#         if not r1_trans_df.shape[0] == 0:
#             r1_trans_df = r1_trans_df\
#                 .assign(edit_dist=r1_trans_df
#                         .apply(lambda x:
#                                edlib.align(
#                                    x['r1_transposon'],
#                                    r1_transposon)['editDistance'],
#                                axis=1))

#             summarised_dict = r1_trans_df\
#                 .sort_values('edit_dist')\
#                 .groupby('edit_dist')[['edit_dist']]\
#                 .count()\
#                 .rename(columns={'edit_dist': 'tally'})\
#                 .loc[:10, :]\
#                 .to_dict(orient='index')

#         output_dict = {'edit_dist': [], 'tally': []}
#         for edit_dist in range(5):
#             output_dict['edit_dist'].append(edit_dist)
#             output_dict['tally'].append(
#                 summarised_dict.get(edit_dist, {}).get('tally', 0))

#         return pd.DataFrame(output_dict)

#     def add_read_qc(self, barcode_details: BarcodeParser, id_to_bc_map_path: str, split_char: str = 'x') -> None:
#         """_summary_

#         Args:
#             batch (str): _description_
#             barcode_details_json (BarcodeParser): _description_
#             id_to_bc_map_path (str): _description_
#             split_char (str): the character on which to split the TF from the replicate in the barcode_details

#         Raises:
#             FileNotFoundError: _description_
#         """
#         if not os.path.exists(id_to_bc_map_path):
#             raise FileNotFoundError(f'{id_to_bc_map_path} DNE')
#         if not id_to_bc_map_path.endswith('tsv'):
#             logging.warning('%s does not end with tsv' %
#                             id_to_bc_map_path)  # pylint:disable=W1201,C0209

#         batch = barcode_details.barcode_dict['batch']
#         batch_sql = f"SELECT * FROM batch WHERE batch = '{batch}'"

#         cur = self.con.cursor()
#         batch_records = self.db_execute(cur, batch_sql).fetchall()

#         id_to_bc_df = pd.read_csv(id_to_bc_map_path, sep='\t')

#         tf_to_bc_dict = {v.split(split_char)[
#             0]: k for k, v in barcode_details.barcode_dict['components']['tf']['map'].items()}

#         # TODO  address hardcoding
#         # possibly extend BarcodeParser with organism specific settings
#         r1_primer_indices = [0, 5]
#         r2_transposon_indicies = [5, 13]

#         r1_transposon = barcode_details.barcode_dict['components']['r1_transposon']['map'][0]

#         for record in batch_records:
#             # TODO address hardcoding in extracting from records
#             if not record['tf'] == 'undetermined':
#                 tf_seq = tf_to_bc_dict[record['tf']]
#                 r1_primer = tf_seq[r1_primer_indices[0]:r1_primer_indices[1]]
#                 r2_transposon = tf_seq[r2_transposon_indicies[0]                                       :r2_transposon_indicies[1]]

#                 # TODO there is a lot of repeated code in each of the summarize
#                 # fucntions -- unify
#                 qc_summaries = {
#                     'qc_r1_to_r2_tf': self._summarize_r1_primer(
#                         id_to_bc_df,
#                         r1_primer,
#                         r2_transposon),
#                     'qc_r2_to_r1_tf': self._summarize_r2_transposon(
#                         id_to_bc_df,
#                         r1_primer,
#                         r2_transposon),
#                     'qc_tf_to_transposon': self._summarize_tf_to_r1_transposon(
#                         id_to_bc_df,
#                         r1_primer,
#                         r2_transposon,
#                         r1_transposon)
#                 }

#                 for k, v in qc_summaries.items():
#                     self._update_qc_table(
#                         k,
#                         'batch_id',
#                         record['id'],
#                         v.to_dict(orient='records'))

#     def consolidate_tables(self, table_class: Literal['background', 'experiment']) -> bool:
#         """_summary_

#         Args:
#             table_class (Literal[&#39;background&#39;,&#39;experiment&#39;]): _description_

#         Raises:
#             NotImplementedError: _description_

#         Returns:
#             bool: _description_
#         """
#         # take all tables which start with the prefix 'background_' and collapse them
#         # into a single table indexed by the batch_id column
#         raise NotImplementedError

#     def _get_passing_batch_ids(self, batch_id: str, **kwargs) -> list:
#         """_summary_

#         Args:
#             batch_id (str): _description_

#         Raises:
#             NotImplementedError: _description_

#         Returns:
#             list: _description_
#         """
#         raise NotImplementedError()

#     def peak_caller(self, replicate_handling: Literal['separate', 'sum'] = 'separate', poisson_pseudocount: float = 0.2, if_exists: str = 'fail', *args, **kwargs) -> None:
#         """Call Peaks and add the result to the database.

#         Args:
#             replicate_handling (str, ['separate','sum'], optional): How to handle replicates. 'sum' 
#             will add hops by region for all replicates. 'separate' will call peaks for each 
#             replicate separately. Defaults to separate.
#             *args (list): additional positional arguments -- currently unused
#             poisson_pseudocount (float, optional): peudocount to add to poisson pvalue calculation. Defaults to 0.2.
#             **kwargs (dict): Use the following combination of keyword arguments 
#             to direct the peak_caller method:
#                 regions='regions_tbl',background='background_tbl',experiment='experiment_tbl' -- Call peaks with background. 

#         Raises:
#             AttributeError: _description_
#         """
#         # note that the idea of using the kwargs argument is to simulate
#         # overloading, so that a different set of kwargs could be used to call
#         # a different peak caller. If another peak caller is added at some point,
#         # this will need to turn into an if, elif, ... statement with this
#         # error handling ahppening in the final 'else'
#         if {'regions', 'background', 'experiment'} != set(kwargs):
#             raise KeyError(f'The combination of tables {kwargs} '
#                            'does not match any expected combinations -- cannot find an '
#                            'appropriate peak calling method')

#         # TODO handle undetermined case better!
#         if kwargs.get('experiment') != 'undetermined':

#             # create the tablename
#             join_sql = self._regions_background_expr_sql(
#                 kwargs.get('regions'),
#                 kwargs.get('background'),
#                 kwargs.get('experiment'))

#             quant_df = pd.read_sql_query(join_sql, self.con)

#             total_hops_dict = \
#                 {'background': self.get_total_hops(kwargs.get('background'))}

#             # group by batch_id (replicates)
#             batch_grouped_df = quant_df.groupby('batch_id')
#             # add total hops for each replicate to the total_hop_dict
#             for group, df in batch_grouped_df:
#                 total_hops_dict[group] = \
#                     self.get_total_hops(kwargs.get('experiment'), group)

#             if replicate_handling == 'separate':
#                 # call peaks
#                 output_df = \
#                     call_peaks_with_background(
#                         batch_grouped_df,
#                         total_hops_dict, poisson_pseudocount)
#                 # send the table to the database

#             elif replicate_handling == 'sum':
#                 passing_batch_ids = \
#                     self._get_passing_batch_ids(
#                         quant_df.loc[0, 'batch_id'],
#                         kwargs.get('include_unreviewed', False))
#                 summed_passing_replicate_df = \
#                     quant_df[quant_df.batch_id in passing_batch_ids]
#                 summed_passing_replicate_df['group'] = \
#                     ['all'] * len(summed_passing_replicate_df)
#                 summed_passing_replicate_df = \
#                     summed_passing_replicate_df\
#                     .groupby('group')
#                 # SUM OVER POSITIONS!

#                 # extract the max hops for a given passing replicate
#                 remove_keys = [
#                     k for k in total_hops_dict if k in passing_batch_ids]
#                 max_expr_hops = max([total_hops_dict.pop(k)
#                                     for k in remove_keys])
#                 total_hops_dict['all'] = max_expr_hops

#                 # call peaks
#                 output_df = call_peaks_with_background(
#                     summed_passing_replicate_df,
#                     total_hops_dict,
#                     poisson_pseudocount)
#             else:
#                 raise IOError(f'replicate handling method '
#                               f'{replicate_handling} not recognized.')

#             sig_tablename = kwargs.get('regions') + '_' + kwargs.get('background') + \
#                 "_" + kwargs.get('experiment') + '_' + \
#                 replicate_handling + "_sig"

#             return output_df
