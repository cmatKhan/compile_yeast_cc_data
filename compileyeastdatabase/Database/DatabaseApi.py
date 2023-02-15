# pylint: disable=W1201,W1203,C0103
import os
import sys
import re
import logging
from typing import Callable, Literal
import warnings
from pathlib import Path, PosixPath
from io import StringIO

from sqlalchemy import create_engine
import pandas as pd
import numpy as np

from compileyeastdatabase.Resources import Resources

# TODO switch to SQLalchemy and remove this
warnings.simplefilter(action='ignore', category=UserWarning)

__all__ = ['DatabaseApi', 'DatabaseAPINew']

logging.getLogger(__name__).addHandler(logging.NullHandler())

# https://github.com/pandas-dev/pandas/issues/14553#issuecomment-778599042
# class for SqlUpsert. Also consider SQLalchemy or that sql package i sent to
# daniel

# discovered the 'row_factory' attribute in writing this. The connection
# row factory is set to sqlite.Row. This is another interesting one:
# see https://docs.python.org/3/library/sqlite.html#connection-objects
# self.con.row_factory = lambda cursor, row: row[0]


class DatabaseAPINew(Resources):

    def __init__(self, db_loc: PosixPath) -> None:

        super().__init__()

        # get absolute path
        try:
            db_loc = db_loc.resolve()
        except AttributeError:
            raise AttributeError("db_loc must be a pathlib.PosixPath object")
        # check that the parent dir exists
        if not os.path.exists(db_loc.parent):
            raise FileNotFoundError(f"Parent directory of {db_loc} "
                                    "does not exist")
        # create connection engine
        self.engine = create_engine(f'sqlite:///{db_loc}', echo=True)

        # read in data from package resources
        chr_map_df = pd.read_csv(StringIO(self.yeast_chr_map))
        background_df_list = [pd.read_csv(StringIO(x)) for
                              x in [self.yeast_background_adh1,
                              self.yeast_background_sir4]]
        regions_df_list = [pd.read_csv(StringIO(x)) for
                           x in [self.yeast_promoters_yiming,
                                 self.yeast_promoters_not_orf]]

        # add data to database
        self.add_chr_map(chr_map_df)
        self.add_background(background_df_list)
        self.add_regions(regions_df_list)
        # create views
        self.create_region_background_view()

    @property
    def engine(self):
        return self._engine

    @engine.setter
    def engine(self, new_engine):
        self._engine = new_engine

    def add_chr_map(self, chr_map_df: pd.DataFrame) -> None:

        sql = """CREATE TABLE "chr_map" (
            "refseq"	TEXT UNIQUE,
            "igenomes"	TEXT UNIQUE,
            "ensembl"	TEXT UNIQUE,
            "ucsc"	TEXT UNIQUE,
            "mitra"	TEXT UNIQUE,
            "seqlength"	INTEGER NOT NULL,
            "numbered"	TEXT UNIQUE,
            PRIMARY KEY("ucsc"));"""

        with self.engine.connect() as conn:
            conn.execute(sql)

        chr_map_df.to_sql('chr_map', self.engine,
                          if_exists='append', index=False)

    def add_background(self, df_list: list) -> None:
        """Add a regions table to the database.

        Args:
            df (pd.DataFrame): A dataframe with columns 'chr', 'start', 'end'
        """

        create_sql = """CREATE TABLE "background" (
            "id"	INTEGER NOT NULL,
            "chr"	TEXT,
            "start"	INTEGER,
            "end"	INTEGER,
            "depth"	INTEGER,
            "strand"	TEXT,
            "annotation"	TEXT,
            "sample"	TEXT,
            PRIMARY KEY("id"),
            FOREIGN KEY("chr") REFERENCES "chr_map"("ucsc"));"""

        index_sql = """CREATE INDEX "background_index" ON "background" (
            "chr",
            "start",
            "strand",
            "sample");"""

        with self.engine.connect() as conn:
            conn.execute(create_sql)
            conn.execute(index_sql)

        for df in df_list:
            df.to_sql('background', self.engine,
                      if_exists='append', index=False)

    def add_regions(self, df_list: list) -> None:

        create_sql = """CREATE TABLE "regions" (
            "id"	INTEGER NOT NULL,
            "chr"	TEXT,
            "start"	INTEGER,
            "end"	INTEGER,
            "name"	TEXT,
            "score"	REAL,
            "strand"	TEXT,
            "sample"	TEXT,
            "common_name"  TEXT DEFAULT 'none',
            PRIMARY KEY("id"),
            FOREIGN KEY("chr") REFERENCES "chr_map"("ucsc"));"""

        index_sql = """CREATE INDEX "regions_index" ON "regions" (
            "chr",
            "start",
            "end",
            "strand",
            "sample");"""

        with self.engine.connect() as conn:
            conn.execute(create_sql)
            conn.execute(index_sql)

        for df in df_list:
            df.to_sql('regions', self.engine,
                      if_exists='append', index=False)

    def create_region_background_view(self) -> None:
        """Create a view that joins the regions and background tables.
        
        The view is called 'region_background_agg' and has the following
        columns:
        chr, start, end, background_hops, regions_sample, background_sample,
        region_id, associated_feature_systematic_id,
        associated_feature_common_name
        
        The view is created by joining the regions and background tables on
        chr and whether background.start is between region.start and region.end 
        on the same chr. The background_hops column is the number of background
        hops in the region. The associated_feature_systematic_id and
        associated_feature_common_name columns are taken from the regions
        table."""

        drop_sql = """DROP VIEW IF EXISTS "main"."region_background_agg";"""
        create_view_sql = \
            """CREATE VIEW region_background_agg AS 
                SELECT r.chr as chr,
                    r.start as start,
                    r.end as end, 
                    COUNT(*) as background_hops, 
                    r.sample as regions_sample,
                    x.sample as background_sample,
                    r.id as region_id,
                    r.name as associated_feature_systematic_id,
                    r.common_name as associated_feature_common_name
                FROM regions as r 
                LEFT JOIN background as x 
                WHERE (x.chr = r.chr AND x.start BETWEEN r.start AND r.end) 
                GROUP BY  regions_sample, background_sample,r.chr,
                          r.start,r.end, r.id
                ORDER BY regions_sample, background_sample, chr, start"""

        bg_total_hops_drop_sql = """DROP VIEW IF EXISTS "main"."total_bg_hops";"""  # noqa
        bg_total_hops_create_sql = \
            """CREATE VIEW "total_bg_hops" AS SELECT sample, COUNT(*) as hops
           FROM background
           GROUP BY sample"""

        with self.engine.connect() as conn:
            conn.execute(drop_sql)
            conn.execute(create_view_sql)
            conn.execute(bg_total_hops_drop_sql)
            conn.execute(bg_total_hops_create_sql)
