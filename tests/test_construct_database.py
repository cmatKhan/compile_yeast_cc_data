import os
import pathlib

from compileyeastdatabase.Database.DatabaseApi import DatabaseAPINew


def test_database_constructor(tmpdir):

    test_db = pathlib.Path(os.path.join(tmpdir, 'test.db'))

    DatabaseAPINew(test_db)

    assert os.path.exists(test_db)
