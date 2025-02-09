
# Running Python tests

First, a development build of the duckdb Python package that matches the commit
of the duckdb submodule used to build the extension.

```shell
cd duckdb/tools/pythonpkg
pip install -r requirements-dev.txt
pip install .
```

Running the `test_*.py` files requires `pytest`:

```shell
# pip install pytest pyarrow
pytest -vv tests/python
```
