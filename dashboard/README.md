# TerraScreen Dashboard

This dashboard serves as an interactive visualization for the output files produced by the TerraScreen program. The dashboard uses (WebDash)[https://github.com/ibdafna/webdash/]. The majority of the files inside `TerraScreen/dashboard` are copied from the original WebDash GitHub repository for ease of replication. Copying occurred on 2025 09 11. 

# Python files

To develop and test the dashboard locally, a Python file `src/dash_app/app.py` is included. This file can be converted into a TypeScript-friendly wrapper by using `src/dash_app/ts_embed_from_py.py`. When this code runs, it writes out the contents of the new typescript file. It should be run as follows: `python src/dash_app/ts_embed_from_py.py src/dash_app/app.py > src/dash_app/cross_filtering_app.ts` to overwrite the contents of `cross_filtering_app.ts` with the latest version of the Python app.

`pyproject.toml` and `poetry.lock` are provided to make it easier to install the dependencies required to run `src/dash_app/app.py` locally. Simply install Poetry and then run `poetry install .` inside this directory to install the dependencies. `app.py` can then be run as follows: `poetry run python src/dash_app/app.py`.




