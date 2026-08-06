"""Put the package dir on sys.path so tests can use flat imports
(`import schema`, `import adapters`, `import aggregate`) without requiring
`__init__.py` files up the `applications/` tree or repo-wide pytest config.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
