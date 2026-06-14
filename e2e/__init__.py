"""Browser end-to-end smoke tests (Playwright).

Kept OUT of the stdlib `tests/` suite so the core CI stays dependency-free.
Run separately (after `pip install -r requirements-dev.txt` and
`python -m playwright install chromium`):

    python -m unittest discover -s e2e -v
"""
