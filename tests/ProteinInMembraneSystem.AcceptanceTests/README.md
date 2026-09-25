# Browser acceptance environment

Build the product with `scripts/build-local.sh`. Keep browser automation separate from the scientific worker's `out/python` environment:

```bash
python3.11 -m venv out/browser-test-python
out/browser-test-python/bin/python -m pip install -r tests/ProteinInMembraneSystem.AcceptanceTests/requirements.txt
out/browser-test-python/bin/python -m playwright install chromium
PIM_BROWSER_REMOTE=1 out/browser-test-python/bin/python -m unittest discover \
  -s tests/ProteinInMembraneSystem.AcceptanceTests/SelectAndPrepareProtein -p 'test_*.py' -v
```

If Chromium is already installed elsewhere, set `PIM_BROWSER_CHROMIUM` to its executable path. Omit `PIM_BROWSER_REMOTE=1` to run the upload and consequential-decision journeys without a live source provider. Screenshots for Product Vision comparison are written to the ignored `out/browser-acceptance/` directory.
