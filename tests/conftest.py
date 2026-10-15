import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

# Fresh clone: put the deposited intermediates in place (copies missing files only) so the
# reproduction tests run without the raw screens or the ATAC-seq pipeline.
from cdtest import intermediates  # noqa: E402

if intermediates.available():
    intermediates.install(verbose=False)
