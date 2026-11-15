"""Checkpoint text dtypes survive pandas 2/3 inference and empty event tables."""
import unittest

import numpy as np
import pandas as pd

from rustbears.uncertainty.checkpoint import frame_arrays, restore_frame


class FrameDtypeTests(unittest.TestCase):
    def test_text_dtype_and_missing_values_round_trip(self):
        dtypes = [object, pd.StringDtype(storage="python")]
        # The NumPy-NaN string dtype is new in pandas 3.
        try:
            dtypes.append(pd.StringDtype(storage="python", na_value=np.nan))
        except TypeError:
            pass
        for dtype in dtypes:
            for values in (["d", "e", None], []):
                with self.subTest(dtype=repr(dtype), count=len(values)):
                    frame = pd.DataFrame({"event_type": pd.Series(values, dtype=dtype)})
                    metadata, arrays = frame_arrays(frame, "events")
                    self.assertFalse(any(a.dtype.hasobject for a in arrays.values()))
                    restored = restore_frame(metadata, arrays)
                    pd.testing.assert_frame_equal(frame, restored)


if __name__ == "__main__":
    unittest.main()
