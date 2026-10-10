# Third-party code

The files `core.py`, `ckm.py`, and `results.py` in `icpdag/_vendor/cscn` are unchanged source from [CSCN](https://github.com/mghanw/CSCN), commit `d579c03d77a6435c570957188de356b72e369f00`.

Copyright (c) 2026 Menghan Wang. The upstream MIT license is preserved in [LICENSE](icpdag/_vendor/cscn/LICENSE). Original Git blob identifiers are recorded in [upstream.json](icpdag/_vendor/cscn/upstream.json).

The local adapter uses the official NMF, skeleton, collider, Meek, DAG-completion, and Katz routines. Intervention response selection, compatibility checks, and the factor-CKM-to-gene projection are implemented outside the vendored files. The projection follows CSCN Supplement S1.1.2 and Fig. S2.
