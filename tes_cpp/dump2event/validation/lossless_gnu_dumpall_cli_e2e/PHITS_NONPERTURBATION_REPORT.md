# PHITS non-perturbation report

**Not tested with a repaired writer in this task.** The Oct 21 2026 historical baseline and diagnostic PHITS reruns used the same seed and had byte-identical (SHA-256 `c11c51299c6c77a7e8a531a640ed72b31c8455fb2d0a17f2b03d3716a2815060`) but malformed 2,736,869,598-byte ASCII Dumpall files; their `[t-deposit]` output hashes matched (`0d87d8b7041c54cfd10e610ac0653ff1eda1844a31ca0691126ce82ff0c5f110`). This supports non-perturbation only for that prior pair, whose writer was still lossy. It does not establish non-perturbation after a writer patch.

No new PHITS run was started because neither a verified lossless writer build nor a reproducible internal diagnostic source/binary pair was available.
