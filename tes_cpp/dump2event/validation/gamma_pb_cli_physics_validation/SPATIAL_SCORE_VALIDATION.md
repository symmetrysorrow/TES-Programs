# Spatial score and region validation

For both 100-history mono-gamma runs, all HDF5 additions were joined to internal PRE/POST scoring rows through raw source-record identity, Dumpall sequence, and repeated scorer occurrence.

| Beam | Addition positions compared | Unavailable | Region mismatch | Event mismatch | Particle mismatch | Maximum coordinate difference |
|---|---:|---:|---:|---:|---:|---:|
| 0.511 MeV | 250,938 | 0 | 0 | 0 | 0 | 0 cm |
| 2 MeV | 549,639 | 0 | 0 | 0 | 0 | 0 cm |

Point deposits are compared with the scorer representative coordinate recorded by PHITS and its corresponding raw child/endpoint. Continuous transport contributions retain their segment endpoints; a segment is not silently represented as a physical point. All compared scoring regions are Pb region 1000. The result is evidence for the current point/segment semantics in the observed slice; it is not a claim that every segment should be converted to a point deposit.
