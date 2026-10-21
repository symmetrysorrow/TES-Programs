# Raw Dumpall, PHITS internal scoring, and HDF5 addition mapping

A reusable join script is `scripts/compare_hdf5_to_phits_trace.py`. It joins HDF5 `source_record_id` (raw physical line) to the per-record `DUMP` sequence sidecar using event, particle, type, NCOL, energy, coordinates, and region; then it joins the resulting serial and occurrence to `TALLYSEQ` and PRE/POST rows. Duplicate serials are occurrence-matched because a serialized record can produce multiple scoring additions. It does not infer record correspondence from histograms.

| Run | PHITS internal rows | HDF5 additions | raw transport rows uniquely mapped | duplicate scorer serials | energy max abs error vs internal tlv | event max abs sum error | region / event / particle mismatches | position rows compared / max error |
|---|---:|---:|---:|---:|---:|---:|---|---|
| R100 | 38,077 | 38,077 | 38,245 | 30 | 5.7976626971489864e-15 MeV | 3.26405569239796e-14 MeV | 0 / 0 / 0 | 38,077 / 0 cm |
| R1000 | 195,185 | 195,185 | 197,060 | 126 | 5.830839283627043e-15 MeV | 5.129230373768223e-14 MeV | 0 / 0 / 0 | 195,185 / 0 cm |

For R1000, the PRE tally increment agrees with the corresponding PHITS `tlv` to a maximum of `1.1102230246251565e-16 MeV`. HDF5 energy differs from the actual internal increment by at most `5.818913059729702e-15 MeV`; 0 rows exceed `1e-12 MeV`. Every HDF5 source line maps; no overflow/duplicate/missing occurrence exists. All region IDs are 1000.

Continuous contributions are compared as segment start/end coordinates. Point contributions with a valid scorer point are compared to the corresponding raw endpoint. The R100/R1000 score position or segment endpoints agree with the internal trace coordinates exactly at the printed precision. This is position agreement with the recorded PHITS trace; it does not claim that a segment has a unique point score coordinate.

The row-level comparison CSVs and internal traces are private artifacts under `D:\temp\normal_cli_pb_egs5_n100_20261021` and `D:\temp\phits_n1000_diag_20261021\run`; they are not committed.
