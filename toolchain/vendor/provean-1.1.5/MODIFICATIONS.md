# Local compatibility changes (2026-10-20)

Source: https://sourceforge.net/projects/provean/files/provean-1.1.5.tar.gz

Original archive SHA-256:
`f543a4446f239bbe57ab16bb0ea2be5f6aeebddb58f8025dc38781ad4c8e007c`

Only `src/SequenceDB.cpp` is modified. This is a local patched build, not an
official PROVEAN release. The original GPL v3 license and copyright apply.

- Request `sseqid` from PSI-BLAST and `%i` from blastdbcmd so search results,
  retrieved FASTA, CD-HIT records and saved supporting sets share complete typed
  sequence IDs. This preserves PIR/PRF namespaces (for example `pir|A32634|`
  and `prf||1106162A`), which cannot reliably be retrieved by bare accession.
  Normalize PSI-BLAST local IDs by restoring `lcl|`, matching blastdbcmd `%i`.
  Full-database testing exposed this issue after ordinary accession tests passed.
- Accept bare accessions in batch retrieval; remove the requirement for `|`.
- Use CD-HIT `-d 0` to retain complete IDs, compare them exactly, and check array
  bounds when consuming cluster records. Invalid mappings fail explicitly.
- Preserve BLOSUM62, gap penalties, search E-value 0.1, one search iteration,
  maximum 100,000 hits, 75% CD-HIT identity, length filters, 30 clusters and
  delta alignment score calculations.

Supporting sets saved by this build use full BLAST sequence IDs, including their
namespace. Recompute sets saved by the earlier bare-accession patch or old
GI-based builds before reusing them; saved FASTA remains usable with the explicit
`--subject_sequences` interface, which reclusters those sequences. A different BLAST or NR
release may change retrieved homologs and resulting scores.

From the MutaBind2-GPU repository root, build with
`sh toolchain/scripts/build_provean.sh`. The public workflow is
`./mutabind2 scan provean`; its adapter stages user filenames and only passes
safe internal paths to the native program's legacy shell commands.
