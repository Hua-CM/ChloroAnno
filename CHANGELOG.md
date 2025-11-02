## 2025.11
- **Blast Compatibility**: Removed deprecated `Bio.Blast.Applications` implementation to resolve compatibility issues introduced by Biopython 1.86. Migrated to native `subprocess` calls for BLAST command execution.
- **Bio.Seqfeature Compatibility**: Seqfeature object dose not support the `strand` attribute from Biopython 1.82.
- **Bio.SeqRecord.SeqRecord**: The `seq` attribute of `SeqRecord` now must be a `Seq` object instead of a string.