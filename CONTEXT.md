# ViHostPlot Domain

ViHostPlot visualizes host-virus integration events across host genomes and viral genomes. The domain language favors generic integration concepts over virus-specific examples or HPV/HBV-specific column names.

## Language

**Integration Event**:
A single observed connection between one host genomic position and one viral genomic position, optionally carrying sample, support, strand, and method metadata.
_Avoid_: insertion, insert record, HPV breakpoint

**Host Genome**:
The set of host genomic sequences used as the coordinate system for host positions.
_Avoid_: chromosome file when referring to the domain object

**Virus Genome**:
The viral sequence used as the coordinate system for viral positions. It defines the viral sequence name and length, but does not by itself imply a specific feature annotation version.
_Avoid_: virus name only, HPV genome

**Virus Feature**:
A named interval on a virus genome, such as a gene, regulatory region, or functional segment. Virus features are supplied explicitly, with helper functions allowed for common virus/version annotations.
_Avoid_: hard-coded gene track, implicit genome annotation

**Integration Track**:
A circos layer that summarizes or displays integration events on the host-virus coordinate system.
_Avoid_: layout task, plot layer

**Virus Density**:
A summary of integration event counts across bins on the virus genome. By default it counts events, not read/support weights.
_Avoid_: read depth unless explicitly weighted

**Linear Integration Plot**:
A non-circos view that places virus and host coordinates on linear axes and connects paired integration positions.
_Avoid_: strudel as the canonical domain term
