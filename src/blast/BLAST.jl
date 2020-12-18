# Module BLAST
# ============
# 
# A module for running command line BLAST and parsing BLAST output files. 

module BLAST

export
    blastn,
    blastp,
    readblastXML,
    BLASTResult

using BioAlignments
using EzXML
using BioSequences 

include("blastcommandline.jl")

end # module BLAST
