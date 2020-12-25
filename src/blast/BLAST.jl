# Module BLAST
# ============
# 
# A module for running command line BLAST and parsing BLAST output files. 

module BLAST

export
    blastn,
    blastp,
    readblastXML,
    readblastXML2,
    BLASTResult

using BioAlignments
using EzXML
using LightXML
using BioSequences 

include("blastcommandline.jl")

end # module BLAST
