# BLAST+ Wrapper
# ==============
#
# Wrapper for BLAST+ command line functions.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

struct BLASTResult
    bitscore::Float64
    expect::Float64
    queryname::String
    hitname::String
    hittaxid::String 
    hitid::String
    hitaccession::String
    hitlen::Int64
    hit::LongSequence
    identity::Int64
    positive::Int64
    gaps::Int64
    queryfrom::Int64
    queryto::Int64
    hitfrom::Int64
    hitto::Int64 
    alignment::AlignedSequence
end



"""
    readblastXML2(blastrun::AbstractString)
Parse XML2 output of a blast run (blastn format #16). Input is an XML string eg:
```julia
results = read(open("blast_results.xml"), String)
readblastXML(results)
```
Returns Vector{BLASTResult} with the sequence of the hit, the Alignment with query sequence, bitscore and expect value
"""
function readblastXML2(blastrun::AbstractString; seqtype="nucl")
    results = BLASTResult[]
    xdoc = parse_string(blastrun)
    xroot = LightXML.root(xdoc)
    blast_outputs = get_elements_by_tagname(xroot, "BlastOutput2")
    for blast_output in blast_outputs
        report = get_elements_by_tagname(blast_output, "report")[1]
        Report = get_elements_by_tagname(report, "Report")[1]
        Search = Report["results"][1]["Results"][1]["search"][1]["Search"][1]
        queryname = content(Search["query-title"][1]) 
        for hit in Search["hits"][1]["Hit"]
            #println(hit["num"])
            hitdescr = hit["description"][1]["HitDescr"][1]
            hitname = content(hitdescr["title"][1])
            hitid = content(hitdescr["id"][1])
            hitaccession = content(hitdescr["accession"][1])
            hittaxid = content(hitdescr["taxid"][1])
            hitlen = parse(Int64,content(hit["len"][1]))
            for hsp in hit["hsps"][1]["Hsp"]
                gaps = parse(Int64,content(hsp["gaps"][1]))
                qseq = content(hsp["qseq"][1])
                hseq = content(hsp["hseq"][1])
                if seqtype == "nucl"
                    qseq = LongSequence{DNAAlphabet{4}}(qseq)
                    hseq = LongSequence{DNAAlphabet{4}}(hseq)
                    positive = 0
                elseif seqtype == "prot"
                    qseq = LongSequence{AminoAcidAlphabet}(qseq)
                    hseq = LongSequence{AminoAcidAlphabet}(hseq)
                    positive = parse(Int64,content(hsp["positive"][1]))
                else
                    throw(error("Please use \"nucl\" or \"prot\" for seqtype"))
                end
                aln = AlignedSequence(qseq, hseq)
                queryfrom = parse(Int64,content(hsp["query-from"][1]))
                queryto = parse(Int64,content(hsp["query-to"][1]))
                hitfrom = parse(Int64,content(hsp["hit-from"][1]))
                hitto = parse(Int64,content(hsp["hit-to"][1]))
                bitscore = parse(Float64,content(hsp["bit-score"][1]))
                score = parse(Int64,content(hsp["score"][1]))
                expect = parse(Float64,content(hsp["evalue"][1]))
                identity = parse(Int64,content(hsp["identity"][1]))
                gaps = parse(Int64,content(hsp["gaps"][1]))
                push!(results, BLASTResult(bitscore, expect, queryname, hitname,hittaxid, hitid, hitaccession, hitlen, hseq, identity, positive, gaps, queryfrom, queryto, hitfrom, hitto, aln))
            end

        end
    end
    return results
end

"""
    readblastXML(blastrun::AbstractString)
Parse XML output of a blast run. Input is an XML string eg:
```julia
results = read(open("blast_results.xml"), String)
readblastXML(results)
```
Returns Vector{BLASTResult} with the sequence of the hit, the Alignment with query sequence, bitscore and expect value
"""
function readblastXML(blastrun::AbstractString; seqtype="nucl")
    dc = EzXML.parsexml(blastrun)
    rt = EzXML.root(dc)
    results = BLASTResult[]
    for iteration in findall("/BlastOutput/BlastOutput_iterations/Iteration", rt)
        queryname = EzXML.nodecontent(findfirst("Iteration_query-def", iteration))
        for hit in findall("Iteration_hits", iteration)
            if EzXML.countelements(hit) > 0
                hitname = EzXML.nodecontent(findfirst("./Hit/Hit_def", hit))
                hitid = EzXML.nodecontent(findfirst("./Hit/Hit_id", hit))
                hitaccession = EzXML.nodecontent(findfirst("./Hit/Hit_accession", hit))
                hitlen = parse(Int64,EzXML.nodecontent(findfirst("./Hit/Hit_len", hit)))  
                hsps_section = EzXML.nodecontent(findfirst("./Hit/Hit_hsps", hit))
                hsps = findall("./Hit/Hit_hsps/Hsp",hit)
                for hsp in hsps
                    if seqtype == "nucl"
                        qseq = LongSequence{DNAAlphabet{4}}(EzXML.nodecontent(findfirst("./Hsp_qseq", hsp)))
                        hseq = LongSequence{DNAAlphabet{4}}(EzXML.nodecontent(findfirst("./Hsp_hseq", hsp)))
                    elseif seqtype == "prot"
                        qseq = LongSequence{AminoAcidAlphabet}(EzXML.nodecontent(findfirst("./Hsp_qseq", hsp)))
                        hseq = LongSequence{AminoAcidAlphabet}(EzXML.nodecontent(findfirst("./Hsp_hseq", hsp)))
                    else
                        throw(error("Please use \"nucl\" or \"prot\" for seqtype"))
                    end
                    aln = AlignedSequence(qseq, hseq)
                    queryfrom = parse(Int64, EzXML.nodecontent(findfirst("./Hsp_query-from", hsp)))
                    queryto = parse(Int64, EzXML.nodecontent(findfirst("./Hsp_query-to", hsp)))
                    hitfrom = parse(Int64, EzXML.nodecontent(findfirst("./Hsp_hit-from", hsp)))
                    hitto = parse(Int64, EzXML.nodecontent(findfirst("./Hsp_hit-to", hsp)))
                    bitscore = parse(Float64, EzXML.nodecontent(findfirst("./Hsp_bit-score", hsp)))
                    expect = parse(Float64, EzXML.nodecontent(findfirst("./Hsp_evalue", hsp)))
                    identity = parse(Int64, EzXML.nodecontent(findfirst("./Hsp_identity", hsp)))
                    positive = parse(Int64, EzXML.nodecontent(findfirst("./Hsp_positive", hsp)))
                    gaps = parse(Int64, EzXML.nodecontent(findfirst("./Hsp_gaps", hsp)))
                    push!(results, BLASTResult(bitscore, expect, queryname,hitname,"0",hitid, hitaccession, hitlen, hseq, identity, positive, gaps, queryfrom, queryto, hitfrom, hitto, aln))
                end 
            end
        end
    end
    return results
end

"""
`readblastXML(blastrun::Cmd)`
Parse command line blast query with XML output. Input is the blast command line command, eg:
```julia
blastresults = `blastn -query seq1.fasta -db some_database -outfmt 5`
readblastXML(blastresults)
```
Returns Vector{BLASTResult} with the sequence of the hit, the Alignment with query sequence, bitscore and expect value
"""
function readblastXML(blastrun::Cmd; seqtype="nucl")
    return readblastXML(read(blastrun, String), seqtype=seqtype)
end


"""
`blastn(query, subject, flags...)``
Runs blastn on `query` against `subject`.
    Subjects and queries may be file names (as strings), LongSequence{DNAAlphabet{4}} type or
    Array of LongSequence{DNAAlphabet{4}}.
    May include optional `flag`s such as `["-perc_identity", 95,]`. Do not use `-outfmt`.
"""
function blastn(query::AbstractString, subject::AbstractString, flags=[]; db::Bool=false)
    if db
        results = readblastXML(`blastn -query $query -db $subject $flags -outfmt 5`)
    else
        results = readblastXML(`blastn -query $query -subject $subject $flags -outfmt 5`)
    end
    return results
end

function blastn(query::LongSequence{DNAAlphabet{4}}, subject::LongSequence{DNAAlphabet{4}}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastn(querypath, subjectpath, flags)
end

function blastn(query::LongSequence{DNAAlphabet{4}}, subject::Vector{LongSequence{DNAAlphabet{4}}}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    blastn(querypath, subjectpath, flags)
end

function blastn(query::LongSequence{DNAAlphabet{4}}, subject::AbstractString, flags=[]; db::Bool=false)
    querypath = makefasta(query)
    if db
        return blastn(querypath, subject, flags, db=true)
    else
        return blastn(querypath, subject, flags)
    end
end

function blastn(query::Vector{LongSequence{DNAAlphabet{4}}}, subject::Vector{LongSequence{DNAAlphabet{4}}}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastn(querypath, subjectpath, flags)
end

function blastn(query::Vector{LongSequence{DNAAlphabet{4}}}, subject::AbstractString, flags=[]; db::Bool=false)
    querypath = makefasta(query)
    if db
        return blastn(querypath, subject, flags, db=true)
    else
        return blastn(querypath, subject, flags)
    end
end

function blastn(query::AbstractString, subject::Vector{LongSequence{DNAAlphabet{4}}}, flags=[])
    subjectpath = makefasta(subject)
    return blastn(query, subjectpath, flags)
end

"""
`blastp(query, subject, flags...)``
Runs blastn on `query` against `subject`.
    Subjects and queries may be file names (as strings), `BioSequence{LongSequence{AminoAcidAlphabet}}` type or
    Array of `BioSequence{LongSequence{AminoAcidAlphabet}}`.
    May include optional `flag`s such as `["-perc_identity", 95,]`. Do not use `-outfmt`.
"""
function blastp(query::AbstractString, subject::AbstractString, flags=[]; db::Bool=false)
    if db
        results = readblastXML(`blastp -query $query -db $subject $flags -outfmt 5`, seqtype = "prot")
    else
        results = readblastXML(`blastp -query $query -subject $subject $flags -outfmt 5`, seqtype = "prot")
    end
    return results
end

function blastp(query::LongSequence{AminoAcidAlphabet}, subject::LongSequence{AminoAcidAlphabet}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastp(querypath, subjectpath, flags)
end

function blastp(query::LongSequence{AminoAcidAlphabet}, subject::Vector{LongSequence{AminoAcidAlphabet}}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastp(querypath, subjectpath, flags)
end

function blastp(query::LongSequence{AminoAcidAlphabet}, subject::AbstractString, flags=[]; db::Bool=false)
    querypath = makefasta(query)
    if db
        return blastp(querypath, subject, flags, db=true)
    else
        return blastp(querypath, subject, flags)
    end
end

function blastp(query::Vector{LongSequence{AminoAcidAlphabet}}, subject::Vector{LongSequence{AminoAcidAlphabet}}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastp(querypath, subjectpath, flags)
end

function blastp(query::Vector{LongSequence{AminoAcidAlphabet}}, subject::AbstractString, flags=[]; db::Bool=false)
    querypath = makefasta(query)
    if db
        return blastp(querypath, subject, flags, db=true)
    else
        return blastp(querypath, subject, flags)
    end
end

function blastp(query::AbstractString, subject::Vector{LongSequence{AminoAcidAlphabet}}, flags=[])
    subjectpath = makefasta(subject)
    return blastp(query, subjectpath, flags)
end

# Create temporary fasta-formated file for blasting.
function makefasta(sequence::BioSequence)
    path, io = mktemp()
    write(io, ">$path\n$(convert(String, sequence))\n")
    close(io)
    return path
end

# Create temporary multi fasta-formated file for blasting.
function makefasta(sequences::Vector{T}) where T <: BioSequence
    path, io = mktemp()
    counter = 1
    for sequence in sequences
        write(io, ">$path$counter\n$(convert(String, sequence))\n")
        counter += 1
    end
    close(io)
    return path
end
