# Convert a fasta formatted file of sequences to a csv file.

using CSV
using DataFrames
using BioSequences


df = DataFrame(id=[], seq=[])
reader = FASTA.Reader(open("/home/dillon/data/uniref50/uniref50.fasta", "r"))
count = 0
for (i, record) in enumerate(reader)
    if i < 30*10^6
        continue
    elseif i == 35*10^6
        break
    end
    try
        if FASTA.hasidentifier(record) && FASTA.hassequence(record)
            push!(df, [FASTA.identifier(record), FASTA.sequence(record)])
        end
        if i % 10^5 == 0
            @printf("%d read\n", i)
        end
    catch exc
        print("Exception encountered: ", exc)
        print(record)
        gc()
    end
end
close(reader)
print("Writing csv.")
CSV.write("uniref50_7.csv", df)
