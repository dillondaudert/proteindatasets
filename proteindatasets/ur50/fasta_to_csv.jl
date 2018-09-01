# Convert a fasta formatted file of sequences to a csv file.

using CSV
using DataFrames
using BioSequences


df = DataFrame(id=[], seq=[])
reader = FASTA.Reader(open("/home/dillon/data/cUR50/cur50.fasta", "r"))
count = 0
for (i, record) in enumerate(reader)
    try
        if FASTA.hasidentifier(record) && FASTA.hassequence(record)
            push!(df, [FASTA.identifier(record), FASTA.sequence(record)])
        end
        if i % 10^5 == 0
            print("$i read\n")
        end
    catch exc
        print("Exception encountered: ", exc)
        print(record)
    end
end
close(reader)
print("Writing csv.")
CSV.write("/home/dillon/data/cUR50/cur50.csv", df)
