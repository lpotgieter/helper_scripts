#outputs set range from fasta file

from Bio import SeqIO

with open("input_fasta","w") as seq:
        for seq_record in SeqIO.parse("output_fasta", "fasta"):
                seq.write(str(seq_record.format("fasta")) + "\n")
                seq.write(str(seq_record.seq[16071:16339]) + "\n")  #range
