import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq



def generate_paired_reads(record, min_HiFi_Length, window_size, paired_distance, step_size, phred_score, beginning):
    sequence = str(record.seq)
    len_sequence = len(sequence)
    read1_list = []
    read2_list = []

    if len_sequence >= min_HiFi_Length: # only consider reads long enough

        for i in range(beginning, len_sequence - window_size - paired_distance, step_size):
            start_read1 = i
            end_read1 = i + window_size
            start_read2 = i + window_size + paired_distance
            end_read2 = start_read2 + window_size

            if end_read2 < len_sequence:
                read1_list.append(SeqRecord(Seq(sequence[start_read1:end_read1]), id=record.id+'_'+str(i), description="1:N", letter_annotations={"phred_quality": [phred_score]*window_size}))
                read2_list.append(SeqRecord(Seq(sequence[start_read2:end_read2]).reverse_complement(), id=record.id+'_'+str(i), description="2:N", letter_annotations={"phred_quality": [phred_score]*window_size}))

    return read1_list, read2_list


def write_Separated(input_file, output_file_prefix, min_HiFi_Length, window_size, paired_distance, step_size, phred_score):
    R1 = []
    R2 = []
    for record in SeqIO.parse(input_file, "fasta"):
        read1_list, read2_list = generate_paired_reads(record, min_HiFi_Length, window_size, paired_distance, step_size, phred_score)
        R1.append(read1_list)
        R2.append(read2_list)

    SeqIO.write(R1, output_file_prefix+".R1.fastq", "fastq")
    SeqIO.write(R2, output_file_prefix+".R2.fastq", "fastq")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract paired-end reads from PacBIO HiFi reads (or any long reads). The script extract substrings of size window_size \
        with insert size of size paired_distance. This is performed in a sliding window manner, with windows of size window_size.")
    parser.add_argument("-w", "--window_size", type=int, default=200, help="Window Size. Default: 200")
    parser.add_argument("-s", "--step_size", type=int, default=150, help="Step size. Default: 150")
    parser.add_argument("-d", "--paired_distance", type=int, default=1600, help="Distance between paired ends. Default: 1600")
    parser.add_argument("-p", "--phred_score", type=int, default=41, help="phredScore to associate to Reads. Default: 41 (maximum score for Illumina 1.8+ fastq format)")
    parser.add_argument("-l", "--minlength", type=int, default=5000, help="minimum HiFi read length. Default: 50000")
    parser.add_argument("-t", "--type_file", type=str, default="Interlaced", help="Type of fastq file. [Interlaced, Separated] Default: Interlaced")
    parser.add_argument("-b", "--beginning", type=int, default=0, help="Starting position to extract the reads. Default: 0 (from beginning of read)")
    parser.add_argument("input_file", type=str, help="Input fasta file containing HiFi reads")
    parser.add_argument("out_prefix", type=str, help="Output prefix for the fastq file. For Separated files, R1 and R2 are gonne be added")

    #Verify the argument
    args = parser.parse_args()

    if args.type_file not in ['Interlaced', 'Separated']:
        raise ValueError("Argument must be either 'Interlaced' or 'Separated'.")
    
    if args.type_file=="Interlaced":
        print("OK")
        read_records = []
        for record in SeqIO.parse(args.input_file, "fastq"):
            read1_list, read2_list = generate_paired_reads(record, args.minlength, args.window_size, args.paired_distance, args.step_size, args.phred_score, args.beginning)
            # Merge read1_list and read2_list into a single list, while preserving order
            for read1, read2 in zip(read1_list, read2_list):
                read_records.append(read1)
                read_records.append(read2)

        # writing the paired end reads to files
        SeqIO.write(read_records, args.out_prefix+".fastq", "fastq")

    elif args.type_file=="Separated":
        R1 = []
        R2 = []
        for record in SeqIO.parse(args.input_file, "fastq"):
            read1_list, read2_list = generate_paired_reads(record, args.minlength, args.window_size, args.paired_distance, args.step_size, args.phred_score, args.beginning)
            R1 += read1_list
            R2 += read2_list

        SeqIO.write(R1, args.out_prefix+".R1.fastq", "fastq")
        SeqIO.write(R2, args.out_prefix+".R2.fastq", "fastq")





    


