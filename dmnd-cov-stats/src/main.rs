use std::io::Write;
use std::io::BufRead;
use std::io::BufWriter;
use std::fs::File; 
use std::io;

// Struct to store the bins
struct Bin {
    lower: f64,
    upper: f64,
    middle: f64,
    count: u64,
}

// Struct to contain the vector of bins, sequence length and a boolean to indicate if the sequence was found in the diamond output
struct Sequence {
    bins: Vec<Bin>,
    sequence_length: u64,
    found: bool,
}

// Function to prepare the bins
// Returns a dictionary with the sequence names as keys, and the Sequence struct as values
fn prepare_bins(sequence_lengths_file: &str, number_of_bins: u64) -> std::collections::HashMap<String, Sequence> {
    // Read the sequence lengths file
    let mut coverage_stats: std::collections::HashMap<String, Sequence> = std::collections::HashMap::new();
    let file: std::fs::File = std::fs::File::open(sequence_lengths_file).unwrap();
    let reader: std::io::BufReader<std::fs::File> = std::io::BufReader::new(file);
    for line in reader.lines() {
        let line: String = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();
        let sequence_name: String = fields[0].to_string();
        let sequence_length: u64 = fields[1].parse().unwrap();
        let mut bins_vec: Vec<Bin> = Vec::new();
        for i in 0..number_of_bins {
            let lower: f64 = i as f64 * sequence_length as f64 / number_of_bins as f64;
            let upper: f64 = (i + 1) as f64 * sequence_length as f64 / number_of_bins as f64;
            let middle: f64 = lower + (upper - lower) / 2.0;         
            bins_vec.push(Bin { lower, upper, middle, count: 0 });
        }
        coverage_stats.insert(sequence_name, Sequence { bins: bins_vec, sequence_length, found: false });
    }
    coverage_stats
}

fn main() {

    // Command line arguments
    // We have the following arguments:
    // 1. Path to the sequence lengths file
    // 2. Path to the diamond output file, can also be stdin, if - is provided
    // 3. Output prefix
    // 4. Number of bins
    // Parse the command line arguments
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: {} <sequence_lengths_file> <diamond_output_file> <output_prefix> <number_of_bins>", args[0]);
        std::process::exit(1);
    }

    // Read the sequence lengths file and prepare the bins
    let sequence_lengths_file = &args[1];
    let number_of_bins = args[4].parse().unwrap();
    let mut sequence_bins = prepare_bins(sequence_lengths_file, number_of_bins);

    // Process the diamond output file line by line. If - is given as the diamond output file, read from stdin
    // The diamond output file is a tab-separated file with the following columns:
    // 'seqid', 'sseqid', 'pident', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'cigar'
    // We need the sseqid, sstart and send columns to calculate the coverage 
    let diamond_output_file = &args[2];
    let reader: Box<dyn std::io::BufRead> = if diamond_output_file == "-" {
        Box::new(std::io::BufReader::new(std::io::stdin()))
    } else {
        Box::new(std::io::BufReader::new(std::fs::File::open(diamond_output_file).unwrap()))
    };
    for line in reader.lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();
        let sequence_name = fields[1];
        let start = fields[5].parse::<u64>().unwrap();
        let end = fields[6].parse::<u64>().unwrap();
        // Check if end is less than start, if so, swap the values
        let (start, end) = if start < end {
            (start, end)
        } else {
            (end, start)
        };
        let sequence_length = sequence_bins.get_mut(sequence_name).unwrap().sequence_length;
        let sequence_bins: &mut Sequence = sequence_bins.get_mut(sequence_name).unwrap();
        let bins = &mut sequence_bins.bins;
        // We count it as a hit in a bin if the start or end of the alignment falls within the bin
        for bin in bins {
            if (start as f64 >= bin.lower && start as f64 <= bin.upper) || (end as f64 >= bin.lower as f64 && end as f64 <= bin.upper) {
                bin.count += 1;
                sequence_bins.found = true;
            }
        }
    }

    // Write the output to a file. If the output prefix is -, set the output to stdout
    let output_prefix = &args[3];
    
    let mut output: Box<dyn Write> = if output_prefix == "-" {
        Box::new(BufWriter::new(io::stdout()))
    } else {
        Box::new(BufWriter::new(File::create(format!("{}.tsv", output_prefix)).unwrap()))
    };



    // Write the header
    // It contains the following columns:
    // sequence_name, found, sequence_length, bin1, bin2, ..., binN
    // The found column is a boolean indicating if the sequence was found in the diamond output
    //The bin names are bin numbers normalized between 0 and 1, e.g if we have 10 bins, the bin names will be 0.05, 0.015, ..., 0.95
    writeln!(output, "sequence_name\tfound\tsequence_length\t{}", (1..=number_of_bins).map(|x| format!("{}", (x as f64 - 0.5) / number_of_bins as f64)).collect::<Vec<String>>().join("\t")).unwrap();
    
    // Write the data
    for (sequence_name, sequence) in sequence_bins.iter() {
        writeln!(output, "{}\t{}\t{}\t{}", sequence_name, sequence.found as u64, sequence.sequence_length, sequence.bins.iter().map(|bin| bin.count.to_string()).collect::<Vec<String>>().join("\t")).unwrap();
    }

    
    

    // Not used alternative
    // bin1, bin2, ..., binN are the counts of the bins, and named as the middle value of the bin.
    // The middle value is contained in the sequence_bins struct, we just take the first entry
    //writeln!(output, "sequence_name\tfound\tsequence_length\t{}", sequence_bins.values().next().unwrap().bins.iter().map(|bin| bin.middle.to_string()).collect::<Vec<String>>().join("\t")).unwrap();

}
