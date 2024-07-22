# Open the input and output files
with open("ancestral_alleles.out", "r") as f_input, open("ancestral_alleles.snpid.txt", "w") as f_output:
    # Iterate through the lines in the input file
    for line in f_input:
        # If the line starts with ">", modify the chromosome (CHR)
        if line.startswith(">"):
            # Remove ">" and split by ":" to get the chromosome and range
            parts = line.strip().replace(">", "").split(":")
            CHR = parts[0] + ":" + parts[1].split("-")[1]
        else:
            # Remove trailing newline characters
            line = line.strip()
            # Write to output file if line is not "N"
            if line != "N":
                f_output.write(f"{CHR}\t{line}\n")
