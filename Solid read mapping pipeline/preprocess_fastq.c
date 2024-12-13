#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

// Maximum line length. Adjust if necessary.
#define MAX_LINE_LENGTH 1048576  // 1 MB

// Function to remove a character at a specific position in a string
void remove_char_at_position(char *str, int pos) {
    int len = strlen(str);
    if (pos < 0 || pos >= len - 1) return;  // Exclude the newline character
    for (int i = pos; i < len; i++) {
        str[i] = str[i + 1];
    }
}

// Function to remove the first 'n' occurrences of a specific character from a string
void remove_first_n_chars(char *str, char c, int n) {
    int count = 0;
    char *read_ptr = str;
    char *write_ptr = str;

    while (*read_ptr) {
        if (*read_ptr == c && count < n) {
            count++;
            read_ptr++;  // Skip this character
            if (count == n) {
                // Copy the rest of the string after removing 'n' characters
                strcpy(write_ptr, read_ptr);
                return;
            }
        } else {
            *write_ptr++ = *read_ptr++;
        }
    }
    // If fewer than 'n' characters were found, ensure the string is null-terminated
    *write_ptr = '\0';
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: preprocess_fastq input.fastq.gz output.fastq.gz\n");
        return EXIT_FAILURE;
    }

    const char *input_filename = argv[1];
    const char *output_filename = argv[2];

    // Open input and output gzip files
    gzFile infile = gzopen(input_filename, "rb");
    if (infile == NULL) {
        fprintf(stderr, "Error: Cannot open input file %s\n", input_filename);
        return EXIT_FAILURE;
    }

    gzFile outfile = gzopen(output_filename, "wb");
    if (outfile == NULL) {
        fprintf(stderr, "Error: Cannot open output file %s\n", output_filename);
        gzclose(infile);
        return EXIT_FAILURE;
    }

    // Allocate buffers for the four lines of a FASTQ record
    char *header = (char *)malloc(MAX_LINE_LENGTH);
    char *seq = (char *)malloc(MAX_LINE_LENGTH);
    char *separator = (char *)malloc(MAX_LINE_LENGTH);
    char *qual = (char *)malloc(MAX_LINE_LENGTH);

    if (header == NULL || seq == NULL || separator == NULL || qual == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        gzclose(infile);
        gzclose(outfile);
        free(header);
        free(seq);
        free(separator);
        free(qual);
        return EXIT_FAILURE;
    }

    // Process the FASTQ file record by record
    while (gzgets(infile, header, MAX_LINE_LENGTH) != NULL) {
        // Read the sequence line
        if (gzgets(infile, seq, MAX_LINE_LENGTH) == NULL) {
            fprintf(stderr, "Error: Unexpected end of file while reading sequence line\n");
            break;
        }

        // Read the separator line
        if (gzgets(infile, separator, MAX_LINE_LENGTH) == NULL) {
            fprintf(stderr, "Error: Unexpected end of file while reading separator line\n");
            break;
        }

        // Read the quality line
        if (gzgets(infile, qual, MAX_LINE_LENGTH) == NULL) {
            fprintf(stderr, "Error: Unexpected end of file while reading quality line\n");
            break;
        }

        // Process the sequence and quality lines based on the conditions
        if (strncmp(seq, "TG", 2) == 0) {
            // Remove 'G' from the sequence line (second character)
            remove_char_at_position(seq, 1);

            // Remove the first two '!' characters from the quality line
            remove_first_n_chars(qual, '!', 2);
        }
        else if (seq[0] == 'T') {
            // Remove the first '!' character from the quality line
            remove_first_n_chars(qual, '!', 1);
        }

        // Write the modified lines to the output file
        if (gzputs(outfile, header) < 0 ||
            gzputs(outfile, seq) < 0 ||
            gzputs(outfile, separator) < 0 ||
            gzputs(outfile, qual) < 0) {
            fprintf(stderr, "Error: Failed to write to output file %s\n", output_filename);
            break;
        }
    }

    // Clean up
    gzclose(infile);
    gzclose(outfile);
    free(header);
    free(seq);
    free(separator);
    free(qual);

    return EXIT_SUCCESS;
}
