#!/usr/bin/python
from sys import argv

script, input_file = argv
fastq = open(input_file, 'r')

output_file = 'bar.' + input_file
outq = open(output_file, 'w')

count = 0
for line in fastq:
    if count == 0:
        headend = line.find(':')
        header = line[0:headend]
        outq.write(line)
        count += 1
    else:
        if len(line) > 3:
            if line[0:len(header)] == header:
                outq.write(line)
            else:
                outq.write(line[6:])
        else:
            outq.write(line)

fastq.close()
outq.close()
