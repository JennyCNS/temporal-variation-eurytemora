
import gzip

vcf_file_path = '/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/PoolSNP/final-maf0.01-mincov50-mincount10-allsites0-excluding20072011.vcf.gz'
output_file_path = '/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/finalfile.vcf'

# Open the VCF file in read mode
with gzip.open(vcf_file_path, 'rt', encoding='utf-8') as vcf_file:
    # Iterate through each line in the file
    with open(output_file_path, 'a', encoding='utf-8') as output_file:
        line_number = 0

        for line in vcf_file:
            if line.startswith('#'):
                continue

            # Increment the line number (Move this line outside the loop)
            # line_number += 1
            # Split the line into fields based on tabs
            fields = line.strip().split('\t')
            # 9th column forward
            fsplit = fields[8].split(':')
            modified_fields = ':'.join([fsplit[0], fsplit[2], fsplit[3], fsplit[4]])
            fields[8] = modified_fields

            idx = 9
            for fcolumn in fields[9:]:
                tmpsplit = fcolumn.split(':')
                tmpjoin = ':'.join([tmpsplit[0], ','.join([tmpsplit[1], tmpsplit[2]]), ':'.join([tmpsplit[3], tmpsplit[4]])])
                fields[idx] = tmpjoin
                idx += 1

            # Write the modified line to the output file
            output_file.write('\t'.join(fields) + '\n')

            line_number += 1
            #print(line_number)

            # Print line number every 50,000 lines
            if line_number % 50000 == 0:
                print(f"Processed {line_number} lines.")




