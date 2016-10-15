import urllib

output = open('COGs.txt', 'w')

counter = 0
for line in urllib.urlopen('ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/COG0303/cogs.csv').read().split('\n'):
    if not line or line.startswith('#'):
        continue

    fields = line.split(',')

    COG = fields[0]
    CATs = ', '.join(list(fields[1]))
    annotation = ','.join(fields[2:]).strip()

    output.write('\t'.join([COG, CATs, annotation]) + '\n')
    counter += 1

print 'Information for %d COG entries stored.' % counter

output.close()
