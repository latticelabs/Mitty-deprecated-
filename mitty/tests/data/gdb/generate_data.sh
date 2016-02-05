# Run this as v 1.33.1
genomes generate gdb.json --ref ../chimera.fa --db gdb-1.33.1.h5 -v -p
genomes genome-file summary gdb-1.33.1.h5 # To verify stuff works
cp gdb-1.33.1.h5 gdb-1.34.0.h5  # In pre for conversion
# Run as after v1.33.1
migratedb gdb-1.34.0.h5
# Run as v1.33.1 - the DB now works on both versions!
genomes genome-file summary gdb-1.34.0.h5