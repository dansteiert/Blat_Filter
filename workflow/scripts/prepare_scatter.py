import sys
import os

class MafSlicer:

    def __init__(self, in_maf, n, out_dir):

        self.in_maf = in_maf
        self.n = int(n)
        self.out_dir = out_dir

    def slice_maf(self):

        with open(self.in_maf, 'r') as maf:
            header_line_count = 0
            header = maf.readline()
            while header[0] == "#" or not header.strip():
                header_line_count +=1
                header=maf.readline()

            header_line_count += 1
            header_columns = header.strip('\n').split('\t')
            assert header_columns[0] == "Hugo_Symbol"

            variant_count = 0
            shard_count = 0
            maf_shard_file = None
            for line in maf:
                if variant_count % self.n == 0:
                    if maf_shard_file is not None:
                        maf_shard_file.close()
                    shard_count += 1
                    maf_shard_filename = "{0}/shard-{1}.maf".format(self.out_dir, shard_count-1)
                    os.makedirs(os.path.dirname(maf_shard_filename), exist_ok=True)
                    maf_shard_file = open(maf_shard_filename, 'w')
                    maf_shard_file.write(header)
                maf_shard_file.write(line)
                variant_count += 1

            # close the final shard file
            if maf_shard_file is not None:
                maf_shard_file.close()
            
def main():
    maf_slicer = MafSlicer(*sys.argv[1:])
    maf_slicer.slice_maf()

if __name__ == "__main__":
    main()
    
