import os
import sys

class Gatherer:

    def __init__(self, numRejectedMutationsCsv, numPassedMutationsCsv, rejectedMafsCsv, passedMafsCsv, allMafsCsv, pairName):
        # inputs can either be specified as files or as string literals
        args = { k : v for k, v in locals().items() if k != "self" }
        for k, v in args.items(): 
            if os.path.exists(v):
                with open(v) as f:
                    args[k] = f.read()

        self.rejected_count_list = list(map(int, args["numRejectedMutationsCsv"].split(",")))
        self.passed_count_list = list(map(int, args["numPassedMutationsCsv"].split(",")))

        self.rejected_maf_list = args["rejectedMafsCsv"].split(",")
        self.passed_maf_list = args["passedMafsCsv"].split(",")
        self.all_maf_list = args["allMafsCsv"].split(",")
        self.pairName = args["pairName"]

        # read header line from one of the files
        with open(self.all_maf_list[0], 'r') as maf_file:
            self.header_line = maf_file.readline()

    def gather_counts(self):
        rejected_count_total = sum(self.rejected_count_list)
        passed_count_total = sum(self.passed_count_list)
        with open('{0}.count_rejected_mutations.txt'.format(self.pairName), 'w') as rejected_count_file:
            rejected_count_file.write(str(rejected_count_total)+'\n')
        with open('{0}.count_passed_mutations.txt'.format(self.pairName), 'w') as passed_count_file:
            passed_count_file.write(str(passed_count_total)+'\n')

    def gather_mafs(self):

        with open('{0}.blat.rejected.maf'.format(self.pairName), 'w') as rejected_maf_file:
            rejected_maf_file.write(self.header_line)
            for maf_shard_filename in self.rejected_maf_list:
                with open(maf_shard_filename, 'r') as maf_shard_file:
                    maf_shard_file.readline()
                    for line in maf_shard_file:
                        rejected_maf_file.write(line)

        with open('{0}.blat.passed.maf'.format(self.pairName), 'w') as passed_maf_file:
            passed_maf_file.write(self.header_line)
            for maf_shard_filename in self.passed_maf_list:
                with open(maf_shard_filename, 'r') as maf_shard_file:
                    maf_shard_file.readline()
                    for line in maf_shard_file:
                        passed_maf_file.write(line)


        with open('{0}.blat.all.maf'.format(self.pairName), 'w') as all_maf_file:
            all_maf_file.write(self.header_line)
            for maf_shard_filename in self.all_maf_list:
                with open(maf_shard_filename, 'r') as maf_shard_file:
                    maf_shard_file.readline()
                    for line in maf_shard_file:
                        all_maf_file.write(line)
                        
def main():

    gatherer = Gatherer(*sys.argv[1:])

    gatherer.gather_counts()

    gatherer.gather_mafs()

if __name__ == "__main__":
    main()
