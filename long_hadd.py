import logging
import subprocess
import argparse
from os.path import isfile
import sys


if __name__ == "__main__":

    parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter,
        description = "hadd files in batches of N (default 50)",
        epilog = "Example:\n$ python long_hadd.py out.root in*root")

    parser.add_argument("output_filename", help="the name of the output file")
    parser.add_argument("-n", "--number-of-files", type=int, default=20, help="amount of files per 1 run of hadd")
    parser.add_argument("-t", "--test",  help='test run with echo instead of hadd', action="store_true")
    parser.add_argument("-d", "--debug", help='debug level of logging', action="store_true")
    parser.add_argument("-i", "--input-files", nargs='+', help='files to hadd')

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    logging.debug(args) # whatch out -- all input files!

    output_file = args.output_filename
    if isfile(output_file):
         print "the file exists " + output_file
         sys.exit(1)

    N = args.number_of_files
    files = args.input_files

    if args.test:
        com = 'echo'
    else:
        com = 'hadd'

    prev_batch = ''
    for i, batch in enumerate(files[i:i+N] for i in xrange(0, len(files), N)):
        # construct command for processing new batch:
        # it outputs into a separate temp file, which needs to be removed later
        # and adds the file of the previus batch to the hadd
        logging.debug("batch %s" % i)
        batch_out = args.output_filename + str(i)
        com_b = com + ' ' + batch_out + ' ' + ' '.join(batch) + ' ' + prev_batch
        logging.debug(prev_batch)

        out = subprocess.check_output(com_b, shell=True)
        print out

        # clean the previous batch output if it extists
        # and set new brev_batch filename
        if isfile(prev_batch):
            print subprocess.check_output('rm ' + prev_batch, shell=True)
        prev_batch = batch_out

    print subprocess.check_output('mv %s %s' % (prev_batch, output_file), shell=True)


