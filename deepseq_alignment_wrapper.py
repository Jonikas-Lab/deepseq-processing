#!/usr/bin/env python
""" This is mostly a wrapper around the bowtie alignment program, that also reads metadata from a preprocessing info file and adds more metadata to a new *_info.txt file.

Use the -B option to provide a full bowtie option list - it will be passed to bowtie with no parsing or changes except for specifying input/output files and possibly verbosity levels (the full command-line, as ran, is printed to stdout).

Output - there are actually two outfiles generated (with names based on the -o option - calling the value X here):
    - alignment file - <X>.sam or .map, depending on options passed to bowtie
    - metadata file - <X>_info.txt - contains general information from preprocessing AND alignment: 
        exact command used, date/time, system, user, etc
        bowtie version and summary bowtie output (normally printed to stdout), 
        all the preprocessing metadata, if a preprocessing metadata file was found.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, Oct 2011

USAGE: deepseq_alignment_wrapper.py [options] infile -o outfile_basename """

# basic libraries
import sys, os
# my modules
from general_utilities import write_header_data, print_text_from_file, run_command_and_print_info
from seq_count_and_lengths import seq_count_and_lengths


if __name__ == "__main__":
    """ Allows both running and importing of this file. """

    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-B','--full_bowtie_options', default='-f -m1 -v1 --strata --best --tryhard -S --sam-nosq',  
                      metavar='"TEXT"', help="Full options to pass to bowtie, as a quoted string. Set to NONE to not run "
                      +'bowtie. For help on available options, run "bowtie -h" on the command-line. Default "%default".')
    parser.add_option('-I','--bowtie_index', metavar='BASENAME', default='Chlre4nm_chl-mit_cassette-pMJ013b',  
                      help='Bowtie index. (for help, run "bowtie -h" on the command-line. Default %default.)')
    # MAYBE-TODO add options for commonly-used bowtie options; if I do that, make sure the specific options are merged sensibly with the general option (-B)!  Useful options: input type (fasta/fastq), output type (SAM/native), N mismatches allowed... --time option for bowtie to report how long stuff took...

    ### infile/outfile options
    parser.add_option('-m','--input_metadata_file', metavar='FILE', default='AUTO', 
                      help="Metadata file for the infile. Can be a filename, AUTO to infer from infile name "
                          +"(<infile_basename>_info.txt), or NONE when there is no metadata file. "
                          +"Unless the value is NONE, a warning will be raised if not found. Default: %default")
    parser.add_option('-o','--outfile_basename', metavar='X', default='test', 
                      help="Basename for the two outfiles (which will be X.sam/X.map and X_info.txt). Default %default.")
    # MAYBE-TODO add an option to remove unaligned reads? (the unaligned read count will be kept in the metadata file)
    # MAYBE-TODO stdout verbosity levels?  Not really necessary, since I think I always want bowtie output printed to stdout - unless I end up doing read-count checks before and after too, then I may want those printed or not.

    try:
        (options, [infile]) = parser.parse_args()
    except ValueError:
        parser.print_help()
        sys.exit("Error: exactly one infile required!")
        # bowtie could take multiple infiles, but then I'd have to deal with multiple preprocessing metafiles, no thanks!

    ### outfile and tmpfile names
    outfile_suffix = '.sam' if any([x in options.full_bowtie_options for x in ['-S','--sam']]) else '.map'
    outfile = options.outfile_basename + outfile_suffix
    infofile = options.outfile_basename + '_info.txt'
    tmp_bowtie_info = options.outfile_basename + '_tmp_bowtie_info.txt'

    
    with open(infofile,'w') as INFOFILE:

        ### write header data
        write_header_data(INFOFILE,options)
        INFOFILE.write('\n')

        # MAYBE-TODO do I want to check infile readcounts here?  Especially if they're collapsed-to-unique (since in that case bowtie won't report the correct original counts, just the unique counts).  Make an option for that?  See how it's done in deepseq_preprocessing_wrapper.py using check_readcount

        ### run bowtie
        # run 'bowtie --version' to get that data (print to INFOFILE but not stdout)
        command = "bowtie --version > %s"%tmp_bowtie_info
        run_command_and_print_info(command, INFOFILE, printing=False, shell=True)
        print_text_from_file(tmp_bowtie_info, INFOFILE, printing=False, add_newlines=1)
        # run the actual bowtie alignment command; always print output to stdout as well as INFOFILE
        #   (bowtie actually prints the summary to stderr, not stdout, so I need to print it to stdout in case there's 
        #    an error, so I can see the error message!  Or I could try to detect whether there was an error or not
        #    based on the output contents, but that seems like unnecessary work.)
        command = "bowtie %s %s %s %s 2> %s"%(options.full_bowtie_options, options.bowtie_index,
                                              infile, outfile, tmp_bowtie_info)
        run_command_and_print_info(command, INFOFILE, printing=True, shell=True)
        print_text_from_file(tmp_bowtie_info, INFOFILE, printing=True, add_newlines=1)
        # remove tmpfile
        if os.path.exists(tmp_bowtie_info):     os.remove(tmp_bowtie_info)

        # MAYBE-TODO parse the bowtie summary to make my own read-count etc summary?  Not sure if there's a point.

        # MAYBE-TODO do I want to check outfile readcounts here, especially if they're collapsed-to-unique?  How?  I might want to rewrite basic_programs/seq_count_and_lengths.py to take SAM input as well as fasta/fastq...?

        ### copy preprocessing metadata file to the bottom of the new metadata file
        INFOFILE.write("\n################## Metadata from input preprocessing ##################\n\n")
        if options.input_metadata_file == 'NONE':
            INFOFILE.write('Not looking for a metadata input file, as specified by options\n')
        else:
            if options.input_metadata_file == 'AUTO':
                options.input_metadata_file = os.path.splitext(infile)[0] + '_info.txt'
                text = 'Automatically determining metadata input file name: %s\n'%options.input_metadata_file
                print text,
            else:
                text = 'Metadata input file name provided in options: %s\n'%options.input_metadata_file
            INFOFILE.write(text+'\n')
            if os.path.exists(options.input_metadata_file):
                print_text_from_file(options.input_metadata_file, INFOFILE, printing=False)
            else:
                text = 'Metadata input file %s not found!\n'%options.input_metadata_file
                print text,
                INFOFILE.write(text)

        # MAYBE-TODO could actually parse input metadata file and make sure the number of reads found by bowtie matches the one reported at the end of the file, etc...


# MAYBE-TODO write unit-tests?

