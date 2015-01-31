#! /usr/bin/env python2.7

"""
______
 -- Weronika Patena, 2015
USAGE: _____
"""

# standard library
from __future__ import division
import sys
import unittest
# other packages
# my modules

class CassetteStrippingError(Exception): pass

def strip_cassette_seq(input_seq, cassette_seq, before_cassette_lens, min_len=5):
    """ Split input_seq into sections before and after cassette_seq, allowing cassette truncations.

    The input_seq is expected to look like this: xxxxxxxxxCCCCCCCCCCCCCCCyyyyyyyy, 
    where CCCCCCCCCC is either all of cassette_seq or a prefix of it (of at least min_len), 
    and the length of xxxxxxxx is one of before_cassette_lens. 

    For each possible length in before_cassette_lens, find the longest possible prefix of the 
    cassette sequence at that position in input_seq, allowing no errors. 

    If the cassette seq is found for multiple lengths, raise Exception.
    If it's not found at all, return None.
    If it's found once, return a (seq_before_cassette, seq_after_cassette, cassette_length) tuple.
    """
    cassette_seq = cassette_seq.upper()
    input_seq = input_seq.upper()
    results = []
    for before_length in before_cassette_lens:
        if_matched_bases = list(c==i for (c,i) in zip(cassette_seq, input_seq[before_length:]))
        # this would probably be faster if I encoded these as large integers and used bitwise operations
        try:                cassette_len = if_matched_bases.index(False)
        except ValueError:  cassette_len = len(if_matched_bases)
        if cassette_len >= min_len:
            results.append((before_length, cassette_len))
    if len(results) > 1:
        raise CassetteStrippingError("Cassette %s found in seq %s in multiple (pos,length) pairs! %s"%(
                                        cassette_seq, input_seq, results))
    if not len(results):
        return
    before_length, cassette_len = results[0]
    return input_seq[:before_length], input_seq[before_length+cassette_len:], cassette_len
    # TODO implement allowing mismatches - slightly complicated, need to define a trade-off between #mismatches and length
    # TODO implement allowing indels - more complicated, requires full alignment really - use Biopython?

# TODO now implement this for whole fastq files - make it a runnable utility!
# LATER-TODO integrate into preprocessing pipeline

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__strip_cassette_seq(self):
        for cassette in 'CC CCCCCCC cccccccc ccaaggggg CCATG'.split():
            for lens in ([2], [1,2,3], [2,10]):
                self.assertEqual(strip_cassette_seq('AACCTT', cassette, lens, 2), ('AA', 'TT', 2))
        self.assertEqual(strip_cassette_seq('AACCTT', 'CCT', [1,2,3], 2), ('AA', 'T', 3))
        for cassette in 'CTCT, AAA, GACTTTT'.split():
            self.assertIsNone(strip_cassette_seq('AACCTT', cassette, [2], 2))
            #self.assertRaises(CassetteStrippingError
    # LATER-TODO add unit-tests!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
