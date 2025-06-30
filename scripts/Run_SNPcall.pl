#!/usr/bin/perl

print "Usage: perl Run_SNPcall.pl <blast-dir> <snps-dir>\n" if @ARGV != 2;

use StrictUnique;

StrictUnique::SNPs(\@ARGV);
