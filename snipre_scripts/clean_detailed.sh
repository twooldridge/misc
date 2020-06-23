#!/bin/bash

## This file gets rid of transcripts with extra-short (i.e. 3bp) coding sequences that were causing errors for SnIPRE. There were only a handful of these in each annotation.
## These transcripts will have a 3 in the Trepl column and a 0 in the Tsil column. No other transcripts should have either of those values there.

awk '{OFS="\t"} {if ($7 == 3 && $6 == 0) {next} else {print $0}}' $1 |sponge $1
