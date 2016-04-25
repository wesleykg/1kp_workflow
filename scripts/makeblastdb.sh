#!/bin/bash

makeblastdb -in $1-SOAPdenovo-Trans-assembly.fa -dbtype nucl -parse_seqids -out $1