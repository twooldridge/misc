#!/bin/bash

~/Software/tabix/bgzip -f $1
~/Software/tabix/tabix -f $1.gz
