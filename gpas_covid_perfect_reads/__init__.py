#! /usr/bin/env python3

#Use of semantic versioning, MAJOR.MINOR.MAINTAINANCE where MAJOR is not backwards compatible, but MINOR and MAINTAINANCE are
__version__="1.0.0"

from .variant_definitions import VariantGenome
from .variant_definitions import load_variant_definitions
from .variant_definitions import mutate_read
